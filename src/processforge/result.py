import json
import os
import shutil
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import zarr
from loguru import logger

from .result_schema import ResultSchema, SolverUnitSchema, StreamSchema
from .types import RunInfo

_VARIABLE_UNITS = {"T": "K", "P": "Pa", "flowrate": "mol/s", "beta": ""}


def _ensure_parent_dir(path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)


def _ensure_array(value):
    arr = np.asarray(value)
    if arr.ndim == 0:
        arr = arr.reshape((1,))
    return arr


def _is_dynamic(results):
    for stream in results.values():
        if not hasattr(stream, "get"):
            continue
        time_series = stream.get("time")
        if isinstance(time_series, (list, tuple, np.ndarray)):
            return True
    return False


def _store_stream(stream_group, stream_data):
    for key, value in stream_data.items():
        if key == "z" and isinstance(value, dict):
            for comp in sorted(value):
                stream_group.create_array(comp, data=_ensure_array(value[comp]))
            if value:
                stream_group.attrs["composition"] = sorted(value)
            continue
        stream_group.create_array(key, data=_ensure_array(value))


def _store_solver_unit(group, data):
    for k, v in data.items():
        if isinstance(v, (np.ndarray, np.generic)):
            v = v.item()
        group.attrs[k] = v


def _store_engine_output(group, obj):
    """Persist an ``EngineOutput``/``StreamOutput`` to a Zarr group.

    Stores ``status``/``engine``/``sim_type`` as group attrs, each field as a
    unit-bearing array, and a nested ``artifacts`` subgroup whose attrs hold the
    serialized :class:`OutputArtifact` (local path for local runs, S3
    ``remote_uris`` for remote docker providers).
    """
    group.attrs["status"] = getattr(obj, "status", "")
    group.attrs["engine"] = getattr(obj, "engine", "")
    group.attrs["sim_type"] = getattr(obj, "sim_type", "")

    for f in getattr(obj, "fields", []):
        val = _ensure_array(f.quantity.value)
        arr = group.create_array(f.name, data=val)
        arr.attrs["unit"] = f.quantity.unit or ""
        if getattr(f.quantity, "std_dev", None) is not None:
            arr.attrs["std_dev"] = f.quantity.std_dev
        if getattr(f, "source", ""):
            arr.attrs["source"] = f.source

    arts = getattr(obj, "artifacts", [])
    if arts:
        ag = group.create_group("artifacts")
        for art in arts:
            ag.attrs[art.name] = json.dumps(art.model_dump())


def _normalize_run_info(run_info: RunInfo | dict) -> dict:
    """Normalize dataclass or dict run_info input to a plain dict."""
    if isinstance(run_info, RunInfo):
        return run_info.as_dict()
    return run_info


def _store_run_info(root, run_info: RunInfo | dict) -> None:
    """Persist provenance metadata in a ``run_info`` Zarr sub-group.

    Layout::

        run_info/
            .attrs  → git_hash, timestamp, python_version, platform,
                       processforge_version, mode, backend
            pkg_versions/
                .attrs  → {package: version, ...}
            initial_guess/
                x0          – float64 array  (length = n_vars)
                .attrs      → var_names  (JSON list of string labels)
    """
    run_info_dict = _normalize_run_info(run_info)
    ri_group = root.create_group("run_info")

    scalar_keys = [
        "git_hash",
        "timestamp",
        "python_version",
        "platform",
        "processforge_version",
        "mode",
        "backend",
    ]
    ri_group.attrs.update(
        {k: str(run_info_dict[k]) for k in scalar_keys if k in run_info_dict}
    )

    pkg_versions = run_info_dict.get("pkg_versions") or {}
    if pkg_versions:
        pkg_group = ri_group.create_group("pkg_versions")
        pkg_group.attrs.update({k: str(v) for k, v in pkg_versions.items()})

    x0 = run_info_dict.get("x0")
    if x0 is not None:
        x0_arr = np.asarray(x0, dtype=float)
        ig_group = ri_group.create_group("initial_guess")
        ig_group.create_array("x0", data=x0_arr)
        var_names = run_info_dict.get("var_names")
        if var_names:
            ig_group.attrs["var_names"] = list(var_names)


def _friendly_dtype(dt: str) -> str:
    dtype = str(dt)
    if dtype.startswith(("<U", "|U")) or dtype.startswith(("<S", "|S")):
        return "str"
    return dtype


def _build_schema(root) -> ResultSchema:
    from . import __version__ as _pf_version

    mode = root.attrs.get("mode", "steady")
    schema = ResultSchema(
        created=datetime.now(timezone.utc).isoformat(),
        mode=mode,
        processforge_version=_pf_version,
    )

    for name in sorted(root.group_keys()):
        if name == "run_info":
            continue
        group = root[name]
        arrays = sorted(group.array_keys())
        if arrays:
            dtypes = {k: _friendly_dtype(group[k].dtype) for k in arrays}
            units = {k: _VARIABLE_UNITS.get(k, "") for k in arrays}
            n_rows = max(group[k].shape[0] for k in arrays) if arrays else 0
            schema.streams[name] = StreamSchema(
                variables=arrays,
                dtypes=dtypes,
                units=units,
                shape=[n_rows, len(arrays)],
                has_time="time" in arrays,
                has_phase="phase" in arrays,
            )
        else:
            attrs_keys = sorted(k for k in group.attrs if not k.startswith("_"))
            if attrs_keys:
                dtypes = {k: type(group.attrs[k]).__name__ for k in attrs_keys}
                schema.solver_units[name] = SolverUnitSchema(
                    variables=attrs_keys, dtypes=dtypes
                )

    if "run_info" in root:
        ri = root["run_info"]
        prov = {}
        for key in ("git_hash", "timestamp", "mode", "backend",
                     "python_version", "platform", "processforge_version"):
            if key in ri.attrs:
                prov[key] = ri.attrs[key]
                if key == "processforge_version":
                    schema.processforge_version = ri.attrs[key]
        if "pkg_versions" in ri:
            prov["pkg_versions"] = dict(ri["pkg_versions"].attrs)
        if prov:
            schema.provenance = prov

    return schema


def _write_schema(store_path):
    store = zarr.storage.LocalStore(store_path)
    root = zarr.open_group(store=store, mode="r")
    schema = _build_schema(root)
    schema_path = store_path + ".schema.json"
    with open(schema_path, "w", encoding="utf-8") as f:
        f.write(schema.model_dump_json(indent=2))
    logger.debug(f"Wrote schema to {schema_path}")


def _write_schema_s3(s3_uri, storage_options):
    import s3fs

    store = zarr.storage.FsspecStore.from_url(s3_uri, storage_options=storage_options)
    root = zarr.open_group(store=store, mode="r")
    schema = _build_schema(root)
    fs = s3fs.S3FileSystem(**storage_options)
    schema_uri = s3_uri + ".schema.json"
    with fs.open(schema_uri, "w") as f:
        f.write(schema.model_dump_json(indent=2))
    logger.debug(f"Wrote schema to {schema_uri}")


def save_results_zarr(results, fname="results.zarr", run_info: RunInfo | dict | None = None):
    """Persist simulation results in a Zarr directory.

    Args:
        results:  Stream result dict as returned by ``Flowsheet.run()``
                  or ``EOFlowsheet.run()``.
        fname:    Output path for the Zarr store directory.
        run_info: Optional provenance metadata from
                  :func:`processforge.provenance.build_run_info`.  When
                  supplied, a ``run_info`` sub-group is written inside the
                  store containing the git hash, package versions, and
                  initial-guess vector for full reproducibility.
    """
    store_path = os.path.abspath(fname)
    if os.path.exists(store_path):
        if os.path.isdir(store_path):
            shutil.rmtree(store_path)
        else:
            os.remove(store_path)
    _ensure_parent_dir(store_path)
    store = zarr.storage.LocalStore(store_path)
    root = zarr.group(store=store)
    mode = "dynamic" if _is_dynamic(results) else "steady"
    root.attrs["mode"] = mode

    streams = {}
    solver_units = {}
    engine_outputs = {}
    for name, data in results.items():
        if hasattr(data, "fields"):
            engine_outputs[name] = data
        elif isinstance(data, dict) and "status" in data:
            solver_units[name] = data
        else:
            streams[name] = data

    for name, data in streams.items():
        group = root.create_group(name)
        _store_stream(group, data)

    for name, data in solver_units.items():
        group = root.create_group(name)
        _store_solver_unit(group, data)

    for name, data in engine_outputs.items():
        group = root.create_group(name)
        _store_engine_output(group, data)

    if run_info is not None:
        _store_run_info(root, run_info)

    _write_schema(store_path)
    logger.info(f"Saved Zarr results to {store_path}")
    return store_path


def relink_latest_results(archive_path: str, run_id: str) -> None:
    """Point ``<archive>/results.zarr`` at the latest per-run zarr.

    Each run now writes its zarr to ``<archive>/results/<run_id>/results.zarr``
    (see `save_results_zarr`). To keep backward compatibility with tooling that
    reads a single ``results.zarr`` at the archive root, we make that path a
    symlink to the newest run. Falls back to a copy on platforms/filesystems
    without symlink support (or when the existing target is a real directory).
    """
    latest = os.path.join(archive_path, "results", run_id, "results.zarr")
    if not os.path.exists(latest):
        return
    link_path = os.path.join(archive_path, "results.zarr")

    # Remove any existing target (real dir or stale/broken symlink) without
    # following the link.
    if os.path.islink(link_path):
        os.unlink(link_path)
    elif os.path.exists(link_path):
        shutil.rmtree(link_path)

    try:
        os.symlink(latest, link_path)
    except (OSError, NotImplementedError):
        # Non-POSIX / restricted filesystem: keep a physical copy instead.
        shutil.copytree(latest, link_path)


def save_results_zarr_s3(results, s3_uri: str, run_info: RunInfo | dict | None = None):
    """Write simulation results directly to an S3-backed Zarr store.

    Requires ``s3fs`` (``pip install s3fs``) and AWS credentials available via
    the standard environment variables (``AWS_ACCESS_KEY_ID``,
    ``AWS_SECRET_ACCESS_KEY``, ``AWS_DEFAULT_REGION``) or an attached IAM role.

    Args:
        results:  Stream result dict as returned by ``Flowsheet.run()``
                  or ``EOFlowsheet.run()``.
        s3_uri:   Destination URI, e.g. ``s3://my-bucket/prefix/results.zarr``.
        run_info: Optional provenance metadata (same as :func:`save_results_zarr`).
    """
    import s3fs  # noqa: F401 – registers the s3:// protocol with fsspec

    storage_options = {
        "key": os.environ.get("S3_ACCESS_KEY"),
        "secret": os.environ.get("S3_SECRET_KEY"),
    }

    client_kwargs = {}
    if endpoint := os.environ.get("S3_ENDPOINT_URL"):
        client_kwargs["endpoint_url"] = endpoint
    if region := os.environ.get("S3_REGION_NAME", "ams3"):
        client_kwargs["region_name"] = region

    if client_kwargs:
        storage_options["client_kwargs"] = client_kwargs

    store = zarr.storage.FsspecStore.from_url(s3_uri, storage_options=storage_options)
    root = zarr.group(store=store, overwrite=True)
    mode = "dynamic" if _is_dynamic(results) else "steady"
    root.attrs["mode"] = mode

    streams = {}
    solver_units = {}
    engine_outputs = {}
    for name, data in results.items():
        if hasattr(data, "fields"):
            engine_outputs[name] = data
        elif isinstance(data, dict) and "status" in data:
            solver_units[name] = data
        else:
            streams[name] = data

    for name, data in streams.items():
        group = root.create_group(name)
        _store_stream(group, data)

    for name, data in solver_units.items():
        group = root.create_group(name)
        _store_solver_unit(group, data)

    for name, data in engine_outputs.items():
        group = root.create_group(name)
        _store_engine_output(group, data)

    if run_info is not None:
        _store_run_info(root, run_info)

    _write_schema_s3(s3_uri, storage_options)
    logger.info(f"Saved Zarr results to {s3_uri}")


def plot_results(results, fname="results.png"):
    os.makedirs("outputs", exist_ok=True)
    streams = list(results.keys())
    temp_values = [_scalar_from_sequence(results[s].get("T")) or 0.0 for s in streams]

    plt.figure(figsize=(8, 4))
    plt.bar(streams, temp_values, color="skyblue")
    plt.ylabel("Temperature [K]")
    plt.title("Stream Temperatures (Steady State)")
    plt.tight_layout()
    plt.savefig(os.path.join("outputs", "temps_" + fname))
    plt.close()

    comps = set()
    for s in results.values():
        comps.update(s.get("z", {}).keys())
    comps = sorted(comps)

    bottom = [0.0] * len(streams)
    plt.figure(figsize=(8, 4))
    for comp in comps:
        vals = [
            _scalar_from_sequence(results[s].get("z", {}).get(comp, 0.0)) or 0.0
            for s in streams
        ]
        plt.bar(streams, vals, bottom=bottom, label=comp)
        bottom = [bottom[i] + vals[i] for i in range(len(vals))]
    plt.ylabel("Mole Fraction")
    plt.title("Stream Compositions (Steady State)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join("outputs", "comps_" + fname))
    plt.close()


def plot_timeseries(results, fname="timeseries.png"):
    os.makedirs("outputs", exist_ok=True)
    streams = sorted(results.keys())
    if not streams:
        return

    times = results[streams[0]].get("time", [])
    if not times:
        return

    plt.figure(figsize=(8, 5))
    for s_name in streams:
        if "T" in results[s_name]:
            plt.plot(times, results[s_name]["T"], label=s_name)
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [K]")
    plt.title("Stream Temperatures vs Time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join("outputs", "temps_" + fname))
    plt.close()

    comps = set()
    for s_data in results.values():
        if "z" in s_data and isinstance(s_data["z"], dict):
            comps.update(s_data["z"].keys())
    comps = sorted(comps)

    for s_name in streams:
        if "z" not in results[s_name]:
            continue
        plt.figure(figsize=(8, 5))
        for comp in comps:
            if comp in results[s_name]["z"]:
                plt.plot(times, results[s_name]["z"][comp], label=comp)
        plt.xlabel("Time [s]")
        plt.ylabel("Mole Fraction")
        plt.title(f"Compositions vs Time ({s_name})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join("outputs", f"comps_{s_name}_" + fname))
        plt.close()


def _convert_value(value):
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return ""
        return _convert_value(value.item())
    if isinstance(value, np.generic):
        return value.item()
    return value


def _scalar_from_sequence(value):
    if value is None:
        return None
    if isinstance(value, (list, tuple, np.ndarray)):
        arr = np.asarray(value)
        if arr.size == 0:
            return None
        candidate = arr[-1]
    else:
        candidate = value
    return _convert_value(candidate)


def _load_dataframe_from_zarr(store_path):
    store = zarr.storage.LocalStore(store_path)
    root = zarr.open(store=store, mode="r")
    rows = []

    for stream in sorted(k for k in root.group_keys() if k != "run_info"):
        group = root[stream]
        arrays = list(group.array_keys())
        if not arrays:
            continue

        has_time = "time" in group
        time_arr = group["time"][:] if has_time else None
        has_phase = "phase" in group
        phase_arr = group["phase"][:] if has_phase else None

        data_vars = [k for k in arrays if k not in ("time", "phase")]
        if not data_vars:
            continue
        n = max(group[k].shape[0] for k in data_vars)

        for idx in range(n):
            row = {"stream": stream}
            if has_time:
                row["time"] = float(time_arr[idx])
            for var in data_vars:
                row[var] = _convert_value(group[var][idx])
            if has_phase:
                row["phase"] = str(phase_arr[idx])
            rows.append(row)

    return pd.DataFrame(rows)


def summarize_zarr_store(store_path: str) -> dict:
    """Return a structured summary of a Zarr result store.

    The returned dict contains enough provider-agnostic information for the
    ``pf runs`` command to render a one-line summary or a detailed result view.
    """
    summary: dict = {"present": False}
    if not os.path.isdir(store_path):
        return summary

    store = zarr.storage.LocalStore(store_path)
    root = zarr.open_group(store=store, mode="r")

    schema = None
    schema_path = store_path + ".schema.json"
    if os.path.isfile(schema_path):
        with open(schema_path, "r", encoding="utf-8") as f:
            schema = json.load(f)

    summary = {
        "present": True,
        "mode": root.attrs.get("mode", "unknown"),
        "schema": schema,
        "streams": {},
        "solver_units": {},
        "engine_outputs": {},
        "artifacts": [],
    }

    for name in root.group_keys():
        if name == "run_info":
            continue
        group = root[name]
        arrays = sorted(group.array_keys())
        attrs = dict(group.attrs)

        if arrays:
            fields = {}
            for arr_name in arrays:
                arr = group[arr_name]
                fields[arr_name] = {
                    "value": _convert_value(arr[:]),
                    "unit": arr.attrs.get("unit", ""),
                    "std_dev": arr.attrs.get("std_dev", None),
                }
            if attrs.get("engine"):
                summary["engine_outputs"][name] = {
                    "engine": attrs.get("engine", ""),
                    "sim_type": attrs.get("sim_type", ""),
                    "status": attrs.get("status", ""),
                    "fields": fields,
                }
            else:
                summary["streams"][name] = {
                    "variables": arrays,
                    "units": {k: group[k].attrs.get("unit", "") for k in arrays},
                    "fields": fields,
                }
        else:
            if attrs.get("engine"):
                summary["engine_outputs"][name] = {
                    "engine": attrs.get("engine", ""),
                    "sim_type": attrs.get("sim_type", ""),
                    "status": attrs.get("status", ""),
                    "fields": {},
                    "attrs": attrs,
                }
            else:
                summary["solver_units"][name] = {
                    "variables": sorted(k for k in attrs if not k.startswith("_")),
                    "attrs": attrs,
                }

        # Extract artifacts from an engine-output or solver-unit subgroup.
        if "artifacts" in group.group_keys():
            art_group = group["artifacts"]
            for art_name, art_json in art_group.attrs.items():
                try:
                    summary["artifacts"].append(json.loads(art_json))
                except json.JSONDecodeError:
                    continue

    return summary


def _one_line_summary(summary: dict) -> str:
    """Pick a short, representative result string from a Zarr summary."""
    if not summary.get("present"):
        return "no results.zarr"

    for name, eo in summary.get("engine_outputs", {}).items():
        engine = eo.get("engine", "").lower()
        fields = eo.get("fields", {})

        if engine == "openmc":
            k_eff = fields.get("k_eff")
            if k_eff is not None:
                val = k_eff["value"]
                std = k_eff.get("std_dev")
                if std is not None:
                    return f"k_eff={val:.5f}±{std:.5f}"
                return f"k_eff={val:.5f}"
            return "openmc completed"

        if engine == "festim":
            for fname, fdata in fields.items():
                if "flux" in fname.lower():
                    unit = fdata.get("unit", "")
                    return f"{fname}={_fmt_scientific(fdata['value'])} {unit}".strip()
            if fields:
                first = next(iter(fields.values()))
                unit = first.get("unit", "")
                return f"{next(iter(fields))}={_fmt_scientific(first['value'])} {unit}".strip()
            return "festim completed"

    n_streams = len(summary.get("streams", {}))
    if n_streams:
        return f"{n_streams} streams"

    return "completed"


def _fmt_scientific(value) -> str:
    """Format a scalar numerically; fall back to its string representation."""
    if value is None:
        return ""
    try:
        f = float(value)
        if abs(f) >= 1e4 or (abs(f) < 1e-3 and abs(f) > 0):
            return f"{f:.3e}"
        return f"{f:.4f}"
    except (TypeError, ValueError):
        return str(value)

