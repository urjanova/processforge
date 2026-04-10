import os
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import zarr
from loguru import logger

from .types import RunInfo


def _ensure_parent_dir(path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)


def _ensure_array(value):
    arr = np.asarray(value)
    if arr.ndim == 0:
        arr = arr.reshape((1,))
    return arr


def _is_dynamic(results):
    for stream in results.values():
        time_series = stream.get("time")
        if isinstance(time_series, (list, tuple, np.ndarray)):
            return True
    return False


def _store_stream(stream_group, stream_data):
    for key, value in stream_data.items():
        if key == "z" and isinstance(value, dict):
            comp_group = stream_group.create_group("__composition__")
            for comp, comp_values in sorted(value.items()):
                comp_arr = _ensure_array(comp_values)
                comp_group.create_dataset(comp, data=comp_arr, shape=comp_arr.shape)
            continue
        arr = _ensure_array(value)
        stream_group.create_dataset(key, data=arr, shape=arr.shape)


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
        ig_group.create_dataset("x0", data=x0_arr, shape=x0_arr.shape)
        var_names = run_info_dict.get("var_names")
        if var_names:
            ig_group.attrs["var_names"] = list(var_names)


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
    for stream_name, stream_data in results.items():
        stream_group = root.create_group(stream_name)
        _store_stream(stream_group, stream_data)
    if run_info is not None:
        _store_run_info(root, run_info)
    logger.info(f"Saved Zarr results to {store_path}")
    return store_path


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


def _get_zarr_value(dataset, idx):
    if dataset is None:
        return ""
    length = dataset.shape[0]
    if length == 0:
        return ""
    position = min(idx, length - 1)
    return _convert_value(dataset[position])


def _build_dataframe_row(group, stream_name, idx, comp_names, include_time):
    row = {"stream": stream_name}
    if include_time and "time" in group:
        row["time"] = _get_zarr_value(group.get("time"), idx)
    row["T [K]"] = _get_zarr_value(group.get("T"), idx)
    row["P [Pa]"] = _get_zarr_value(group.get("P"), idx)
    row["Phase"] = _get_zarr_value(group.get("phase"), idx)
    row["VaporFrac"] = _get_zarr_value(group.get("beta"), idx)
    row["flowrate"] = _get_zarr_value(group.get("flowrate"), idx)
    comp_group = group.get("__composition__")
    if comp_group is not None:
        for comp in comp_names:
            row[comp] = _get_zarr_value(comp_group.get(comp), idx)
    return row


def _load_dataframe_from_zarr(store_path):
    store = zarr.storage.LocalStore(store_path)
    root = zarr.open(store=store, mode="r")
    streams = sorted(k for k in root.group_keys() if k != "run_info")
    rows = []
    components = set()
    mode = root.attrs.get("mode", "steady")
    for stream in streams:
        group = root[stream]
        comp_group = group.get("__composition__")
        comp_names = sorted(comp_group.keys()) if comp_group is not None else []
        components.update(comp_names)
        has_time = "time" in group and mode == "dynamic"
        length = group["time"].shape[0] if has_time else 1
        for idx in range(length):
            rows.append(_build_dataframe_row(group, stream, idx, comp_names, has_time))
    df = pd.DataFrame(rows)
    for comp in sorted(components):
        if comp in df:
            df[comp] = df[comp].fillna(0.0)
        else:
            df[comp] = 0.0
    return df

