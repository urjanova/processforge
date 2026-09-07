"""``pf run`` — run a process simulation from a flowsheet JSON file."""
from __future__ import annotations

import os
from datetime import datetime, timezone

import typer
from loguru import logger

from ..eo import EOFlowsheet
from ..flowsheet import Flowsheet
from ..persistence.archive import ProcessStateArchive
from ..provenance import build_dynamic_x0, build_run_info
from .common import (
    check_providers,
    display_backend,
    flowsheet_basename,
    output_root,
    require_existing_file,
    validate_runtime_flowsheet,
    _resolve_provider_url,
)
from .persist import flowsheet_hash, make_run_id, persist_run


def _log_active_providers(config: dict) -> None:
    """Log a one-line summary of the provider(s) a run will dispatch to.

    Highlights containerized providers (which do their real work inside a Docker
    container and can take a while) so the user sees what is running instead of
    a silent wait. Pip-installable providers are mentioned for completeness.
    """
    from ..providers.registry import is_containerized

    providers = config.get("providers", {})
    if not providers:
        return

    containers = []
    pips = []
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = _resolve_provider_url(cfg, ptype)
            containers.append(f"{name} [{ptype}] @ {url}")
        else:
            pips.append(f"{name} [{ptype}]")

    if containers:
        logger.info("=== Active providers (containerized) ===")
        for entry in containers:
            logger.info(f"  → {entry}")
    if pips:
        logger.info("=== Active providers (local) ===")
        for entry in pips:
            logger.info(f"  → {entry}")


def run(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
    export_images: bool = typer.Option(
        False,
        "--export-images",
        help="Generate PNG plots for simulation outputs",
    ),
) -> None:
    """Run a process simulation from a flowsheet JSON file."""
    require_existing_file(flowsheet)
    config = validate_runtime_flowsheet(flowsheet)

    # Check provider availability (assumes any containers are already running)
    check_providers(config, flowsheet)

    # Surface which containerized provider(s) this run will dispatch to, so the
    # command line isn't silent while the container does the heavy compute.
    _log_active_providers(config)

    base_name = flowsheet_basename(flowsheet)
    outputs_dir = output_root()

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    is_dynamic = mode == "dynamic"

    # Standardized run id, shared with container-side artifact uploads.
    run_id = make_run_id()
    flowsheet_hash_value = flowsheet_hash(config)

    if is_dynamic:
        # Load .pfarchive snapshot as t=0 if available.
        archive_path = os.path.join(outputs_dir, f"{base_name}.pfarchive")
        archive = ProcessStateArchive(archive_path)
        state = archive.load_snapshot()
        if state is not None:
            stream_inits = archive.state_to_stream_dicts(state)
            streams_cfg = config.get("streams", {})
            for s_name, s_vals in stream_inits.items():
                if s_name in streams_cfg:
                    streams_cfg[s_name].update(s_vals)
                else:
                    logger.debug(
                        f"Stream '{s_name}' from snapshot not in current flowsheet — skipped."
                    )
            logger.info("Using .pfarchive converged state as dynamic t=0.")
        else:
            logger.debug("No .pfarchive found — starting dynamic run from flowsheet defaults.")

        fs = Flowsheet(config)
        fs.set_run_context(run_id, flowsheet_hash_value)
        logger.info("=== Dynamic Results ===")
        results = fs.run()

        _check_for_failed_units(fs)

        if hasattr(fs, "converged"):
            if fs.converged:
                logger.info("Dynamic simulation converged.")
            else:
                logger.warning("Dynamic simulation did NOT converge. Results may be unreliable.")

        try:
            x0, var_names = build_dynamic_x0(config)
        except Exception as e:
            logger.error(f"Failed to build initial state vector: {type(e).__name__}: {e}")
            raise SystemExit(1)

        run_info = build_run_info(config, x0=x0, var_names=var_names)
        snap_x, snap_vn, snap_backend, snap_ok = x0, var_names, "dynamic", True
    else:
        # Pass None so EOFlowsheet resolves backend from config (with scipy default).
        fs = EOFlowsheet(config, backend=None)
        fs.set_run_context(run_id, flowsheet_hash_value)
        logger.info("=== Steady-State EO Results ===")
        results = fs.run()

        _check_for_failed_units(fs)

        if hasattr(fs, "converged"):
            if fs.converged:
                logger.info("Steady-state simulation converged.")
            else:
                logger.warning("Steady-state simulation did NOT converge. Results may be unreliable.")

        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
        snap_x, snap_vn = fs.x_converged, fs.var_names
        snap_backend = getattr(fs, "backend", "scipy")
        snap_ok = bool(getattr(fs, "converged", False))

    # Build + persist the standardized run manifest.
    archive_path = os.path.join(outputs_dir, f"{base_name}.pfarchive")
    archive = ProcessStateArchive(archive_path)
    persist_run(archive, fs, run_id, results, run_info, config, base_name)

    # Persist a StateManager snapshot so `pf plan` (and `pf apply`) can diff
    # against this run as a baseline. Without it, `pf plan` always reports
    # "No prior state found" and never surfaces flowsheet edits.
    if (is_dynamic or snap_ok) and snap_x is not None:
        try:
            from .common import build_run_metadata, save_snapshot

            meta = build_run_metadata(config, 1e-6, 50, snap_backend)
            save_snapshot(
                archive, config, snap_x, snap_vn,
                metadata=meta, parent_snapshot_id=None, label="run state",
            )
        except Exception as e:
            logger.warning(f"Failed to save state snapshot: {type(e).__name__}: {e}")

    logger.info(f"Backend      : {display_backend(config, getattr(fs, 'backend', 'dynamic'))}")
    logger.info(f"Run saved    : {os.path.join(archive_path, 'runs', run_id + '.json')}")

    # Summarize standardized outputs.
    manifest = archive.load_run(run_id)
    for unit_name, out in manifest.units.items():
        for f in out.fields:
            logger.info(f"  [{unit_name}] {f.name} = {f.quantity.value} {f.quantity.unit}")
        for art in out.artifacts:
            if art.remote_uris:
                logger.info(f"  [{unit_name}] artifact {art.name} → {art.remote_uris[0]}")
    for stream_name, out in manifest.streams.items():
        for f in out.fields:
            logger.info(f"  [{stream_name}] {f.name} = {f.quantity.value} {f.quantity.unit}")

    if export_images:
        try:
            from ..result import plot_results, plot_timeseries

            plot_results(stream_results, fname=f"{base_name}_results.png")
            plot_timeseries(stream_results, fname=f"{base_name}_timeseries.png")
            logger.info(f"Plots saved: {base_name}_results.png, {base_name}_timeseries.png")
        except Exception as e:
            logger.warning(f"Failed to generate plots: {type(e).__name__}: {e}")


def _check_for_failed_units(fs):
    """Fail loudly if any SolverUnit run returned ``status="failed"``.

    Containerized providers return HTTP 200 with a structured ``EngineOutput``
    (status="failed") rather than raising, so the failure would otherwise be
    swallowed and recorded as a successful run. Surface it clearly and exit
    non-zero.
    """
    from ..types import EngineOutput

    failed = {
        name: out
        for name, out in getattr(fs, "engine_outputs", {}).items()
        if isinstance(out, EngineOutput) and out.status == "failed"
    }
    if not failed:
        return

    for name, out in failed.items():
        err = out.error
        category = getattr(err, "category", "unknown") if err else "unknown"
        message = getattr(err, "message", "") if err else ""
        logger.error(f"Unit '{name}' simulation FAILED [{category}]: {message}")
        hint = getattr(err, "hint", "") if err else ""
        if hint:
            logger.error(f"  hint: {hint}")
    raise SystemExit(1)
