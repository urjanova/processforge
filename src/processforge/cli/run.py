"""``pf run`` — run a process simulation from a flowsheet JSON file."""
from __future__ import annotations

import hashlib
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
)


def _flowsheet_hash(config: dict) -> str:
    return hashlib.sha256(
        json_dumps(config).encode("utf-8")
    ).hexdigest()[:16]


def json_dumps(obj) -> str:
    import json

    return json.dumps(obj, sort_keys=True, default=str)


def _run_id() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + os.urandom(3).hex()


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

    base_name = flowsheet_basename(flowsheet)
    outputs_dir = output_root()

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    is_dynamic = mode == "dynamic"

    # Standardized run id, shared with container-side artifact uploads.
    run_id = _run_id()
    flowsheet_hash = _flowsheet_hash(config)

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
        logger.info("=== Dynamic Results ===")
        results = fs.run()

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
    else:
        # Pass None so EOFlowsheet resolves backend from config (with scipy default).
        fs = EOFlowsheet(config, backend=None)
        logger.info("=== Steady-State EO Results ===")
        results = fs.run()

        if hasattr(fs, "converged"):
            if fs.converged:
                logger.info("Steady-state simulation converged.")
            else:
                logger.warning("Steady-state simulation did NOT converge. Results may be unreliable.")

        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)

    # Propagate run context to any containerized providers so their S3 uploads
    # are keyed by the same run_id/flowsheet_hash as this archive run.
    for provider in getattr(fs, "_provider_map", {}).values():
        if hasattr(provider, "set_run_context"):
            provider.set_run_context(run_id, flowsheet_hash)

    # Build + persist the standardized run manifest.
    archive_path = os.path.join(outputs_dir, f"{base_name}.pfarchive")
    archive = ProcessStateArchive(archive_path)
    manifest = fs.collect_outputs(run_id, mode, base_name, provenance=run_info.model_dump())
    stream_results = {
        k: v for k, v in results.items()
        if not hasattr(v, "fields")  # exclude EngineOutput objects
    }
    archive.save_run(manifest, stream_results=stream_results)

    logger.info(f"Backend      : {display_backend(config, getattr(fs, 'backend', 'dynamic'))}")
    logger.info(f"Run saved    : {os.path.join(archive_path, 'runs', run_id + '.json')}")

    # Summarize standardized outputs.
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
