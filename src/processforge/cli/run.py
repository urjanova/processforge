"""``pf run`` — run a process simulation from a flowsheet JSON file."""

from __future__ import annotations

import os

import typer
from loguru import logger

from ..flowsheet import Flowsheet
from ..eo import EOFlowsheet
from ..provenance import build_dynamic_x0, build_run_info
from ..result import plot_results, plot_timeseries, save_results_zarr
from ..state import StateManager
from .common import (
    check_providers,
    output_root,
    require_existing_file,
    validate_runtime_flowsheet,
)


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

    # Check provider availability
    check_providers(config, flowsheet)

    base_name = os.path.splitext(os.path.basename(flowsheet))[0]
    outputs_dir = output_root()

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    is_dynamic = mode == "dynamic"

    if is_dynamic:
        # Load .pfstate as t=0 if available
        state_path = os.path.join(outputs_dir, f"{base_name}.pfstate")
        sm = StateManager(state_path)
        state = sm.load_state()
        if state is not None:
            stream_inits = sm.state_to_stream_dicts(state)
            streams_cfg = config.get("streams", {})
            for s_name, s_vals in stream_inits.items():
                if s_name in streams_cfg:
                    streams_cfg[s_name].update(s_vals)
                else:
                    logger.debug(
                        f"Stream '{s_name}' from snapshot not in current flowsheet — skipped."
                    )
            logger.info("Using .pfstate converged state as dynamic t=0.")
        else:
            logger.debug("No .pfstate found — starting dynamic run from flowsheet defaults.")

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

    os.makedirs(outputs_dir, exist_ok=True)
    zarr_path = os.path.join(outputs_dir, f"{base_name}_results.zarr")
    save_results_zarr(results, zarr_path, run_info=run_info)
    logger.info(f"Results saved to {zarr_path}")

    if export_images:
        try:
            plot_results(results, fname=f"{base_name}_results.png")
            plot_timeseries(results, fname=f"{base_name}_timeseries.png")
            logger.info(f"Plots saved: {base_name}_results.png, {base_name}_timeseries.png")
        except Exception as e:
            logger.warning(f"Failed to generate plots: {type(e).__name__}: {e}")
