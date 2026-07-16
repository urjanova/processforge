"""``pf apply`` — state-based warm start with homotopy fallback."""

from __future__ import annotations

import os
import time
from typing import Literal

import typer
from loguru import logger

from ..eo import EOFlowsheet
from ..provenance import build_run_info
from ..result import save_results_zarr
from ..state import StateManager
from .common import (
    build_divergence_report,
    build_run_metadata,
    check_providers,
    load_state_manager,
    log_residual_breakdown,
    output_root,
    require_existing_file,
    save_snapshot,
    validate_runtime_flowsheet,
    write_divergence_report,
)
from .display import print_structural_diff


def apply(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
    backend: Literal["scipy", "pyomo", "casadi"] | None = typer.Option(
        None,
        "--backend",
        help="Override the flowsheet's simulation.backend",
    ),
    tolerance: float = typer.Option(
        1e-6,
        "--tolerance",
        help="Newton solver convergence tolerance (default: 1e-6)",
    ),
    max_iter: int = typer.Option(
        50,
        "--max-iter",
        help="Max Newton iterations (default: 50)",
    ),
    skip_homotopy: bool = typer.Option(
        False,
        "--skip-homotopy",
        help="Disable homotopy fallback; cold-start only",
    ),
) -> None:
    """Apply flowsheet: drift detection, warm-start, homotopy fallback, convergence guardrails."""
    require_existing_file(flowsheet)
    config = validate_runtime_flowsheet(flowsheet)

    # Check provider availability
    check_providers(config, flowsheet)

    base_name = os.path.splitext(os.path.basename(flowsheet))[0]

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    if mode != "steady":
        logger.error("pf apply is only supported for steady-state EO flowsheets.")
        raise SystemExit(1)

    outputs_dir = output_root()
    os.makedirs(outputs_dir, exist_ok=True)
    sm, state = load_state_manager(outputs_dir, base_name)

    # Structural diff: detect topology changes
    topology_changed = False
    if state is not None:
        diff = sm.detect_structural_diff(config, state)
        print_structural_diff(diff)
        topology_changed = bool(diff.get("added") or diff.get("removed"))
        if topology_changed:
            logger.warning(
                "Topology changed (units added/removed). "
                "Homotopy requires identical topology — falling back to cold start."
            )

    # Parameter drift (only meaningful when topology is unchanged)
    drifted: list[str] = []
    current_metadata = build_run_metadata(config, tolerance, max_iter, backend or "scipy")
    if state is not None and not topology_changed:
        mismatches = sm.validate_metadata(current_metadata, state)
        if mismatches:
            logger.warning(f"Metadata mismatch: {mismatches}")
        drifted = sm.detect_drift(config, state)
        if not drifted:
            logger.info("No drift detected. System is already at the desired state.")
            return
        stream_drifts = [d for d in drifted if d.startswith("streams.")]
        unit_drifts = [d for d in drifted if d.startswith("units.")]
        logger.warning("Drift detected:")
        if stream_drifts:
            logger.warning(f"  Stream drift : {stream_drifts}")
        if unit_drifts:
            logger.warning(f"  Unit drift   : {unit_drifts}")

    # Build flowsheet; attach saved state for warm-start unless topology changed
    fs = EOFlowsheet(config, backend=backend)
    fs.saved_state = state if not topology_changed else None
    fs.solver_tol = tolerance
    fs.solver_max_iter = max_iter

    logger.info("=== Running Apply (Steady-State EO) ===")
    t0 = time.perf_counter()
    results = fs.run()
    elapsed = time.perf_counter() - t0

    logger.info(
        f"Direct solve completed in {elapsed:.2f}s (converged={fs.converged})."
    )

    if fs.converged:
        snapshot_id = save_snapshot(
            sm, config, fs.x_converged, fs.var_names,
            metadata=current_metadata,
            parent_snapshot_id=state.snapshot_id if state is not None and not topology_changed else None,
            label="converged state",
        )
        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
        zarr_path = os.path.join(outputs_dir, f"{base_name}_results.zarr")
        save_results_zarr(results, zarr_path, run_info=run_info)
        logger.info("=== Apply Summary ===")
        logger.info("  Status       : CONVERGED")
        logger.info(f"  Final ||F||  : {fs.solver_stats.get('final_norm', '?'):.3e}")
        logger.info(f"  Iterations   : {fs.solver_stats.get('iterations', '?')}")
        logger.info(f"  Backend      : {fs.backend}")
        logger.info(f"  Snapshot ID  : {snapshot_id}")
        logger.info(f"  Results zarr : {zarr_path}")
        logger.info(f"  Elapsed (s)  : {elapsed:.2f}")
        return

    # Direct solve failed — try homotopy (only when topology is same and state exists)
    if state is not None and not topology_changed and drifted and not skip_homotopy:
        logger.warning("Direct solve failed. Attempting homotopy continuation...")
        from ..eo.solver import EOSolver, solve_with_homotopy
        from ..eo.flowsheet import EOFlowsheet as _EO
        from ..providers.manager import teardown_providers

        solver = EOSolver(backend=fs.backend, tol=tolerance, max_iter=max_iter)
        tmp_fs = _EO(config, backend=backend)
        manager = tmp_fs._build()
        try:
            x_hom, hom_converged, hom_stats = solve_with_homotopy(
                tmp_fs, manager, solver, state, drifted
            )
        finally:
            teardown_providers(tmp_fs._provider_map)

        if hom_converged:
            logger.info(
                f"Homotopy converged: ||F||={hom_stats.get('final_norm', '?'):.3e}, "
                f"iterations={hom_stats.get('iterations', '?')}"
            )
            save_snapshot(
                sm, config, x_hom, fs.var_names,
                metadata=current_metadata,
                parent_snapshot_id=state.snapshot_id if state is not None and not topology_changed else None,
                label="homotopy solution",
            )
            run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
            save_results_zarr(
                results,
                os.path.join(outputs_dir, f"{base_name}_results.zarr"),
                run_info=run_info,
            )
            logger.info("Homotopy apply succeeded. New snapshot saved.")
            return

        # Both failed — auto-revert and write divergence report
        logger.error("Homotopy also failed to converge.")
        prev_id = state.snapshot_id if state is not None else "unknown"
        sm.rollback(1)
        logger.warning(f"Reverted .pfstate to snapshot before {prev_id}.")
        breakdown = log_residual_breakdown(fs)
        divergence = build_divergence_report(
            drifted_params=drifted,
            solver_stats=hom_stats,
            x_last=x_hom,
            var_names=fs.var_names,
            breakdown=breakdown,
        )
    else:
        # Cold start also failed; no rollback (nothing to revert to)
        if not drifted and state is None:
            logger.error("Cold-start solve failed to converge (no prior snapshot).")
        elif skip_homotopy:
            logger.error("Cold-start solve failed to converge (--skip-homotopy).")
        else:
            logger.error("Cold-start solve failed to converge.")
        breakdown = log_residual_breakdown(fs)
        divergence = build_divergence_report(
            drifted_params=[],
            solver_stats=fs.solver_stats,
            x_last=getattr(fs, "x_converged", []),
            var_names=getattr(fs, "var_names", []),
            breakdown=breakdown,
        )

    write_divergence_report(outputs_dir, base_name, divergence)
    raise SystemExit(1)
