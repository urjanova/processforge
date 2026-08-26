"""``pf plan`` — validate, DOF analysis, structural diff, and Mermaid diagram."""

from __future__ import annotations

import json
import os

import typer
from loguru import logger

from ..analysis.dof import analyze_dof
from ..utils.validate_flowsheet import validate_flowsheet_dict
from ..utils.unit_consistency import check_unit_consistency, strip_units_annotations
from ..pcl import load_pcl, PCLCompileError
from ..utils.mermaid_diagram import generate_mermaid
from .common import (
    flowsheet_basename,
    load_state_manager,
    output_root,
    require_existing_file,
    validate_snapshot_config,
)
from .display import (
    print_dof_report,
    print_param_drift,
    print_provider_health,
    print_structural_diff,
    print_unit_mismatches,
)


def _find_processforge_root(start: str) -> str | None:
    """Walk up from *start*'s directory to find the dir holding `.processforge`."""
    d = os.path.dirname(os.path.abspath(start))
    while True:
        if os.path.isdir(os.path.join(d, ".processforge")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def plan(
    flowsheet: str = typer.Argument(help="Path to .pcl or .json flowsheet file"),
    output_dir: str = typer.Option(
        "diagrams",
        "--output-dir",
        "-o",
        help="Output directory for the Mermaid diagram (default: diagrams)",
    ),
    no_diagram: bool = typer.Option(
        False,
        "--no-diagram",
        help="Skip Mermaid diagram generation",
    ),
    no_dof: bool = typer.Option(
        False,
        "--no-dof",
        help="Skip Degrees-of-Freedom analysis",
    ),
    no_health: bool = typer.Option(
        False,
        "--no-health",
        help="Skip provider/container health check (faster, for CI/quick plans)",
    ),
    strict: bool = typer.Option(
        False,
        "--strict",
        help="Treat missing pip providers as hard failures (in addition to unreachable containers)",
    ),
) -> None:
    """Parse a PCL or JSON flowsheet, run validators, DOF analysis, structural diff, and emit a Mermaid diagram."""
    require_existing_file(flowsheet, label="File")

    # Step 1: Load config
    if flowsheet.endswith(".pcl"):
        try:
            json_config = load_pcl(flowsheet)
        except PCLCompileError as e:
            logger.error(f"PCL compile error: {e}")
            raise typer.Exit(code=1)
    else:
        try:
            with open(flowsheet, "r", encoding="utf-8") as f:
                json_config = json.load(f)
        except json.JSONDecodeError as e:
            logger.error(f"JSON parse error in '{flowsheet}': {e}")
            raise typer.Exit(code=1)
        except OSError as e:
            logger.error(f"Failed to read '{flowsheet}': {e}")
            raise typer.Exit(code=1)

    # Step 2: Pint unit consistency check (before _units are stripped)
    mismatches = check_unit_consistency(json_config)
    print_unit_mismatches(mismatches)

    # Step 3: Strip _units annotations before schema validation
    strip_units_annotations(json_config)

    # Step 4: Schema + connectivity validation
    try:
        config = validate_flowsheet_dict(json_config, source_name=flowsheet)
    except Exception as e:
        logger.error(f"Validation failed: {type(e).__name__}: {e}")
        raise typer.Exit(code=1)

    # Step 5: DOF analysis
    if no_dof:
        logger.info("=== Degrees of Freedom Analysis ===")
        logger.info("  Skipped (--no-dof).")
    else:
        report = analyze_dof(config)
        print_dof_report(report)

        # Step 5b: DOF fix suggestions
        if report.system_dof != 0:
            logger.warning("=== DOF Fix Suggestions ===")
            for r in report.per_unit:
                if r.dof > 0:
                    for issue in r.issues:
                        logger.warning(
                            f"  Unit '{r.unit_name}' [{r.unit_type}]: {issue}"
                        )
                elif r.dof < 0:
                    logger.warning(
                        f"  Unit '{r.unit_name}' [{r.unit_type}]: "
                        f"possibly over-specified by {abs(r.dof)} equation(s)"
                    )

    # Step 6: Structural diff vs. saved state
    base_name = flowsheet_basename(flowsheet)
    outputs_dir = output_root()
    sm, state = load_state_manager(outputs_dir, base_name)
    diff = None
    if state is not None:
        validate_snapshot_config(state, base_name)
        diff = sm.detect_structural_diff(config, state)
        print_structural_diff(diff)
        # Step 6b: Parameter drift (non-structural value changes)
        old_config = state.config
        drifted = sm.detect_drift(config, state)
        print_param_drift(drifted, old_config, config)
    else:
        logger.info("=== Structural Diff vs. Saved State ===")
        logger.info("  No prior state found — this will be a cold start.")

    # Step 7: Warm-start and homotopy eligibility
    logger.info("=== Warm-Start Status ===")
    if state is not None:
        snap_id = state.snapshot_id
        snap_ts = state.timestamp
        logger.info(f"  Warm-start available : Yes  (snapshot {snap_id}, {snap_ts})")
        topology_ok = diff is None or not bool(diff.get("added") or diff.get("removed"))
        logger.info(
            f"  Homotopy eligible    : {'Yes' if topology_ok else 'No (topology changed)'}"
        )
    else:
        logger.info("  Warm-start available : No snapshot found")
        logger.info("  Homotopy eligible    : No (cold start)")

    # Step 8: Provider / container health check
    health_failures: list[str] = []
    if no_health:
        logger.info("=== Provider / Container Health ===")
        logger.info("  Skipped (--no-health).")
    else:
        health_failures = print_provider_health(config, strict=strict)

    # Step 9: Mermaid diagram
    if not no_diagram:
        diagram_output_dir = output_dir or "diagrams"
        os.makedirs(diagram_output_dir, exist_ok=True)
        out_path = os.path.join(diagram_output_dir, f"{base_name}_plan.mmd")
        try:
            src = generate_mermaid(config)
            existing = None
            if os.path.exists(out_path):
                try:
                    with open(out_path, encoding="utf-8") as f:
                        existing = f.read()
                except OSError:
                    existing = None
            if existing == src:
                logger.info(f"Mermaid diagram unchanged: {out_path}")
            else:
                with open(out_path, "w", encoding="utf-8") as f:
                    f.write(src)
                logger.info(f"Mermaid diagram -> {out_path}")
        except Exception as e:
            logger.warning(
                f"Failed to generate Mermaid diagram: {type(e).__name__}: {e}"
            )

    # Exit non-zero on hard errors
    hard_errors = [m for m in mismatches if not m.compatible]
    if hard_errors or (not no_dof and report.system_dof < 0) or health_failures:
        raise typer.Exit(code=1)
