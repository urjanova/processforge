"""``pf plan`` — validate, DOF analysis, structural diff, and Mermaid diagram."""

from __future__ import annotations

import argparse
import json
import os

from loguru import logger

from ..utils.validate_flowsheet import validate_flowsheet_dict
from ..state import StateManager
from ..utils.mermaid_diagram import generate_mermaid
from .common import (
    load_state_manager,
    output_root,
    require_existing_file,
    validate_snapshot_config,
)
from .display import print_dof_report, print_structural_diff, print_unit_mismatches


def add_plan_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the ``plan`` subcommand."""
    parser.add_argument("flowsheet", help="Path to .pcl or .json flowsheet file")
    parser.add_argument(
        "--output-dir", "-o", default="diagrams",
        help="Output directory for the Mermaid diagram (default: diagrams)",
    )
    parser.add_argument(
        "--no-diagram", action="store_true",
        help="Skip Mermaid diagram generation",
    )


def cmd_plan(args: argparse.Namespace) -> None:
    """Parse a PCL or JSON flowsheet, run validators, DOF analysis, structural diff, and emit a Mermaid diagram."""
    fname = args.flowsheet
    require_existing_file(fname, label="File")

    # Step 1: Load config
    if fname.endswith(".pcl"):
        from ..pcl import load_pcl, PCLCompileError
        try:
            json_config = load_pcl(fname)
        except PCLCompileError as e:
            logger.error(f"PCL compile error: {e}")
            raise SystemExit(1)
    else:
        try:
            with open(fname, "r", encoding="utf-8") as f:
                json_config = json.load(f)
        except json.JSONDecodeError as e:
            logger.error(f"JSON parse error in '{fname}': {e}")
            raise SystemExit(1)

    # Step 2: Pint unit consistency check (before _units are stripped)
    from ..utils.unit_consistency import check_unit_consistency, strip_units_annotations
    mismatches = check_unit_consistency(json_config)
    print_unit_mismatches(mismatches)

    # Step 3: Strip _units annotations before schema validation
    strip_units_annotations(json_config)

    # Step 4: Schema + connectivity validation
    try:
        config = validate_flowsheet_dict(json_config, source_name=fname)
    except SystemExit:
        raise
    except Exception as e:
        logger.error(f"Validation failed unexpectedly: {type(e).__name__}: {e}")
        raise SystemExit(1)

    # Step 5: DOF analysis
    from ..analysis.dof import analyze_dof
    report = analyze_dof(config)
    print_dof_report(report)

    # Step 5b: DOF fix suggestions
    if report.system_dof != 0:
        logger.warning("=== DOF Fix Suggestions ===")
        for r in report.per_unit:
            if r.dof > 0:
                for issue in r.issues:
                    logger.warning(f"  Unit '{r.unit_name}' [{r.unit_type}]: {issue}")
            elif r.dof < 0:
                logger.warning(
                    f"  Unit '{r.unit_name}' [{r.unit_type}]: "
                    f"possibly over-specified by {abs(r.dof)} equation(s)"
                )

    # Step 6: Structural diff vs. saved state
    base_name = os.path.splitext(os.path.basename(fname))[0]
    outputs_dir = output_root()
    sm, state = load_state_manager(outputs_dir, base_name)
    diff = None
    if state is not None:
        validate_snapshot_config(state, base_name)
        diff = sm.detect_structural_diff(config, state)
        print_structural_diff(diff)
    else:
        logger.info("=== Structural Diff vs. Saved State ===")
        logger.info("  No prior state found — this will be a cold start.")

    # Step 6b: Warm-start and homotopy eligibility
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

    # Step 7: Mermaid diagram
    if not args.no_diagram:
        output_dir = args.output_dir or "diagrams"
        os.makedirs(output_dir, exist_ok=True)
        out_path = os.path.join(output_dir, f"{base_name}_plan.mmd")
        try:
            generate_mermaid(config, output_path=out_path)
            logger.info(f"Mermaid diagram -> {out_path}")
        except Exception as e:
            logger.warning(f"Failed to generate Mermaid diagram: {type(e).__name__}: {e}")

    # Exit non-zero on hard errors
    hard_errors = [m for m in mismatches if not m.compatible]
    if hard_errors or report.system_dof < 0:
        raise SystemExit(1)
