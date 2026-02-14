import argparse
import json
import os

from loguru import logger

from .utils.validate_flowsheet import validate_flowsheet
from .flowsheet import Flowsheet
from .result import (
    save_results_csv,
    save_timeseries_csv,
    save_results_json,
    save_timeseries_json,
    generate_validation_excel,
)
from .utils.flowsheet_diagram import draw_flowsheet


def _cmd_run(args):
    """Run a process simulation from a flowsheet JSON file."""
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    try:
        config = validate_flowsheet(fname)
    except SystemExit:
        raise
    except Exception as e:
        logger.error(f"Failed to validate flowsheet file '{fname}': {e}")
        raise SystemExit(1)

    fs = Flowsheet(config)
    results = fs.run()

    is_dynamic = config.get("simulation", {}).get("mode") == "dynamic"
    base_name = os.path.splitext(os.path.basename(fname))[0]

    if is_dynamic:
        logger.info("=== Dynamic Results ===")
        save_timeseries_json(results, f"{base_name}_timeseries.json")
        save_results_json(results, f"{base_name}_results.json")
        save_timeseries_csv(results, f"{base_name}_timeseries.csv")
        generate_validation_excel(
            data_source=os.path.join("outputs", f"{base_name}_timeseries.csv"),
            output_filename=os.path.join("outputs", f"{base_name}_validation.xlsx"),
        )
    else:
        logger.info("=== Steady-State Results ===")
        save_results_json(results, f"{base_name}_results.json")
        save_results_csv(results, f"{base_name}_results.csv")
        save_timeseries_json(results, f"{base_name}_timeseries.json")
        generate_validation_excel(
            data_source=os.path.join("outputs", f"{base_name}_results.csv"),
            output_filename=os.path.join("outputs", f"{base_name}_validation.xlsx"),
        )
        logger.info(f"Saved {base_name}_results.csv and {base_name}_validation.xlsx")


def _cmd_validate(args):
    """Validate a flowsheet JSON file."""
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    try:
        validate_flowsheet(fname)
        logger.info(f"Flowsheet '{fname}' is valid.")
    except SystemExit:
        raise
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        raise SystemExit(1)


def _cmd_diagram(args):
    """Generate a flowsheet diagram from a JSON file."""
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    with open(fname, "r") as f:
        flowsheet_schema = json.load(f)

    output_dir = args.output_dir or "."
    fmt = args.format or "png"

    base_name = os.path.splitext(os.path.basename(fname))[0]

    output_path = draw_flowsheet(
        flowsheet_schema,
        output_directory=output_dir,
        output_filename=base_name,
        file_format=fmt,
    )
    logger.info(f"Diagram saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="ProcessForge - Chemical Process Simulation",
        prog="processforge",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # processforge run
    run_parser = subparsers.add_parser("run", help="Run a process simulation")
    run_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")

    # processforge validate
    validate_parser = subparsers.add_parser("validate", help="Validate a flowsheet JSON file")
    validate_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")

    # processforge diagram
    diagram_parser = subparsers.add_parser("diagram", help="Generate a flowsheet diagram")
    diagram_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    diagram_parser.add_argument("--output-dir", "-o", default=".", help="Output directory (default: current directory)")
    diagram_parser.add_argument("--format", "-f", default="png", choices=["png", "svg", "pdf"], help="Output format (default: png)")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        raise SystemExit(1)

    commands = {
        "run": _cmd_run,
        "validate": _cmd_validate,
        "diagram": _cmd_diagram,
    }
    commands[args.command](args)


if __name__ == "__main__":
    main()
