import argparse
import json
import os

from loguru import logger

from .utils.validate_flowsheet import validate_flowsheet
from .flowsheet import Flowsheet
from .eo import EOFlowsheet
from .result import (
    generate_validation_excel,
    plot_results,
    plot_timeseries,
    save_results_zarr,
    upload_zarr_to_s3,
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

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    is_dynamic = mode == "dynamic"

    if is_dynamic:
        fs = Flowsheet(config)
        logger.info("=== Dynamic Results ===")
    else:
        backend = sim_cfg.get("backend", "scipy")
        fs = EOFlowsheet(config, backend=backend)
        logger.info("=== Steady-State EO Results ===")

    results = fs.run()
    base_name = os.path.splitext(os.path.basename(fname))[0]
    zarr_path = save_results_zarr(
        results,
        os.path.join("outputs", f"{base_name}_results.zarr"),
    )
    validation_path = os.path.join("outputs", f"{base_name}_validation.xlsx")
    generate_validation_excel(
        data_source=zarr_path,
        output_filename=validation_path,
    )
    if args.upload_to_s3:
        s3_uri = upload_zarr_to_s3(zarr_path)
        if s3_uri:
            logger.info(f"Zarr results available at {s3_uri}")
    if args.export_images:
        plot_results(results, fname=f"{base_name}_results.png")
        plot_timeseries(results, fname=f"{base_name}_timeseries.png")


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
    run_parser.add_argument(
        "--export-images",
        action="store_true",
        help="Generate PNG plots for simulation outputs",
    )
    run_parser.add_argument(
        "--upload-to-s3",
        action="store_true",
        help=(
            "Upload Zarr results to S3 and remove the local copy. "
            "Requires env vars: S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_KEY. "
            "Optional: S3_ENDPOINT_URL, S3_REGION_NAME."
        ),
    )

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
