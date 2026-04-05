import argparse
import json
import os

from loguru import logger

from .utils.validate_flowsheet import validate_flowsheet
from .flowsheet import Flowsheet
from .eo import EOFlowsheet
from .provenance import build_dynamic_x0, build_run_info
from .result import (
    plot_results,
    plot_timeseries,
    save_results_zarr,
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

    # Stash the config file path for providers that need it (e.g. ModelicaProvider).
    # Added after schema validation so additionalProperties:false doesn't reject it.
    config["_config_path"] = fname

    base_name = os.path.splitext(os.path.basename(fname))[0]

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    is_dynamic = mode == "dynamic"

    if is_dynamic:
        fs = Flowsheet(config)
        logger.info("=== Dynamic Results ===")
        results = fs.run()
        x0, var_names = build_dynamic_x0(config)
        run_info = build_run_info(config, x0=x0, var_names=var_names)
    else:
        # Pass None so EOFlowsheet resolves backend from config (with scipy default).
        fs = EOFlowsheet(config, backend=None)
        logger.info("=== Steady-State EO Results ===")
        results = fs.run()
        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)

    zarr_path = save_results_zarr(
        results,
        os.path.join("outputs", f"{base_name}_results.zarr"),
        run_info=run_info,
    )
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


def _cmd_export_fmu(args):
    """Export a flowsheet as an FMI 2.0 co-simulation FMU."""
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    from .fmu import build_fmu  # local import — pythonfmu is optional

    output_dir = args.output_dir or "outputs"
    backend = args.backend or "scipy"

    try:
        fmu_path = build_fmu(fname, output_dir=output_dir, backend=backend)
        logger.info(f"FMU written to: {fmu_path}")
    except Exception as e:
        logger.error(f"FMU export failed: {e}")
        raise SystemExit(1)


def _cmd_export_modelica(args):
    """Transpile a flowsheet to Modelica and optionally compile via OMPython."""
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    from .modelica import transpile, compile_modelica
    from .modelica.transpiler import _derive_model_name
    from .utils.validate_flowsheet import validate_flowsheet as _vf

    output_dir = args.output_dir or "outputs"

    try:
        mo_path = transpile(fname, output_dir=output_dir)
    except Exception as e:
        logger.error(f"Modelica transpilation failed: {e}")
        raise SystemExit(1)

    if not args.no_compile:
        try:
            config = _vf(fname)
            model_name = _derive_model_name(config, fname)
            fmu_path = compile_modelica(mo_path, model_name, output_dir=output_dir)
            logger.info(f"Model Exchange FMU written to: {fmu_path}")
        except Exception as e:
            logger.error(f"OMPython compilation failed: {e}")
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
    # processforge validate
    validate_parser = subparsers.add_parser(
        "validate", help="Validate a flowsheet JSON file"
    )
    validate_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")

    # processforge diagram
    diagram_parser = subparsers.add_parser(
        "diagram", help="Generate a flowsheet diagram"
    )
    diagram_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    diagram_parser.add_argument(
        "--output-dir",
        "-o",
        default=".",
        help="Output directory (default: current directory)",
    )
    diagram_parser.add_argument(
        "--format",
        "-f",
        default="png",
        choices=["png", "svg", "pdf"],
        help="Output format (default: png)",
    )

    # processforge export-modelica
    mo_parser = subparsers.add_parser(
        "export-modelica",
        help="Transpile flowsheet to Modelica .mo and compile via OMPython",
    )
    mo_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    mo_parser.add_argument(
        "--output-dir",
        "-o",
        default="outputs",
        help="Directory for the .mo and .fmu outputs (default: outputs/)",
    )
    mo_parser.add_argument(
        "--no-compile",
        action="store_true",
        help="Write the .mo file only; skip OMPython/omc compilation",
    )

    # processforge export-fmu
    fmu_parser = subparsers.add_parser(
        "export-fmu",
        help="Export flowsheet as FMI 2.0 co-simulation FMU",
    )
    fmu_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    fmu_parser.add_argument(
        "--output-dir",
        "-o",
        default="outputs",
        help="Directory for the output FMU (default: outputs/)",
    )
    fmu_parser.add_argument(
        "--backend",
        choices=["scipy", "pyomo", "casadi"],
        default="scipy",
        help="EO solver backend for steady-state mode (default: scipy)",
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        raise SystemExit(1)

    commands = {
        "run": _cmd_run,
        "validate": _cmd_validate,
        "diagram": _cmd_diagram,
        "export-fmu": _cmd_export_fmu,
        "export-modelica": _cmd_export_modelica,
    }
    commands[args.command](args)


if __name__ == "__main__":
    main()
