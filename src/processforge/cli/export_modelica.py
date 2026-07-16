"""``pf export-modelica`` — transpile flowsheet to Modelica .mo and compile via OMPython."""

from __future__ import annotations

import argparse

from loguru import logger

from .common import require_existing_file


def add_export_modelica_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the ``export-modelica`` subcommand."""
    parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    parser.add_argument(
        "--output-dir",
        "-o",
        default="outputs",
        help="Directory for the .mo and .fmu outputs (default: outputs/)",
    )
    parser.add_argument(
        "--no-compile",
        action="store_true",
        help="Write the .mo file only; skip OMPython/omc compilation",
    )


def cmd_export_modelica(args: argparse.Namespace) -> None:
    """Transpile a flowsheet to Modelica and optionally compile via OMPython."""
    fname = args.flowsheet
    require_existing_file(fname)

    from ..modelica import transpile, compile_modelica
    from ..modelica.transpiler import _derive_model_name
    from ..utils.validate_flowsheet import validate_flowsheet as _vf

    output_dir = args.output_dir or "outputs"

    try:
        mo_path = transpile(fname, output_dir=output_dir)
    except Exception as e:
        logger.error(f"Modelica transpilation failed: {type(e).__name__}: {e}")
        logger.debug("Transpile traceback:", exc_info=True)
        raise SystemExit(1)

    if not args.no_compile:
        try:
            config = _vf(fname)
            model_name = _derive_model_name(config, fname)
            fmu_path = compile_modelica(mo_path, model_name, output_dir=output_dir)
            logger.info(f"Model Exchange FMU written to: {fmu_path}")
        except Exception as e:
            logger.error(f"OMPython compilation failed: {type(e).__name__}: {e}")
            logger.debug("Compilation traceback:", exc_info=True)
            raise SystemExit(1)
