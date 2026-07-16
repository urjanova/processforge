"""``pf export-modelica`` — transpile flowsheet to Modelica .mo and compile via OMPython."""

from __future__ import annotations

import typer
from loguru import logger

from .common import require_existing_file


def export_modelica(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
    output_dir: str = typer.Option(
        "outputs",
        "--output-dir",
        "-o",
        help="Directory for the .mo and .fmu outputs (default: outputs/)",
    ),
    no_compile: bool = typer.Option(
        False,
        "--no-compile",
        help="Write the .mo file only; skip OMPython/omc compilation",
    ),
) -> None:
    """Transpile a flowsheet to Modelica and optionally compile via OMPython."""
    require_existing_file(flowsheet)

    from ..modelica import transpile, compile_modelica
    from ..modelica.transpiler import _derive_model_name
    from ..utils.validate_flowsheet import validate_flowsheet as _vf

    output_dir = output_dir or "outputs"

    try:
        mo_path = transpile(flowsheet, output_dir=output_dir)
    except Exception as e:
        logger.error(f"Modelica transpilation failed: {type(e).__name__}: {e}")
        logger.debug("Transpile traceback:", exc_info=True)
        raise SystemExit(1)

    if not no_compile:
        try:
            config = _vf(flowsheet)
            model_name = _derive_model_name(config, flowsheet)
            fmu_path = compile_modelica(mo_path, model_name, output_dir=output_dir)
            logger.info(f"Model Exchange FMU written to: {fmu_path}")
        except Exception as e:
            logger.error(f"OMPython compilation failed: {type(e).__name__}: {e}")
            logger.debug("Compilation traceback:", exc_info=True)
            raise SystemExit(1)
