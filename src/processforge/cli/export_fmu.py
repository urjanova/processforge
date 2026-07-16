"""``pf export-fmu`` — export flowsheet as FMI 2.0 co-simulation FMU."""

from __future__ import annotations

from typing import Literal

import typer
from loguru import logger

from .common import require_existing_file


def export_fmu(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
    output_dir: str = typer.Option(
        "outputs",
        "--output-dir",
        "-o",
        help="Directory for the output FMU (default: outputs/)",
    ),
    backend: Literal["scipy", "pyomo", "casadi"] = typer.Option(
        "scipy",
        "--backend",
        help="EO solver backend for steady-state mode (default: scipy)",
    ),
) -> None:
    """Export a flowsheet as an FMI 2.0 co-simulation FMU."""
    require_existing_file(flowsheet)

    from ..fmu import build_fmu  # local import — pythonfmu is optional

    output_dir = output_dir or "outputs"
    backend = backend or "scipy"

    try:
        fmu_path = build_fmu(flowsheet, output_dir=output_dir, backend=backend)
        logger.info(f"FMU written to: {fmu_path}")
    except Exception as e:
        logger.error(f"FMU export failed: {type(e).__name__}: {e}")
        logger.debug("FMU export traceback:", exc_info=True)
        raise SystemExit(1)
