"""``pf export-fmu`` — export flowsheet as FMI 2.0 co-simulation FMU."""

from __future__ import annotations

import argparse

from loguru import logger

from .common import require_existing_file


def add_export_fmu_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the ``export-fmu`` subcommand."""
    parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    parser.add_argument(
        "--output-dir",
        "-o",
        default="outputs",
        help="Directory for the output FMU (default: outputs/)",
    )
    parser.add_argument(
        "--backend",
        choices=["scipy", "pyomo", "casadi"],
        default="scipy",
        help="EO solver backend for steady-state mode (default: scipy)",
    )


def cmd_export_fmu(args: argparse.Namespace) -> None:
    """Export a flowsheet as an FMI 2.0 co-simulation FMU."""
    fname = args.flowsheet
    require_existing_file(fname)

    from ..fmu import build_fmu  # local import — pythonfmu is optional

    output_dir = args.output_dir or "outputs"
    backend = args.backend or "scipy"

    try:
        fmu_path = build_fmu(fname, output_dir=output_dir, backend=backend)
        logger.info(f"FMU written to: {fmu_path}")
    except Exception as e:
        logger.error(f"FMU export failed: {type(e).__name__}: {e}")
        logger.debug("FMU export traceback:", exc_info=True)
        raise SystemExit(1)
