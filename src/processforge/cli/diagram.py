"""``pf diagram`` — generate a flowsheet diagram."""

from __future__ import annotations

import argparse
import json
import os

from loguru import logger

from ..utils.flowsheet_diagram import draw_flowsheet
from .common import require_existing_file


def add_diagram_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the ``diagram`` subcommand."""
    parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    parser.add_argument(
        "--output-dir",
        "-o",
        default="diagrams",
        help="Output directory (default: diagrams)",
    )
    parser.add_argument(
        "--format",
        "-f",
        default="png",
        choices=["png", "svg", "pdf"],
        help="Output format (default: png)",
    )


def cmd_diagram(args: argparse.Namespace) -> None:
    """Generate a flowsheet diagram from a JSON file."""
    fname = args.flowsheet
    require_existing_file(fname)

    try:
        with open(fname, "r", encoding="utf-8") as f:
            flowsheet_schema = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.error(f"Failed to read flowsheet '{fname}': {e}")
        raise SystemExit(1)

    output_dir = args.output_dir or "diagrams"
    fmt = args.format or "png"

    base_name = os.path.splitext(os.path.basename(fname))[0]

    try:
        output_path = draw_flowsheet(
            flowsheet_schema,
            output_directory=output_dir,
            output_filename=base_name,
            file_format=fmt,
        )
    except Exception as e:
        logger.error(f"Failed to generate diagram: {type(e).__name__}: {e}")
        raise SystemExit(1)

    logger.info(f"Diagram saved to {output_path}")
