"""``pf diagram`` — generate a flowsheet diagram."""

from __future__ import annotations

import json
import os
from typing import Literal

import typer
from loguru import logger

from ..utils.flowsheet_diagram import draw_flowsheet
from .common import require_existing_file


def diagram(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
    output_dir: str = typer.Option(
        "diagrams",
        "--output-dir",
        "-o",
        help="Output directory (default: diagrams)",
    ),
    format: Literal["png", "svg", "pdf"] = typer.Option(
        "png",
        "--format",
        "-f",
        help="Output format (default: png)",
    ),
) -> None:
    """Generate a flowsheet diagram from a JSON file."""
    require_existing_file(flowsheet)

    try:
        with open(flowsheet, "r", encoding="utf-8") as f:
            flowsheet_schema = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.error(f"Failed to read flowsheet '{flowsheet}': {e}")
        raise SystemExit(1)

    output_dir = output_dir or "diagrams"
    fmt = format or "png"

    base_name = os.path.splitext(os.path.basename(flowsheet))[0]

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
