"""``pf run`` — run a process simulation from a flowsheet JSON file."""
from __future__ import annotations

import typer
from loguru import logger

from ..runner import ProcessforgeRunError, run_flowsheet


def run(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
    export_images: bool = typer.Option(
        False,
        "--export-images",
        help="Generate PNG plots for simulation outputs",
    ),
) -> None:
    """Run a process simulation from a flowsheet JSON file."""
    try:
        result = run_flowsheet(flowsheet, export_images=export_images)
    except ProcessforgeRunError as exc:
        logger.error(str(exc))
        raise SystemExit(1) from exc

    logger.info("=== Run Summary ===")
    logger.info(f"  Status      : {result.status}")
    logger.info(f"  Run ID      : {result.run_id}")
    logger.info(f"  Archive     : {result.archive_path}")
    if result.remote_uris:
        logger.info(f"  Remote URIs : {len(result.remote_uris)} object(s) in S3")
