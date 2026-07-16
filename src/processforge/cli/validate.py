"""``pf validate`` — check providers and environment are ready for a flowsheet."""

from __future__ import annotations

import typer
from loguru import logger

from .common import check_providers, require_existing_file, validate_runtime_flowsheet


def validate(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
) -> None:
    """Check that providers and environment are ready for a flowsheet."""
    require_existing_file(flowsheet, label="Flowsheet file")

    config = validate_runtime_flowsheet(flowsheet)
    check_providers(config, flowsheet, fail_fast=False, verbose=True)

    logger.info("All providers ready.")
