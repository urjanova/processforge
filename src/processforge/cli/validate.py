"""``pf validate`` — check providers and environment are ready for a flowsheet."""

from __future__ import annotations

import argparse

from loguru import logger

from .common import check_providers, require_existing_file, validate_runtime_flowsheet


def add_validate_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the ``validate`` subcommand."""
    parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")


def cmd_validate(args: argparse.Namespace) -> None:
    """Check that providers and environment are ready for a flowsheet."""
    fname = args.flowsheet
    require_existing_file(fname, label="Flowsheet file")

    config = validate_runtime_flowsheet(fname)
    check_providers(config, fname, fail_fast=False, verbose=True)

    logger.info("All providers ready.")
