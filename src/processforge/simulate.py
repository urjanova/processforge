"""CLI entry point for ``pf`` / ``processforge`` commands."""

from __future__ import annotations

import argparse

from loguru import logger

from .cli import register_commands

__all__ = ["main"]


def main() -> None:
    """Parse arguments and dispatch to the appropriate subcommand."""
    parser = argparse.ArgumentParser(
        description="Processforge - Chemical Process Simulation",
        prog="processforge",
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Enable detailed debug tracebacks for errors",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    register_commands(subparsers)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        raise SystemExit(1)

    try:
        args.func(args)
    except SystemExit:
        raise
    except Exception as e:
        logger.debug("Full traceback:", exc_info=True)
        if getattr(args, "debug", False):
            logger.exception("An unhandled exception occurred during command execution:")
        else:
            logger.error(f"{type(e).__name__}: {str(e)}")
            logger.info("Run with --debug for the full traceback.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
