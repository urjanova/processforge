"""CLI entry point for ``pf`` / ``processforge`` commands."""

from __future__ import annotations

import typer
from loguru import logger

from .cli import register_commands

__all__ = ["main"]

app = typer.Typer(
    help="Processforge - Chemical Process Simulation",
    add_completion=False,
)

_debug: bool = False


@app.callback(invoke_without_command=True)
def _cli(
    ctx: typer.Context,
    debug: bool = typer.Option(
        False,
        "--debug",
        "-d",
        help="Enable detailed debug tracebacks",
    ),
) -> None:
    global _debug
    _debug = debug
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise SystemExit(1)


register_commands(app)


def main() -> None:
    """Parse arguments and dispatch to the appropriate subcommand."""
    try:
        app()
    except SystemExit:
        raise
    except Exception as e:
        logger.debug("Full traceback:", exc_info=True)
        if _debug:
            logger.exception(
                "An unhandled exception occurred during command execution:"
            )
        else:
            logger.error(f"{type(e).__name__}: {str(e)}")
            logger.info("Run with --debug for the full traceback.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
