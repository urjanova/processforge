"""CLI entry point for ``pf`` / ``processforge`` commands."""

from __future__ import annotations

import sys

import typer
from loguru import logger

from .cli import register_commands

__all__ = ["main"]


def configure_logging(level: str = "INFO") -> None:
    """Install a concise stderr sink (no timestamp, no module/line).

    Replaces loguru's default handler — which prints
    ``2026-10-22 17:54:26.021 | ERROR | module:function:line - ...`` — with a
    compact ``LEVEL | message`` line. Colours are emitted only when stderr is a
    TTY, so redirected/CI output stays clean.
    """
    logger.remove()
    logger.add(
        sys.stderr,
        level=level,
        format="<level>{level: <8}</level> | <level>{message}</level>",
    )

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
    # ``--debug`` is parsed by typer inside the callback, which runs after this
    # point — scan argv so debug logging is enabled for the whole invocation.
    debug = "--debug" in sys.argv or "-d" in sys.argv
    configure_logging("DEBUG" if debug else "INFO")
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
