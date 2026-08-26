"""``pf runs`` — list runs for a flowsheet, or show one run's full manifest."""

from __future__ import annotations

import json
import os

import arrow
import typer
from loguru import logger

from .common import flowsheet_basename, output_root
from ..persistence.archive import ProcessStateArchive


def _humanize_timestamp(ts: str) -> str:
    """Render a stored timestamp as a human-readable relative string.

    Falls back to the raw string for values that aren't parseable ISO-8601
    (e.g. ``20260101T000000Z`` used in some test fixtures).
    """
    try:
        return arrow.get(ts).humanize()
    except (arrow.parser.ParserError, ValueError, TypeError):
        return ts


def runs(
    flowsheet: str = typer.Argument(..., help="Path to the flowsheet JSON"),
    run_id: str = typer.Argument(
        None,
        help="Show this run's full manifest instead of listing all runs",
    ),
):
    """List runs stored in a flowsheet's ``.pfarchive`` (or show one run)."""
    base = flowsheet_basename(flowsheet)
    archive_path = os.path.join(output_root(), f"{base}.pfarchive")

    if not os.path.isdir(archive_path):
        logger.error(
            f"No archive found for '{flowsheet}' at {archive_path}. "
            "Run 'pf run' first."
        )
        raise typer.Exit(code=1)

    archive = ProcessStateArchive(archive_path)
    latest = archive.latest_run()
    latest_id = latest.run_id if latest is not None else None

    if run_id is not None:
        if run_id not in archive.list_runs():
            logger.error(f"Run '{run_id}' not found in {archive_path}/runs.")
            raise typer.Exit(code=1)
        manifest = archive.load_run(run_id)
        typer.echo(json.dumps(manifest.model_dump(), indent=2, default=str))
        typer.echo("\nArtifacts on disk:")
        for unit_name, out in manifest.units.items():
            for art in out.artifacts:
                present = os.path.exists(art.local_path) if art.local_path else False
                remote = bool(art.remote_uris)
                status = (
                    "local✓" if present else ("remote-only" if remote else "missing")
                )
                loc = art.local_path or "no local path"
                typer.echo(f"  [{unit_name}] {art.name}: {status} ({loc})")
        return

    run_ids = archive.list_runs()
    if not run_ids:
        typer.echo(f"No runs found for '{flowsheet}'.")
        return

    header = f"{'RUN ID':<28} {'TIMESTAMP':<22} LATEST"
    typer.echo(header)
    typer.echo("-" * len(header))
    for rid in run_ids:
        manifest = archive.load_run(rid)
        if manifest is None:
            continue
        marker = "*" if rid == latest_id else " "
        typer.echo(
            f"{rid:<28} {_humanize_timestamp(manifest.timestamp):<22} {marker}"
        )
    typer.echo("")
    typer.echo(
        "'*' marks the latest run. "
        "'pf runs <flowsheet> <run_id>' prints the full manifest."
    )
