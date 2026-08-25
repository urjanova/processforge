"""``pf runs`` — list runs for a flowsheet, or show one run's full manifest."""

from __future__ import annotations

import json
import os

import typer
from loguru import logger

from .common import flowsheet_basename, output_root
from ..persistence.archive import ProcessStateArchive


def _summary_for(manifest) -> str:
    """One-line summary of a run's primary unit fields (e.g. keff)."""
    parts = []
    for unit_name, out in manifest.units.items():
        for f in out.fields:
            val = f.quantity.value
            unit = f.quantity.unit or ""
            parts.append(f"{unit_name}.{f.name}={val}{(' ' + unit) if unit else ''}")
    return "  ".join(parts) if parts else "-"


def _zarr_status(archive_path: str, run_id: str) -> str:
    zarr_path = os.path.join(archive_path, "results", run_id, "results.zarr")
    return "✓" if os.path.exists(zarr_path) else "✗"


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

    header = (
        f"{'RUN ID':<28} {'TIMESTAMP':<22} {'MODE':<8} "
        f"{'ZARR':<5} {'LATEST':<6} SUMMARY"
    )
    typer.echo(header)
    typer.echo("-" * len(header))
    for rid in run_ids:
        manifest = archive.load_run(rid)
        if manifest is None:
            continue
        marker = "*" if rid == latest_id else " "
        typer.echo(
            f"{rid:<28} {manifest.timestamp:<22} {manifest.mode:<8} "
            f"{_zarr_status(archive_path, rid):<5} {marker:<6} "
            f"{_summary_for(manifest)}"
        )
    typer.echo("")
    typer.echo(
        "Zarr: ✓ present on disk, ✗ deleted/missing. '*' marks the latest run. "
        "'pf runs <flowsheet> <run_id>' prints the full manifest."
    )
