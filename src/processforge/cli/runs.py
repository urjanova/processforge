"""``pf runs`` — list runs and inspect their Zarr result summaries."""

from __future__ import annotations

import json
import os

import arrow
import typer
from loguru import logger

from .common import flowsheet_basename, output_root
from ..persistence.archive import ProcessStateArchive
from ..result import _fmt_scientific, _one_line_summary, summarize_zarr_store


def _humanize_timestamp(ts: str) -> str:
    """Render a stored timestamp as a human-readable relative string.

    Falls back to the raw string for values that aren't parseable ISO-8601
    (e.g. ``20260101T000000Z`` used in some test fixtures).
    """
    try:
        return arrow.get(ts).humanize()
    except (arrow.parser.ParserError, ValueError, TypeError):
        return ts


def _zarr_path(archive_path: str, run_id: str) -> str:
    """Return the per-run Zarr store path."""
    return os.path.join(archive_path, "results", run_id, "results.zarr")


def _artifact_status(artifact: dict) -> str:
    """Return a short presence marker for an artifact dict."""
    local_path = artifact.get("local_path")
    remote_uris = artifact.get("remote_uris") or []
    if local_path and os.path.exists(local_path):
        return "local✓"
    if remote_uris:
        return "remote-only"
    return "missing"


def _format_run_details(manifest, summary: dict) -> str:
    """Format a detailed, human-readable result summary for a single run."""
    lines = [
        f"Run: {manifest.run_id}",
        f"Timestamp: {manifest.timestamp}",
        f"Mode: {summary.get('mode', manifest.mode)}",
        "",
        "Results:",
        "────────",
    ]

    engine_outputs = summary.get("engine_outputs", {})
    if engine_outputs:
        for name, eo in engine_outputs.items():
            engine = eo.get("engine", "")
            sim_type = eo.get("sim_type", "")
            status = eo.get("status", "")
            lines.append(f"{name} ({engine} / {sim_type}) [{status}]")
            for fname, fdata in eo.get("fields", {}).items():
                val = _fmt_scientific(fdata.get("value"))
                unit = fdata.get("unit", "")
                std = fdata.get("std_dev")
                if std is not None:
                    lines.append(f"  {fname}: {val} ± {_fmt_scientific(std)} {unit}".strip())
                else:
                    lines.append(f"  {fname}: {val} {unit}".strip())
            if eo.get("attrs") and not eo.get("fields"):
                for k, v in eo["attrs"].items():
                    if not k.startswith("_") and k not in ("engine", "sim_type", "status"):
                        lines.append(f"  {k}: {_fmt_scientific(v)}")
            lines.append("")
    else:
        lines.append("  (no solver/engine outputs)")
        lines.append("")

    streams = summary.get("streams", {})
    if streams:
        lines.append("Streams:")
        all_vars: set[str] = set()
        for sdata in streams.values():
            all_vars.update(sdata.get("variables", []))
        all_vars = sorted(all_vars)

        header_parts = ["stream".ljust(20)] + [v.rjust(14) for v in all_vars]
        lines.append("  " + " ".join(header_parts))
        lines.append("  " + "-" * len("  ".join(header_parts)))
        for sname in sorted(streams):
            sdata = streams[sname]
            fields = sdata.get("fields", {})
            row = [sname.ljust(20)]
            for var in all_vars:
                fdata = fields.get(var, {})
                val = fdata.get("value")
                row.append(_fmt_scientific(val).rjust(14))
            lines.append("  " + " ".join(row))
        lines.append("")

    artifacts = summary.get("artifacts", [])
    if artifacts:
        lines.append("Artifacts:")
        for art in artifacts:
            name = art.get("name", "unknown")
            kind = art.get("kind", "")
            status = _artifact_status(art)
            loc = art.get("local_path") or (art.get("remote_uris") or ["no location"])[0]
            lines.append(f"  [{kind}] {name}: {status} ({loc})")
    else:
        lines.append("Artifacts: none")

    zarr_path = summary.get("_zarr_path")
    if zarr_path:
        lines.append("")
        lines.append(f"Zarr store: {zarr_path}")

    return "\n".join(lines)


def runs(
    flowsheet: str = typer.Argument(..., help="Path to the flowsheet JSON"),
    run_id: str = typer.Argument(
        None,
        help="Show this run's result summary. Use 'latest' for the most recent run.",
    ),
    schema: bool = typer.Option(
        False,
        "--schema",
        help="Print the result schema JSON instead of the summary.",
    ),
):
    """List runs for a flowsheet, or inspect one run's Zarr result summary."""
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
        resolved_id = latest_id if run_id == "latest" and latest_id else run_id
        if resolved_id != run_id and latest_id is None:
            logger.error(f"No latest run found in {archive_path}.")
            raise typer.Exit(code=1)

        if resolved_id not in archive.list_runs():
            logger.error(f"Run '{resolved_id}' not found in {archive_path}/runs.")
            raise typer.Exit(code=1)

        zarr_path = _zarr_path(archive_path, resolved_id)

        if schema:
            schema_path = zarr_path + ".schema.json"
            if not os.path.isfile(schema_path):
                logger.error(f"No schema found for run '{resolved_id}'.")
                raise typer.Exit(code=1)
            with open(schema_path, "r", encoding="utf-8") as f:
                typer.echo(f.read())
            return

        manifest = archive.load_run(resolved_id)
        summary = summarize_zarr_store(zarr_path)
        summary["_zarr_path"] = zarr_path
        typer.echo(_format_run_details(manifest, summary))
        return

    # List all runs with a one-line summary.
    run_ids = archive.list_runs()
    if not run_ids:
        typer.echo(f"No runs found for '{flowsheet}'.")
        return

    header = f"{'RUN ID':<28} {'TIMESTAMP':<22} {'LATEST':<7} {'SUMMARY':<30}"
    typer.echo(header)
    typer.echo("-" * len(header))
    for rid in run_ids:
        manifest = archive.load_run(rid)
        if manifest is None:
            continue
        marker = "*" if rid == latest_id else " "
        zarr_path = _zarr_path(archive_path, rid)
        summary = summarize_zarr_store(zarr_path)
        one_liner = _one_line_summary(summary)
        typer.echo(
            f"{rid:<28} {_humanize_timestamp(manifest.timestamp):<22} {marker:<7} {one_liner:<30}"
        )
    typer.echo("")
    typer.echo(
        "'*' marks the latest run. "
        "'pf runs <flowsheet> <run_id>' shows the result summary. "
        "Use 'latest' as the run id for the most recent run."
    )
