"""Plain-text rendering of a CAE assurance report.

The renderer is deliberately dependency-free (no rich/markdown) so its output
can be piped, diffed, or pasted into a regulatory submission without carrying
terminal formatting.
"""

from __future__ import annotations

import datetime
from typing import Any

from ..evidence import ResolvedEvidence
from ..types import CAECase

WIDTH = 80
LINE = "=" * WIDTH


def render_cae_report(
    case: CAECase,
    evidence: dict[str, ResolvedEvidence],
    *,
    flowsheet_path: str,
    run_id: str,
    run_timestamp: str = "",
    run_mode: str = "",
    zarr_path: str = "",
    flowsheet_hash: str = "",
    jurisdiction: str = "",
) -> str:
    """Render one CAE case as a plain-text report."""
    lines: list[str] = [
        LINE,
        "CAE ASSURANCE REPORT",
        LINE,
        f"Case            : {case.case.title} ({case.case.id} v{case.case.version})",
    ]
    if jurisdiction:
        lines.append(f"Jurisdiction    : {jurisdiction}")
    elif case.jurisdictions():
        lines.append(f"Jurisdiction(s) : {', '.join(case.jurisdictions())}")
    lines.extend(
        [
            f"Flowsheet       : {flowsheet_path}",
            f"Run ID          : {run_id}",
            f"Run timestamp   : {run_timestamp or 'unknown'}",
            f"Run mode        : {run_mode or 'unknown'}",
            f"Zarr store      : {zarr_path or 'unknown'}",
            f"Flowsheet hash  : {flowsheet_hash or 'unavailable'}",
            f"CAE schema      : {case.schema_version}",
            f"Generated (UTC) : {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
            "",
        ]
    )

    lines.extend(_render_case_body(case, evidence))
    return "\n".join(lines)


def _render_case_body(
    case: CAECase, evidence: dict[str, ResolvedEvidence]
) -> list[str]:
    lines: list[str] = []

    if case.requirements:
        lines.append("REQUIREMENTS")
        lines.append("-" * len("REQUIREMENTS"))
        for requirement in case.requirements:
            lines.append(f"[{requirement.id}] {requirement.statement}")
            if requirement.source:
                lines.append(f"  Source : {requirement.source}")
            lines.append("")
    else:
        lines.append("REQUIREMENTS")
        lines.append("-" * len("REQUIREMENTS"))
        lines.append("(none)")
        lines.append("")

    lines.append("CLAIMS")
    lines.append("-" * len("CLAIMS"))
    if not case.claims:
        lines.append("(none)")
    for claim in case.claims:
        header = f"[{claim.id}]"
        if claim.jurisdiction:
            header += f" ({claim.jurisdiction})"
        lines.append(f"{header} {claim.statement}")
        if claim.summary:
            lines.append(f"  Summary      : {claim.summary}")
        lines.append(
            "  Requirements : "
            + (", ".join(claim.requirement_refs) if claim.requirement_refs else "none")
        )
        lines.append(
            "  Arguments    : "
            + (", ".join(claim.argument_refs) if claim.argument_refs else "none")
        )
        lines.append("")

    lines.append("ARGUMENTS")
    lines.append("-" * len("ARGUMENTS"))
    if not case.arguments:
        lines.append("(none)")
    claim_lookup = {c.id: c for c in case.claims}
    for argument in case.arguments:
        supporters = [
            claim.id for claim in claim_lookup.values() if argument.id in claim.argument_refs
        ]
        lines.append(f"[{argument.id}] {argument.statement}")
        lines.append(
            "  Supports : " + (", ".join(supporters) if supporters else "none")
        )
        lines.append(
            "  Evidence : "
            + (", ".join(argument.evidence_refs) if argument.evidence_refs else "none")
        )
        lines.append("")

    lines.append("EVIDENCE")
    lines.append("-" * len("EVIDENCE"))
    if not case.evidence:
        lines.append("(none)")
    for item in case.evidence:
        resolved = evidence.get(item.id)
        lines.append(f"[{item.id}] {item.description}")
        if resolved is None:
            lines.append("  Status : UNRESOLVED")
            lines.append("")
            continue
        ref = item.ref
        source = f"{ref.type} '{ref.name}'"
        if ref.field:
            source += f" field '{ref.field}'"
        source += f" ({ref.extract})"
        lines.append(f"  Source   : {source}")
        lines.append(f"  Value    : {_fmt_value(resolved.value)}")
        if resolved.unit:
            lines.append(f"  Unit     : {resolved.unit}")
        lines.append(f"  Status   : {resolved.status}")
        lines.append(f"  Run      : {resolved.run_id}")
        if resolved.location:
            lines.append(f"  Location : {resolved.location}")
        lines.append("")

    return lines


def _fmt_value(value: Any) -> str:
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_fmt_value(item) for item in value) + "]"
    if value is None:
        return "none"
    return str(value)
