"""High-level orchestration for building CAE reports.

Keeps the CLI thin: locate the repository, resolve links, load cases, resolve
evidence and render each report.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass

from .evidence import resolve_case_evidence
from .file_repository import FileCAERepository
from .loader import CAEError, find_cae_root
from .reporting.cae_renderer import render_cae_report
from .types import CAECase


@dataclass
class CAEReport:
    """A rendered report plus the case it came from."""

    case: CAECase
    text: str


def generate_cae_reports(
    *,
    flowsheet_path: str,
    flowsheet_config: dict,
    run_id: str,
    run_summary: dict,
    manifest,
    zarr_path: str,
    jurisdiction: str | None = None,
    cae_root: str | None = None,
    run_timestamp: str = "",
    run_mode: str = "",
) -> list[CAEReport]:
    """Build every CAE report linked to *flowsheet_path*.

    Raises :class:`~processforge.cae.loader.CAEError` when no case is linked,
    the requested jurisdiction has no claims, or any evidence cannot be
    resolved.
    """
    if cae_root:
        root = os.path.abspath(os.path.expanduser(cae_root))
        if not os.path.isfile(os.path.join(root, "index.yaml")):
            raise CAEError(
                f"No CAE repository found at '{root}' (missing index.yaml)."
            )
    else:
        root = find_cae_root(os.path.dirname(os.path.abspath(flowsheet_path)))
        if not root:
            root = find_cae_root()
        if not root:
            raise CAEError(
                "No CAE repository found. Expected a 'cae/index.yaml' in the current "
                "directory or a parent, or set PROCESSFORGE_CAE_DIR."
            )

    repository = FileCAERepository(root)
    links = repository.links_for_flowsheet(flowsheet_path)
    if not links:
        raise CAEError(f"No CAE case linked to flowsheet '{flowsheet_path}'.")

    flowsheet_hash = _flowsheet_hash(flowsheet_config)
    reports: list[CAEReport] = []

    for link in links:
        case = repository.load_case(link.case_id, link.version)
        if jurisdiction:
            case = case.filter_jurisdiction(jurisdiction)
            if not case.claims:
                continue

        evidence = resolve_case_evidence(
            case,
            flowsheet_config=flowsheet_config,
            run_summary=run_summary,
            manifest=manifest,
            run_id=run_id,
            zarr_path=zarr_path,
        )
        text = render_cae_report(
            case,
            evidence,
            flowsheet_path=flowsheet_path,
            run_id=run_id,
            run_timestamp=run_timestamp,
            run_mode=run_mode,
            zarr_path=zarr_path,
            flowsheet_hash=flowsheet_hash,
            jurisdiction=jurisdiction or "",
        )
        reports.append(CAEReport(case=case, text=text))

    if not reports:
        available = sorted(
            {
                j
                for link in links
                for j in repository.load_case(link.case_id, link.version).jurisdictions()
            }
        )
        raise CAEError(
            f"No CAE claims found for jurisdiction '{jurisdiction}' "
            f"(available: {', '.join(available) or 'none'})"
        )
    return reports


def _flowsheet_hash(config: dict) -> str:
    payload = {k: v for k, v in config.items() if not k.startswith("_")}
    encoded = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]
