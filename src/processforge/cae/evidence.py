"""Resolve CAE evidence references against a concrete simulation run.

Resolution is deliberately strict: every evidence item must point at an element
that exists in the flowsheet *and* at a value present in the run's Zarr results
(or manifest artifacts). A missing reference aborts the report rather than
silently emitting an unsupported claim.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Any, Iterable

from .loader import CAEError
from .types import CAECase, CAEEvidence, ExtractOp


@dataclass
class ResolvedEvidence:
    """An evidence item paired with the value read from a run."""

    evidence: CAEEvidence
    value: Any
    unit: str = ""
    location: str = ""
    run_id: str = ""
    status: str = "resolved"


def flowsheet_elements(config: dict) -> tuple[set[str], set[str]]:
    """Return the unit and stream names declared (or implied) by a flowsheet.

    Streams are only sometimes listed explicitly; most appear as a unit's
    ``in``/``out``/``inlets`` values, so those are collected too.
    """
    units = config.get("units", {}) or {}
    unit_names = set(units) if isinstance(units, dict) else set()
    stream_names: set[str] = set((config.get("streams", {}) or {}).keys())

    for unit_cfg in units.values() if isinstance(units, dict) else []:
        if not isinstance(unit_cfg, dict):
            continue
        for key in ("in", "out", "inlets"):
            _collect_stream_names(unit_cfg.get(key), stream_names)
    return unit_names, stream_names


def _collect_stream_names(value: Any, target: set[str]) -> None:
    if isinstance(value, str):
        target.add(value)
    elif isinstance(value, dict):
        for item in value.values():
            _collect_stream_names(item, target)
    elif isinstance(value, Iterable) and not isinstance(value, (bytes, bytearray)):
        for item in value:
            _collect_stream_names(item, target)


def resolve_case_evidence(
    case: CAECase,
    *,
    flowsheet_config: dict,
    run_summary: dict,
    manifest: Any | None,
    run_id: str,
    zarr_path: str,
) -> dict[str, ResolvedEvidence]:
    """Resolve every evidence item in *case*, raising on the first failure.

    All failures are collected and reported together so a single report run
    surfaces every broken reference.
    """
    unit_names, stream_names = flowsheet_elements(flowsheet_config)
    resolved: dict[str, ResolvedEvidence] = {}
    errors: list[str] = []

    for evidence in case.evidence:
        try:
            resolved[evidence.id] = _resolve_one(
                evidence,
                unit_names=unit_names,
                stream_names=stream_names,
                run_summary=run_summary,
                manifest=manifest,
                run_id=run_id,
                zarr_path=zarr_path,
            )
        except CAEError as exc:
            errors.append(f"  - {evidence.id}: {exc}")

    if errors:
        raise CAEError(
            "Cannot resolve CAE evidence for run "
            f"'{run_id}':\n" + "\n".join(errors)
        )
    return resolved


def _resolve_one(
    evidence: CAEEvidence,
    *,
    unit_names: set[str],
    stream_names: set[str],
    run_summary: dict,
    manifest: Any | None,
    run_id: str,
    zarr_path: str,
) -> ResolvedEvidence:
    ref = evidence.ref

    if ref.type == "unit":
        return _resolve_unit(evidence, unit_names, run_summary, run_id, zarr_path)
    if ref.type == "stream":
        return _resolve_stream(evidence, stream_names, run_summary, run_id, zarr_path)
    return _resolve_artifact(evidence, unit_names, manifest, run_id)


def _resolve_unit(
    evidence: CAEEvidence,
    unit_names: set[str],
    run_summary: dict,
    run_id: str,
    zarr_path: str,
) -> ResolvedEvidence:
    ref = evidence.ref
    if ref.name not in unit_names:
        raise CAEError(
            f"unit '{ref.name}' is not declared in the flowsheet "
            f"(available: {_sorted(unit_names)})"
        )

    fields = _unit_fields(run_summary, ref.name)
    if fields is None:
        raise CAEError(
            f"unit '{ref.name}' produced no engine output in run '{run_id}'"
        )
    field = _require_field(fields, ref.name, ref.field)
    value, unit = _extract(field, ref.extract, f"unit '{ref.name}.{ref.field}'")
    return ResolvedEvidence(
        evidence=evidence,
        value=value,
        unit=unit,
        location=zarr_path,
        run_id=run_id,
    )


def _resolve_stream(
    evidence: CAEEvidence,
    stream_names: set[str],
    run_summary: dict,
    run_id: str,
    zarr_path: str,
) -> ResolvedEvidence:
    ref = evidence.ref
    if ref.name not in stream_names:
        raise CAEError(
            f"stream '{ref.name}' is not declared or referenced in the flowsheet "
            f"(available: {_sorted(stream_names)})"
        )

    stream = (run_summary.get("streams") or {}).get(ref.name)
    if stream is None:
        raise CAEError(f"stream '{ref.name}' produced no results in run '{run_id}'")
    fields = stream.get("fields") or {}
    field = _require_field(fields, f"stream '{ref.name}'", ref.field)
    value, unit = _extract(field, ref.extract, f"stream '{ref.name}.{ref.field}'")
    return ResolvedEvidence(
        evidence=evidence,
        value=value,
        unit=unit,
        location=zarr_path,
        run_id=run_id,
    )


def _resolve_artifact(
    evidence: CAEEvidence,
    unit_names: set[str],
    manifest: Any | None,
    run_id: str,
) -> ResolvedEvidence:
    ref = evidence.ref
    if manifest is None:
        raise CAEError(f"run manifest for '{run_id}' is unavailable")

    owners: list[tuple[str, str]] = []
    if ref.unit:
        if ref.unit not in unit_names:
            raise CAEError(
                f"artifact owner unit '{ref.unit}' is not declared in the flowsheet"
            )
        owners.append(("unit", ref.unit))
    else:
        owners.extend(("unit", name) for name in getattr(manifest, "units", {}))
        owners.extend(("stream", name) for name in getattr(manifest, "streams", {}))

    for kind, owner in owners:
        collection = (
            getattr(manifest, "units", {})
            if kind == "unit"
            else getattr(manifest, "streams", {})
        )
        output = collection.get(owner)
        if output is None:
            continue
        for artifact in getattr(output, "artifacts", []):
            if artifact.name != ref.name:
                continue
            location = artifact.local_path or (
                artifact.remote_uris[0] if artifact.remote_uris else ""
            )
            status = "resolved"
            if artifact.local_path and not os.path.exists(artifact.local_path):
                if artifact.remote_uris:
                    status = "remote-only"
                else:
                    status = "missing"
            return ResolvedEvidence(
                evidence=evidence,
                value=location,
                unit=artifact.kind,
                location=location,
                run_id=run_id,
                status=status,
            )

    scope = f"unit '{ref.unit}'" if ref.unit else "the run"
    raise CAEError(f"artifact '{ref.name}' was not produced by {scope}")


def _unit_fields(run_summary: dict, name: str) -> dict | None:
    engine = (run_summary.get("engine_outputs") or {}).get(name)
    if engine is not None:
        return engine.get("fields") or {}

    solver = (run_summary.get("solver_units") or {}).get(name)
    if solver is not None:
        return {
            key: {"value": value, "unit": "", "std_dev": None}
            for key, value in (solver.get("attrs") or {}).items()
            if not key.startswith("_")
        }
    return None


def _require_field(fields: dict, owner: str, name: str | None) -> dict:
    field = fields.get(name) if name else None
    if field is None:
        raise CAEError(
            f"{owner} has no field '{name}' "
            f"(available: {_sorted(fields.keys())})"
        )
    return field


def _extract(field: dict, op: ExtractOp, label: str) -> tuple[Any, str]:
    unit = field.get("unit", "") or ""
    raw = field.get("value")

    if op == "std_dev":
        std = field.get("std_dev")
        if std is None:
            raise CAEError(f"{label} has no standard deviation")
        return std, unit

    if op in ("value", "final"):
        if raw is None or (isinstance(raw, (list, tuple)) and not raw):
            raise CAEError(f"{label} has no value")
        return _last_scalar(raw), unit

    numbers = _numeric_values(raw)
    if not numbers:
        raise CAEError(f"{label} has no numeric values for '{op}'")
    if op == "mean":
        return sum(numbers) / len(numbers), unit
    if op == "max":
        return max(numbers), unit
    if op == "min":
        return min(numbers), unit
    raise CAEError(f"unsupported extraction '{op}'")


def _last_scalar(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        if not value:
            raise CAEError("empty sequence")
        return value[-1]
    return value


def _numeric_values(value: Any) -> list[float]:
    if value is None:
        return []
    if isinstance(value, bool):
        return []
    if isinstance(value, (int, float)):
        return [] if math.isnan(value) else [float(value)]
    if isinstance(value, (list, tuple)):
        result: list[float] = []
        for item in value:
            result.extend(_numeric_values(item))
        return result
    return []


def _sorted(values: Iterable[str]) -> str:
    return ", ".join(sorted(str(v) for v in values)) or "none"
