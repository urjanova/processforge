"""OutputCollector — assemble a standardized :class:`RunManifest` from a run.

This is the single point where every engine's outputs are gathered into the
canonical, unit-aware record persisted by :class:`ProcessStateArchive`:

* **Solver-unit** outputs (OpenMC, FESTIM, …) arrive as
  :class:`EngineOutput` objects already produced by the providers.
* **Stream thermo** outputs (CoolProp / Cantera) are produced by calling the
  thermo provider's ``get_thermo_properties`` on each final stream state and
  wrapping the ``{H, Cp, K_values}`` dict into a unit-tagged
  :class:`StreamOutput`.  This makes CoolProp/Cantera first-class in the
  standardized store without touching the hot residual paths that consume the
  raw dict.
"""
from __future__ import annotations

from typing import Optional

from loguru import logger

from .types import (
    EngineOutput,
    OutputField,
    OutputProvenance,
    Quantity,
    RunManifest,
    StreamOutput,
)


def collect_stream_outputs(provider, stream_results: dict) -> dict[str, StreamOutput]:
    """Build :class:`StreamOutput` for each stream using *provider*'s thermo."""
    if provider is None:
        return {}
    out: dict[str, StreamOutput] = {}
    for name, sdata in stream_results.items():
        if not isinstance(sdata, dict):
            continue
        z = sdata.get("z")
        T = sdata.get("T")
        P = sdata.get("P")
        if not isinstance(z, dict) or T is None or P is None:
            continue
        try:
            props = provider.get_thermo_properties({"z": z, "T": T, "P": P})
        except Exception as exc:  # noqa: BLE001
            logger.debug(f"thermo properties unavailable for stream '{name}': {exc}")
            continue
        if not isinstance(props, dict):
            continue
        fields: list[OutputField] = []
        if props.get("H") is not None:
            fields.append(OutputField(
                name="H", quantity=Quantity(value=float(props["H"]), unit="J/mol"),
                kind="scalar", source="PropsSI(H)",
            ))
        if props.get("Cp") is not None:
            fields.append(OutputField(
                name="Cp", quantity=Quantity(value=float(props["Cp"]), unit="J/(mol*K)"),
                kind="scalar", source="PropsSI(Cp)",
            ))
        K = props.get("K_values") or {}
        if K:
            fields.append(OutputField(
                name="K_values",
                quantity=Quantity(value=[float(v) for v in K.values()], unit=""),
                kind="vector", source="PropsSI(K)",
            ))
        if fields:
            engine = type(provider).__name__.replace("Provider", "").lower()
            out[name] = StreamOutput(
                stream=name, engine=engine, fields=fields,
                provenance=OutputProvenance(),
            )
    return out


def build_run_manifest(
    *,
    run_id: str,
    mode: str,
    flowsheet_name: str,
    provenance: dict,
    engine_outputs: dict[str, EngineOutput],
    stream_results: dict,
    thermo_provider=None,
) -> RunManifest:
    """Assemble a :class:`RunManifest` from a finished run."""
    streams = collect_stream_outputs(thermo_provider, stream_results)
    return RunManifest(
        run_id=run_id,
        mode=mode,
        flowsheet=flowsheet_name,
        units=dict(engine_outputs),
        streams=streams,
        provenance=provenance or {},
    )
