"""FMU interface description for SolverUnit-driven flowsheets.

This module derives the logical input/output ports that a flowsheet exposes
to PathSim / FMI consumers.  It is deliberately independent of PythonFMU so it
can be reused by the adapter and by manifest generation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


def _sanitize_name(s: str) -> str:
    """Replace characters invalid in Python identifiers with underscores."""
    import re
    return re.sub(r"[^a-zA-Z0-9_]", "_", s)


@dataclass
class SolverUnitPort:
    """A single FMI variable associated with a SolverUnit."""

    unit_name: str
    path: str
    attr_name: str
    causality: str
    variability: str
    initial_value: float
    description: str
    unit: str = ""


def _get_fmu_block(unit_cfg: dict) -> dict:
    """Return the unit's ``fmu`` block, normalised to a dict."""
    block = unit_cfg.get("fmu") or unit_cfg.get("extra", {}).get("fmu")
    if isinstance(block, dict):
        return block
    return {}


def _get_nested_value(root: dict, dotted_path: str, default: Any = None) -> Any:
    """Read ``a.b.c`` from nested dicts."""
    cur = root
    for part in dotted_path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def _derive_openmc_default_outputs(unit_cfg: dict) -> list[str]:
    """Best-effort default output names for an OpenMC SolverUnit."""
    sim_type = unit_cfg.get("sim_type", "")
    outputs: list[str] = []
    if "eigenvalue" in sim_type or sim_type in ("eigenvalue_csg", "reactor_core"):
        outputs.append("k_eff")
        outputs.append("power")

    mesh_tallies = _get_nested_value(unit_cfg, "solver_config.mesh_tallies") or []
    if not mesh_tallies:
        # The OpenMC provider synthesises a default tally id=1 with scores flux/fission.
        mesh_tallies = [{"tally_id": 1, "scores": ["flux", "fission"]}]
    for t in mesh_tallies:
        tid = t.get("tally_id", 1)
        for score in t.get("scores", []):
            outputs.append(f"tally_{tid}_{score}_mean")
            outputs.append(f"tally_{tid}_{score}_integrated")
    return outputs


def _derive_festim_default_outputs(unit_cfg: dict) -> list[str]:
    """Best-effort default output names for a FESTIM SolverUnit.

    Mirrors the naming produced by ``FestimProvider._extract_results``.
    """
    exports = _get_nested_value(unit_cfg, "solver_config.exports") or []
    outputs: list[str] = []
    for ex in exports:
        ex_type = ex.get("type", "")
        key: str | None = None
        if ex_type in ("surface_flux", "total_volume"):
            key = ex.get("filename", "")
            if key:
                key = key.replace(".csv", "")
        elif ex_type == "profile_1d":
            key = ex.get("field", "")
        if key:
            outputs.append(key)
            outputs.append(f"{key}_mean_total")
            outputs.append(f"{key}_std_dev")
    return outputs


def derive_solverunit_outputs(unit_name: str, unit_cfg: dict) -> list[str]:
    """Return the list of output field names for a SolverUnit."""
    explicit = _get_fmu_block(unit_cfg).get("outputs")
    if explicit is not None:
        return list(explicit)

    provider = unit_cfg.get("provider", "").lower()
    if provider == "openmc":
        return _derive_openmc_default_outputs(unit_cfg)
    if provider == "festim":
        return _derive_festim_default_outputs(unit_cfg)
    return []


def derive_solverunit_inputs(unit_name: str, unit_cfg: dict) -> list[str]:
    """Return the list of dotted input paths declared in the ``fmu`` block."""
    return list(_get_fmu_block(unit_cfg).get("inputs", []))


def build_solverunit_interface(
    config: dict,
) -> dict[str, list[SolverUnitPort]]:
    """Return {unit_name: [port, ...]} for every SolverUnit in the config."""
    interface: dict[str, list[SolverUnitPort]] = {}
    for unit_name, unit_cfg in config.get("units", {}).items():
        if unit_cfg.get("type") != "SolverUnit":
            continue

        ports: list[SolverUnitPort] = []
        safe_unit = _sanitize_name(unit_name)

        for dotted in derive_solverunit_inputs(unit_name, unit_cfg):
            safe_path = dotted.replace(".", "_")
            attr = f"in_{safe_unit}_{safe_path}"
            initial = float(_get_nested_value(unit_cfg, dotted, 0.0) or 0.0)
            ports.append(
                SolverUnitPort(
                    unit_name=unit_name,
                    path=dotted,
                    attr_name=attr,
                    causality="input",
                    variability="continuous",
                    initial_value=initial,
                    description=f"Input '{dotted}' for SolverUnit '{unit_name}'",
                )
            )

        for field_name in derive_solverunit_outputs(unit_name, unit_cfg):
            safe_field = _sanitize_name(field_name)
            attr = f"out_{safe_unit}_{safe_field}"
            ports.append(
                SolverUnitPort(
                    unit_name=unit_name,
                    path=field_name,
                    attr_name=attr,
                    causality="output",
                    variability="continuous",
                    initial_value=0.0,
                    description=f"Output '{field_name}' from SolverUnit '{unit_name}'",
                )
            )

        interface[unit_name] = ports
    return interface


def build_fmu_interface_manifest(
    config: dict,
    stream_specs: list[dict],
    solverunit_ports: dict[str, list[SolverUnitPort]],
    parameter_specs: list[dict] | None = None,
) -> dict:
    """Build the JSON-serialisable ``fmu_interface.json`` manifest."""
    manifest: dict = {
        "version": "1.0",
        "flowsheet_name": config.get("metadata", {}).get("name", ""),
        "variables": [],
    }

    for spec in stream_specs:
        manifest["variables"].append(
            {
                "attr_name": spec["attr_name"],
                "logical_name": spec.get("logical_name", spec["attr_name"]),
                "causality": spec["causality"],
                "variability": spec["variability"],
                "initial_value": spec["initial_value"],
                "description": spec["description"],
                "unit": spec.get("unit", ""),
            }
        )

    if parameter_specs:
        for spec in parameter_specs:
            manifest["variables"].append(
                {
                    "attr_name": spec["attr_name"],
                    "logical_name": spec.get("logical_name", spec["attr_name"]),
                    "causality": "parameter",
                    "variability": "fixed",
                    "initial_value": spec["initial_value"],
                    "description": spec["description"],
                    "unit": spec.get("unit", ""),
                }
            )

    for unit_name, ports in solverunit_ports.items():
        for port in ports:
            manifest["variables"].append(
                {
                    "attr_name": port.attr_name,
                    "logical_name": f"{unit_name}.{port.path}",
                    "causality": port.causality,
                    "variability": port.variability,
                    "initial_value": port.initial_value,
                    "description": port.description,
                    "unit": port.unit,
                }
            )

    return manifest
