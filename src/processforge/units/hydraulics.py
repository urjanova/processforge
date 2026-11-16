"""Shared hydraulic property resolution and unit math.

All fluid properties (Cp, molar volume, density, viscosity) are sourced from
the attached thermodynamic provider or from explicit flowsheet unit
parameters — never from hard-coded fluid defaults. Explicit unit parameters
always win over provider values so flowsheets can override any property.
"""
from __future__ import annotations

import math
from typing import Callable, Optional

from loguru import logger


def _finite_float(value) -> Optional[float]:
    """Coerce *value* to a finite float, or return ``None``."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def normalize_efficiency(efficiency, unit_name: str) -> float:
    """Return pump efficiency as a fraction in ``(0, 1]``.

    Values greater than 1 are interpreted as percent (e.g. ``36`` -> ``0.36``)
    with a warning, matching flowsheets that specify efficiency in %.
    """
    eff = _finite_float(efficiency)
    if eff is None or eff <= 0:
        raise ValueError(
            f"Unit '{unit_name}': efficiency must be a positive fraction "
            f"(or percent > 1); got {efficiency!r}."
        )
    if eff > 1.0:
        logger.warning(
            f"Unit '{unit_name}': efficiency {eff} > 1 — interpreting as percent."
        )
        eff /= 100.0
    return min(eff, 1.0)


def resolve_properties(
    unit_name: str,
    stream: dict,
    overrides: dict,
    thermo: Optional[Callable[[dict], dict]] = None,
) -> dict:
    """Resolve mixture hydraulic properties for a unit.

    Precedence: explicit flowsheet unit parameters (*overrides*) first, then
    the provider's ``get_thermo_properties``. Returns a dict with keys
    ``cp`` [J/mol/K], ``molar_volume`` [m³/mol] and ``viscosity`` [Pa·s];
    entries that cannot be resolved are ``None``.
    """
    cp = _finite_float(overrides.get("cp"))
    molar_volume = _finite_float(overrides.get("molar_volume"))
    viscosity = _finite_float(overrides.get("viscosity"))
    if molar_volume is None:
        rho = _finite_float(overrides.get("density"))
        mw = _finite_float(overrides.get("molar_mass"))
        if rho and mw and rho > 0:
            molar_volume = mw / rho

    if thermo is not None:
        try:
            props = thermo(stream)
        except Exception as exc:  # noqa: BLE001 — fall through to overrides
            logger.warning(
                f"Unit '{unit_name}': provider property lookup failed ({exc})."
            )
            props = {}
        if cp is None:
            cp = _finite_float(props.get("Cp"))
        if molar_volume is None:
            molar_volume = _finite_float(props.get("molar_volume"))
            if molar_volume is None:
                rho = _finite_float(props.get("density"))
                mw = _finite_float(props.get("molar_mass"))
                if rho and mw and rho > 0:
                    molar_volume = mw / rho
        if viscosity is None:
            viscosity = _finite_float(props.get("viscosity"))

    return {"cp": cp, "molar_volume": molar_volume, "viscosity": viscosity}


def require_pump_properties(unit_name: str, props: dict) -> None:
    """Raise a clear error when Cp / molar volume are unavailable."""
    missing = []
    if not props.get("cp") or props["cp"] <= 0:
        missing.append('"cp" (J/mol/K)')
    if not props.get("molar_volume") or props["molar_volume"] <= 0:
        missing.append('"molar_volume" (m³/mol) or "density"+"molar_mass"')
    if missing:
        raise RuntimeError(
            f"Unit '{unit_name}' could not resolve mixture properties: "
            f"{', '.join(missing)}. Attach a thermodynamic provider that "
            "supplies them (e.g. coolprop) or set them explicitly on the "
            "unit in the flowsheet."
        )


def pump_power_and_temperature_rise(
    flow_mol, molar_volume, cp, delta_p, efficiency
) -> tuple[float, float]:
    """Shaft power [W] and outlet temperature rise [K] for a displacement pump.

    ``power = F_mol · v_molar · ΔP / η`` and the lost work appears as heat:
    ``dT = power · (1 − η) / (F_mol · Cp_molar)``.
    """
    if flow_mol is None or flow_mol <= 1e-12 or not molar_volume or molar_volume <= 0:
        return 0.0, 0.0
    power = flow_mol * molar_volume * delta_p / efficiency
    dT = power * (1.0 - efficiency) / (flow_mol * cp) if cp and cp > 0 else 0.0
    return power, dT


def laminar_pipe_flow(diameter, delta_p, viscosity, length) -> float:
    """Hagen–Poiseuille volumetric flow [m³/s] for a laminar pipe segment."""
    if diameter <= 0 or viscosity <= 0 or length <= 0 or delta_p <= 0:
        return 0.0
    return math.pi * diameter**4 * delta_p / (128.0 * viscosity * length)
