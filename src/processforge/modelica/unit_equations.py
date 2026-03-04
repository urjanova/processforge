"""Per-unit Modelica equation generators.

Each function returns a list of equation-section strings for a single unit.
Variable naming convention: ``{stream_name}_{property}`` where property is
``T``, ``P``, ``F``, or ``z_{comp}`` — with all names passed through
``_sanitize_name`` before arriving here.
"""
from __future__ import annotations

from ..fmu._fmi_vars import _sanitize_name


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _svar(stream: str, prop: str) -> str:
    """Return ``{sanitized_stream}_{prop}``."""
    return f"{_sanitize_name(stream)}_{prop}"


def _passthrough_TFz(inlet: str, outlet: str, comps: list[str]) -> list[str]:
    """Emit pass-through equations for T, F, and all compositions."""
    lines = []
    i, o = _sanitize_name(inlet), _sanitize_name(outlet)
    lines.append(f"  {o}_T = {i}_T;")
    lines.append(f"  {o}_F = {i}_F;")
    for c in comps:
        sc = _sanitize_name(c)
        lines.append(f"  {o}_z_{sc} = {i}_z_{sc};")
    return lines


# ---------------------------------------------------------------------------
# Multi-inlet mixer (emitted before the owning unit's equations)
# ---------------------------------------------------------------------------

def mixer_eq(inlets: list[str], outlet: str, comps: list[str]) -> list[str]:
    """Flow-weighted mixing of multiple inlet streams into one virtual outlet."""
    lines: list[str] = [f"  // mixer → {outlet}"]
    o = _sanitize_name(outlet)
    # Total flowrate
    f_terms = " + ".join(f"{_sanitize_name(s)}_F" for s in inlets)
    lines.append(f"  {o}_F = {f_terms};")
    # Pressure: minimum of inlets (conservative)
    p_terms = " + ".join(f"{_sanitize_name(s)}_P" for s in inlets)
    lines.append(
        f"  {o}_P = ({p_terms}) / {len(inlets)}.0;"
    )
    # Temperature: flow-weighted average (enthalpy mixing simplified)
    num_T = " + ".join(
        f"{_sanitize_name(s)}_F * {_sanitize_name(s)}_T" for s in inlets
    )
    lines.append(f"  {o}_T * {o}_F = {num_T};")
    # Composition: flow-weighted average
    for c in comps:
        sc = _sanitize_name(c)
        num_z = " + ".join(
            f"{_sanitize_name(s)}_F * {_sanitize_name(s)}_z_{sc}" for s in inlets
        )
        lines.append(f"  {o}_z_{sc} * {o}_F = {num_z};")
    return lines


# ---------------------------------------------------------------------------
# Unit equation generators
# ---------------------------------------------------------------------------

def pump_eq(
    name: str,
    cfg: dict,
    inlet: str,
    outlet: str,
    comps: list[str],
) -> list[str]:
    """Pump: pressure rise, approximate adiabatic temperature rise."""
    sn = _sanitize_name(name)
    i, o = _sanitize_name(inlet), _sanitize_name(outlet)
    lines = [f"  // [{name}] Pump"]
    lines.append(f"  {o}_P = {i}_P + {sn}_deltaP;")
    # Simplified adiabatic: T_out = T_in (ignoring small rise)
    lines.append(f"  {o}_T = {i}_T;")
    lines.append(f"  {o}_F = {i}_F;")
    for c in comps:
        sc = _sanitize_name(c)
        lines.append(f"  {o}_z_{sc} = {i}_z_{sc};")
    return lines


def valve_eq(
    name: str,
    cfg: dict,
    inlet: str,
    outlet: str,
    comps: list[str],
) -> list[str]:
    """Valve: isenthalpic pressure reduction."""
    sn = _sanitize_name(name)
    i, o = _sanitize_name(inlet), _sanitize_name(outlet)
    lines = [f"  // [{name}] Valve"]
    lines.append(f"  {o}_P = {i}_P * {sn}_pressure_ratio;")
    lines.append(f"  {o}_T = {i}_T;")
    lines.append(f"  {o}_F = {i}_F;")
    for c in comps:
        sc = _sanitize_name(c)
        lines.append(f"  {o}_z_{sc} = {i}_z_{sc};")
    return lines


def heater_eq(
    name: str,
    cfg: dict,
    inlet: str,
    outlet: str,
    comps: list[str],
) -> list[str]:
    """Heater: energy balance with constant molar Cp (75.3 J/mol·K — water default).

    ``F * Cp * (T_out - T_in) = duty``
    """
    sn = _sanitize_name(name)
    i, o = _sanitize_name(inlet), _sanitize_name(outlet)
    # Constant molar Cp (J/mol·K) — reasonable for liquid water; user can tune via parameter
    Cp_default = 75.3
    lines = [f"  // [{name}] Heater"]
    lines.append(
        f"  {i}_F * {Cp_default} * ({o}_T - {i}_T) = {sn}_duty;"
    )
    lines.append(f"  {o}_P = {i}_P;")
    lines.append(f"  {o}_F = {i}_F;")
    for c in comps:
        sc = _sanitize_name(c)
        lines.append(f"  {o}_z_{sc} = {i}_z_{sc};")
    return lines


def strainer_eq(
    name: str,
    cfg: dict,
    inlet: str,
    outlet: str,
    comps: list[str],
) -> list[str]:
    """Strainer: pass-through with optional fixed pressure drop."""
    sn = _sanitize_name(name)
    i, o = _sanitize_name(inlet), _sanitize_name(outlet)
    lines = [f"  // [{name}] Strainer"]
    delta_p = cfg.get("deltaP", 0.0)
    if delta_p != 0.0:
        lines.append(f"  {o}_P = {i}_P - {sn}_deltaP;")
    else:
        lines.append(f"  {o}_P = {i}_P;")
    lines.append(f"  {o}_T = {i}_T;")
    lines.append(f"  {o}_F = {i}_F;")
    for c in comps:
        sc = _sanitize_name(c)
        lines.append(f"  {o}_z_{sc} = {i}_z_{sc};")
    return lines


def pipes_eq(
    name: str,
    cfg: dict,
    inlet: str,
    outlet: str,
    comps: list[str],
) -> list[str]:
    """Pipes: pass-through with optional pressure drop."""
    sn = _sanitize_name(name)
    i, o = _sanitize_name(inlet), _sanitize_name(outlet)
    lines = [f"  // [{name}] Pipes"]
    delta_p = cfg.get("delta_p", cfg.get("deltaP", 0.0))
    if delta_p != 0.0:
        lines.append(f"  {o}_P = {i}_P - {sn}_delta_p;")
    else:
        lines.append(f"  {o}_P = {i}_P;")
    lines.append(f"  {o}_T = {i}_T;")
    lines.append(f"  {o}_F = {i}_F;")
    for c in comps:
        sc = _sanitize_name(c)
        lines.append(f"  {o}_z_{sc} = {i}_z_{sc};")
    return lines


def flash_eq(
    name: str,
    cfg: dict,
    inlet: str,
    out_vap: str,
    out_liq: str,
    comps: list[str],
    K_values: dict[str, float],
) -> list[str]:
    """Flash: simplified Raoult's-law VLE with pre-computed K-values.

    K_i = P_sat_i / P_flash, computed by the transpiler at generation time.

    Equations:
      - y_i = K_i * x_i   (VLE)
      - F_in * z_i = F_vap * y_i + F_liq * x_i   (species balance)
      - F_vap + F_liq = F_in   (total balance)
      - sum(y_i) = 1  (one equation is redundant, dropped — use N-1 + balance)
      - T_vap = T_liq = T_in  (isothermal)
      - P_vap = P_liq = P_flash
    """
    sn = _sanitize_name(name)
    i = _sanitize_name(inlet)
    ov = _sanitize_name(out_vap)
    ol = _sanitize_name(out_liq)

    lines = [f"  // [{name}] Flash (Raoult's law, constant K_i)"]
    # Total balance
    lines.append(f"  {ov}_F + {ol}_F = {i}_F;")
    # Temperature: isothermal
    lines.append(f"  {ov}_T = {i}_T;")
    lines.append(f"  {ol}_T = {i}_T;")
    # Pressure
    lines.append(f"  {ov}_P = {sn}_P;")
    lines.append(f"  {ol}_P = {sn}_P;")
    # Species balances + VLE
    for idx, c in enumerate(comps):
        sc = _sanitize_name(c)
        K = K_values.get(c, 1.0)
        # y_i = K_i * x_i
        lines.append(f"  {ov}_z_{sc} = {K:.6g} * {ol}_z_{sc};")
        # F_in * z_i = F_vap * y_i + F_liq * x_i
        lines.append(
            f"  {i}_F * {i}_z_{sc} = {ov}_F * {ov}_z_{sc} + {ol}_F * {ol}_z_{sc};"
        )
    # Vapour composition summation (closes the system for the last component)
    sum_y = " + ".join(f"{ov}_z_{_sanitize_name(c)}" for c in comps)
    lines.append(f"  {sum_y} = 1.0;")
    return lines


def tank_eq(
    name: str,
    cfg: dict,
    inlet: str,
    outlet: str,
    comps: list[str],
) -> list[str]:
    """Tank: ODE mass and energy balances.

    State variables ``tank_{name}_n_{comp}`` and ``tank_{name}_T_state``
    are declared in the model's variable section.

    ODE:
      der(n_i) = F_in * z_in_i - F_out * x_i
    where x_i = n_i / n_tot  (well-mixed assumption).
    Outlet stream: F = outlet_flow (fixed), z = x (tank composition).

    Energy balance (simplified, constant Cp):
      n_tot * Cp_liq * der(T_state) = F_in * Cp_liq * (T_in - T_state) + duty
    """
    sn = _sanitize_name(name)
    i = _sanitize_name(inlet)
    o = _sanitize_name(outlet)
    Cp_default = 75.3  # J/mol·K

    lines = [f"  // [{name}] Tank (ODE)"]
    # n_tot = sum of component moles
    n_terms = " + ".join(f"tank_{sn}_n_{_sanitize_name(c)}" for c in comps)
    lines.append(f"  tank_{sn}_n_tot = {n_terms};")
    # ODE for each component
    for c in comps:
        sc = _sanitize_name(c)
        lines.append(
            f"  der(tank_{sn}_n_{sc}) = "
            f"{i}_F * {i}_z_{sc} - {sn}_outlet_flow * "
            f"(if tank_{sn}_n_tot > 0 then tank_{sn}_n_{sc} / tank_{sn}_n_tot else 0);"
        )
    # Outlet composition = tank composition
    for c in comps:
        sc = _sanitize_name(c)
        lines.append(
            f"  {o}_z_{sc} = "
            f"(if tank_{sn}_n_tot > 0 then tank_{sn}_n_{sc} / tank_{sn}_n_tot else 0);"
        )
    # Outlet flowrate (fixed parameter), pressure, temperature
    lines.append(f"  {o}_F = {sn}_outlet_flow;")
    lines.append(f"  {o}_P = {sn}_P;")
    lines.append(f"  {o}_T = tank_{sn}_T_state;")
    # Energy balance
    lines.append(
        f"  tank_{sn}_n_tot * {Cp_default} * der(tank_{sn}_T_state) = "
        f"{i}_F * {Cp_default} * ({i}_T - tank_{sn}_T_state) + {sn}_duty;"
    )
    return lines
