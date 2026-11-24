"""Shared hydraulic property resolution and unit math.

All fluid properties (Cp, molar volume, density, viscosity) are sourced from
the attached thermodynamic provider or from explicit flowsheet unit
parameters — never from hard-coded fluid defaults. Explicit unit parameters
always win over provider values so flowsheets can override any property.

The pressure–flow laws in this module are the single source of truth for the
native hydraulic network solver (:mod:`processforge.hydraulic_network`):

* Valves follow the IDAES ``Valve`` pressure-flow relation
  ``F**2 = Cv**2 * (P_in - P_out) * f(x)**2``.
* Pipes use Hagen–Poiseuille for laminar flow and Darcy–Weisbach with a
  Colebrook-approximate friction factor for turbulent flow.
* Strainers use a linear design-point resistance.
* Pumps support a fixed ``delta_p`` or an EPANET/WNTR-style three-point
  head-flow curve ``dP = A - B * Q**C``.
"""
from __future__ import annotations

import math
from typing import Callable, Optional

from loguru import logger


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class HydraulicError(Exception):
    """Base class for hydraulic modelling errors."""


class HydraulicConfigError(HydraulicError):
    """A unit or flowsheet is not usable in hydraulic network mode."""


class HydraulicSolveError(HydraulicError):
    """The hydraulic network could not be solved consistently."""


# ---------------------------------------------------------------------------
# Valve laws (IDAES-compatible)
# ---------------------------------------------------------------------------

_VALVE_FUNCTIONS = frozenset({"linear", "quick_opening", "equal_percentage"})


def valve_function_value(
    opening: float, valve_function: str = "linear", alpha: float = 100.0
) -> float:
    """Return ``f(x)`` for an IDAES-style valve.

    * ``linear``:            ``f(x) = x``
    * ``quick_opening``:     ``f(x) = sqrt(x)``
    * ``equal_percentage``:  ``f(x) = alpha**(x - 1)``

    Args:
        opening: Fraction open, from 0 to 1.
        valve_function: One of ``linear``, ``quick_opening`` or
            ``equal_percentage``.
        alpha: Equal-percentage valve parameter (default 100, matching
            IDAES).

    Returns:
        The valve function value, roughly in ``[0, 1]``.
    """
    if opening < 0.0 or opening > 1.0:
        raise HydraulicConfigError(
            f"valve_opening must be a fraction in [0, 1]; got {opening!r}"
        )
    if valve_function not in _VALVE_FUNCTIONS:
        raise HydraulicConfigError(
            f"valve_function must be one of {sorted(_VALVE_FUNCTIONS)}; "
            f"got {valve_function!r}"
        )
    if valve_function == "linear":
        return opening
    if valve_function == "quick_opening":
        return math.sqrt(opening) if opening > 0.0 else 0.0
    # equal_percentage
    if alpha <= 0.0:
        raise HydraulicConfigError(
            f"alpha must be positive for equal_percentage valves; got {alpha!r}"
        )
    return alpha ** (opening - 1.0)


def valve_delta_p(
    flow_mol: float,
    cv: float,
    opening: float = 1.0,
    valve_function: str = "linear",
    alpha: float = 100.0,
) -> float:
    """Signed pressure drop [Pa] across an IDAES-style valve.

    Inverts ``F**2 = Cv**2 * (P_in - P_out) * f(x)**2`` to

    ``P_in - P_out = (F / (Cv * f(x)))**2 * sign(F)``

    so reverse flow produces a negative pressure change (i.e. the outlet is at
    a higher pressure than the inlet).  A *positive* return value means the
    pressure drops from inlet to outlet.

    Args:
        flow_mol: Molar flow through the valve [mol/s]; negative values
            represent reverse flow.
        cv: IDAES valve flow coefficient [mol/s/Pa**0.5].
        opening: Fraction open, 0 to 1.
        valve_function: Valve characteristic — see
            :func:`valve_function_value`.
        alpha: Equal-percentage valve parameter.

    Returns:
        ``P_in - P_out`` [Pa].
    """
    if cv is None or cv <= 0.0:
        raise HydraulicConfigError(
            f"Cv must be positive for valves in hydraulic mode; got {cv!r}"
        )
    f = valve_function_value(opening, valve_function, alpha)
    if f <= 0.0:
        # A fully closed valve is handled by the network active set; treat it
        # as an infinite resistance so no flow is computed here.
        return math.inf if flow_mol != 0.0 else 0.0
    sign = 1.0 if flow_mol >= 0.0 else -1.0
    return (flow_mol / (cv * f)) ** 2 * sign


# ---------------------------------------------------------------------------
# Pipe and strainer laws
# ---------------------------------------------------------------------------

def laminar_pipe_flow(diameter, delta_p, viscosity, length) -> float:
    """Hagen–Poiseuille volumetric flow [m³/s] for a laminar pipe segment."""
    if diameter <= 0 or viscosity <= 0 or length <= 0 or delta_p <= 0:
        return 0.0
    return math.pi * diameter**4 * delta_p / (128.0 * viscosity * length)


def _friction_factor(re: float, roughness: float, diameter: float) -> float:
    """Darcy friction factor from the Swamee–Jain explicit Colebrook fit."""
    if re < 2000.0:
        return 64.0 / re if re > 0.0 else 64.0
    rel_rough = max(roughness / diameter, 0.0) if diameter > 0 else 0.0
    if rel_rough <= 0.0:
        # Smooth-pipe Swamee–Jain correlation
        return 0.25 / (math.log10(5.74 / re ** 0.9)) ** 2
    return 0.25 / (
        math.log10(rel_rough / 3.7 + 5.74 / re ** 0.9)
    ) ** 2


def pipe_delta_p(
    flow_mol: float,
    diameter: float,
    length: float,
    viscosity: float,
    molar_volume: float,
    density: float,
    roughness: float = 0.0,
) -> float:
    """Signed pressure drop [Pa] across a pipe segment.

    Laminar flow uses Hagen–Poiseuille; turbulent flow uses Darcy–Weisbach
    with a Swamee–Jain friction factor.

    Args:
        flow_mol: Molar flow [mol/s]; negative values represent reverse flow.
        diameter: Pipe inside diameter [m].
        length: Pipe length [m].
        viscosity: Dynamic viscosity [Pa·s].
        molar_volume: Fluid molar volume [m³/mol].
        density: Fluid density [kg/m³] (needed for Reynolds number).
        roughness: Absolute pipe roughness [m] (default smooth).

    Returns:
        ``P_in - P_out`` [Pa].
    """
    if diameter is None or diameter <= 0.0:
        raise HydraulicConfigError("Pipes in hydraulic mode need a positive diameter")
    if length is None or length <= 0.0:
        raise HydraulicConfigError("Pipes in hydraulic mode need a positive length")
    if viscosity is None or viscosity <= 0.0:
        raise HydraulicConfigError("Pipes in hydraulic mode need a positive viscosity")
    if molar_volume is None or molar_volume <= 0.0:
        raise HydraulicConfigError(
            "Pipes in hydraulic mode need a positive molar_volume"
        )
    if density is None or density <= 0.0:
        raise HydraulicConfigError("Pipes in hydraulic mode need a positive density")

    if flow_mol == 0.0:
        return 0.0

    q = flow_mol * molar_volume
    area = math.pi * diameter**2 / 4.0
    velocity = q / area
    re = abs(density * velocity * diameter / viscosity)
    sign = 1.0 if flow_mol >= 0.0 else -1.0

    if re < 2000.0:
        return sign * 128.0 * viscosity * length * abs(q) / (math.pi * diameter**4)
    f_d = _friction_factor(re, roughness, diameter)
    return sign * f_d * (length / diameter) * density * velocity**2 / 2.0


def design_point_delta_p(
    flow_mol: float, delta_p: float, design_flowrate: float
) -> float:
    """Linear design-point resistance for elements without full geometry.

    ``dP = delta_p_design * (F / F_design)``, so the element reproduces the
    specified pressure drop at the design flow and scales linearly away from
    it.  Negative flows produce a negative pressure change.
    """
    if delta_p is None or delta_p < 0.0:
        raise HydraulicConfigError(
            "delta_p must be non-negative for design-point hydraulic elements"
        )
    if design_flowrate is None or design_flowrate <= 0.0:
        raise HydraulicConfigError(
            "design_flowrate must be positive for design-point hydraulic elements"
        )
    return delta_p * flow_mol / design_flowrate


# ---------------------------------------------------------------------------
# Pump laws
# ---------------------------------------------------------------------------

_GRAVITY = 9.80665  # m/s²


def fit_pump_curve(points: list[tuple[float, float]]) -> tuple[float, float, float]:
    """Fit an EPANET/WNTR-style head curve to ``(flow, delta_p)`` points.

    Fits ``dP = A - B * Q**C``.  Following WNTR's ``HeadPump``:

    * 1 point:  ``A = 4/3 dP``, ``B = (1/3) dP / Q**2``, ``C = 2``.
    * 2 points: straight line (``C = 1``).
    * 3+ points: nonlinear fit via :func:`scipy.optimize.curve_fit`.

    Args:
        points: Sorted or unsorted list of ``(flow_mol, delta_p)`` pairs.
            Flow is molar flow [mol/s]; delta_p is pressure rise [Pa].

    Returns:
        ``(A, B, C)`` coefficients of ``dP = A - B * F**C``.
    """
    pts = sorted(points, key=lambda p: p[0])
    if not pts:
        raise HydraulicConfigError("pump_curve needs at least one (flow, delta_p) point")
    for f, dp in pts:
        if f < 0.0 or dp < 0.0:
            raise HydraulicConfigError(
                "pump_curve flow and delta_p/head values must be non-negative"
            )

    if len(pts) == 1:
        q0, h0 = pts[0]
        if q0 <= 0.0 or h0 <= 0.0:
            raise HydraulicConfigError(
                "A single pump-curve point needs a positive flow and head/delta_p"
            )
        return (4.0 / 3.0 * h0, (1.0 / 3.0) * h0 / q0**2, 2.0)

    if len(pts) == 2:
        (q0, h0), (q1, h1) = pts
        if q1 <= q0:
            raise HydraulicConfigError("pump_curve points must have distinct flows")
        # Straight line: dP = A - B*Q  (C = 1)
        b = (h0 - h1) / (q1 - q0)
        a = h0 - b * q0
        return (a, b, 1.0)

    # Three or more points: nonlinear power-law fit.
    import numpy as np
    from scipy.optimize import curve_fit

    q = np.array([p[0] for p in pts], dtype=float)
    h = np.array([p[1] for p in pts], dtype=float)

    def model(x, a, b, c):
        return a - b * x**c

    a0 = max(h.max(), 1.0)
    b0 = max((a0 - h.min()) / max(q.max(), 1e-12) ** 2, 1e-12)
    try:
        popt, _ = curve_fit(
            model,
            q,
            h,
            p0=[a0, b0, 2.0],
            bounds=([0.0, 0.0, 0.1], [np.inf, np.inf, 10.0]),
            maxfev=20000,
        )
        a, b, c = (float(v) for v in popt)
    except Exception as exc:  # noqa: BLE001 — fall back to a linear fit
        logger.warning(
            f"pump_curve power-law fit failed ({exc}); using a linear fit"
        )
        coeffs = np.polyfit(q, h, 1)  # h = m*q + k -> A - B*Q with C = 1
        a = float(coeffs[1])
        b = -float(coeffs[0])
        c = 1.0
    if b < 0.0:
        # Enforce a decreasing head curve; a negative slope is meaningless.
        raise HydraulicConfigError(
            "pump_curve fit produced a rising head curve (B < 0); check the points"
        )
    return (a, b, c)


def pump_curve_points_from_config(
    pump_curve: dict, density: Optional[float]
) -> list[tuple[float, float]]:
    """Convert a ``pump_curve`` config dict to ``(flow, delta_p)`` points.

    Accepts either a ``delta_p`` list (pressure rise in Pa) or a ``head``
    list (m).  ``head`` is converted using ``dP = rho * g * head``.
    """
    if not isinstance(pump_curve, dict):
        raise HydraulicConfigError("pump_curve must be an object")
    flows = pump_curve.get("flow")
    if flows is None:
        raise HydraulicConfigError("pump_curve requires a 'flow' array")
    if not isinstance(flows, list) or not flows:
        raise HydraulicConfigError("pump_curve 'flow' must be a non-empty array")

    if "delta_p" in pump_curve:
        dps = pump_curve["delta_p"]
        if not isinstance(dps, list) or len(dps) != len(flows):
            raise HydraulicConfigError(
                "pump_curve 'delta_p' must be an array the same length as 'flow'"
            )
        return [(float(f), float(dp)) for f, dp in zip(flows, dps)]

    if "head" in pump_curve:
        heads = pump_curve["head"]
        if not isinstance(heads, list) or len(heads) != len(flows):
            raise HydraulicConfigError(
                "pump_curve 'head' must be an array the same length as 'flow'"
            )
        if density is None or density <= 0.0:
            raise HydraulicConfigError(
                "A head-based pump_curve needs the fluid density "
                "(from the provider or an explicit 'density' parameter)"
            )
        return [
            (float(f), density * _GRAVITY * float(h))
            for f, h in zip(flows, heads)
        ]

    raise HydraulicConfigError(
        "pump_curve requires either a 'head' or a 'delta_p' array"
    )


def pump_pressure_rise(
    flow_mol: float,
    delta_p: Optional[float] = None,
    curve: Optional[tuple[float, float, float]] = None,
) -> float:
    """Pressure rise [Pa] produced by a pump.

    Args:
        flow_mol: Molar flow [mol/s].
        delta_p: Fixed pressure rise [Pa]; used when ``curve`` is None.
        curve: ``(A, B, C)`` coefficients of ``dP = A - B * F**C``.

    Returns:
        ``P_out - P_in`` [Pa].
    """
    if curve is None:
        if delta_p is None:
            raise HydraulicConfigError(
                "A pump in hydraulic mode needs either delta_p or pump_curve"
            )
        return delta_p
    a, b, c = curve
    return a - b * (abs(flow_mol)) ** c


# ---------------------------------------------------------------------------
# Property resolution for hydraulic elements
# ---------------------------------------------------------------------------

_HYDRAULIC_PROP_KEYS = ("density", "viscosity", "molar_volume", "molar_mass", "cp")


def resolve_hydraulic_props(unit, stream: dict) -> dict:
    """Resolve the hydraulic properties needed by the network solver.

    Provider values are used first; explicit unit parameters (``density``,
    ``viscosity``, ``molar_volume``, ``molar_mass``, ``cp``) override them.

    Returns:
        Dict with ``density`` [kg/m³], ``viscosity`` [Pa·s], ``molar_volume``
        [m³/mol], ``molar_mass`` [kg/mol] and ``cp`` [J/mol/K].  Entries that
        cannot be resolved are ``None``.
    """
    z = stream.get("z", {}) or {}
    t = stream.get("T", 298.15)
    p = stream.get("P", 101325.0)

    props: dict[str, float | None] = {k: None for k in _HYDRAULIC_PROP_KEYS}

    provider = getattr(unit, "_provider", None)
    if provider is not None:
        try:
            raw = provider.get_thermo_properties({"z": z, "T": t, "P": p}) or {}
        except Exception:  # noqa: BLE001 — fall through to explicit overrides
            raw = {}
        for key in _HYDRAULIC_PROP_KEYS:
            val = raw.get(key)
            if val is not None and (not isinstance(val, (int, float)) or math.isfinite(val)):
                props[key] = float(val)

    # Explicit unit parameters always win.
    params = getattr(unit, "params", None) or {}
    for key in _HYDRAULIC_PROP_KEYS:
        val = params.get(key)
        if val is None:
            val = getattr(unit, key, None)
        if val is not None and (not isinstance(val, (int, float)) or math.isfinite(val)):
            props[key] = float(val)

    # Derive density from molar mass / molar volume when possible.
    if props["density"] is None and props["molar_mass"] and props["molar_volume"]:
        if props["molar_volume"] > 0.0:
            props["density"] = props["molar_mass"] / props["molar_volume"]
    # Derive molar volume from density / molar mass when possible.
    if props["molar_volume"] is None and props["molar_mass"] and props["density"]:
        if props["density"] > 0.0:
            props["molar_volume"] = props["molar_mass"] / props["density"]

    return props


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
