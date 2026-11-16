import copy
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from dataclasses import dataclass
from .base import BaseUnitMixin
from .registry import register_unit


# Fallback molar heat capacity [J/mol/K] used when no thermo provider is
# available (or the provider returns a non-positive Cp).  Matches the
# ``Cp_default`` used by the Modelica writer (unit_equations.tank_eq).
# This is not silently used — it is only a last resort when neither the
# thermodynamic provider nor an explicit ``cp`` unit parameter is available.
TANK_CP_FALLBACK = 75.3

def _finite(value):
    """Coerce *value* to a finite float, or return ``None``."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def tank_energy_derivative(unit, z_in: dict, T_in: float, x: dict, T: float,
                           F_in: float, n_tot: float) -> float:
    """``dT/dt`` for a well-mixed tank at the tank's own pressure.

    The tank is a constant-pressure boundary: both the inlet and the outlet
    enthalpy are evaluated at ``unit.P`` so that a difference between the
    inlet stream pressure and the tank pressure cannot appear as a spurious
    ``v·ΔP`` temperature source/sink.

    Well-mixed enthalpy balance (constant tank pressure):

        d(n·H)/dt = F_in·H_in − F_out·H + Q
        =>  dT/dt = (F_in·(H_in − H_out) + Q) / (n_tot·Cp)

    Properties come from the unit's thermodynamic provider; if the provider
    is unavailable or returns unusable H/Cp, an explicit ``cp`` unit
    parameter (J/mol/K) is used with the constant-Cp form
    ``dT/dt = (F_in·(T_in − T) + Q/cp) / n_tot``. Raises when neither source
    is available — no hard-coded fluid default.

    Args:
        unit: Tank-like object with ``P``, ``duty``, ``cp`` (optional) and
            ``_get_thermo_properties``.
        z_in: Inlet mole fractions.
        T_in: Inlet temperature [K].
        x: Tank (well-mixed) mole fractions.
        T: Tank temperature [K].
        F_in: Inlet molar flow [mol/s].
        n_tot: Tank molar inventory [mol].

    Returns:
        Time derivative of the tank temperature [K/s].
    """
    Q = _finite(getattr(unit, "duty", 0.0)) or 0.0
    if n_tot <= 0.0:
        return 0.0

    import math
    print(f"DEBUG tank_energy_derivative: T_in={T_in:.15f}, T={T:.15f}, T_in==T: {T_in==T}, math.isclose: {math.isclose(T_in, T, rel_tol=1e-12, abs_tol=1e-15)}, Q={Q:.15f}")
    print(f"DEBUG tank_energy_derivative: about to check provider path")

    # If inlet temperature equals tank temperature and there is no heat duty,
    # there is no driving force for temperature change.
    if math.isclose(T_in, T, rel_tol=1e-12, abs_tol=1e-15) and Q == 0.0:
        print(f"DEBUG tank_energy_derivative: returning 0.0 (edge case)")
        return 0.0

    provider = getattr(unit, "_provider", None)
    if provider is not None:
        try:
            h_in = _finite(unit._get_thermo_properties(z_in, T_in, unit.P).get("H"))
            props_out = unit._get_thermo_properties(x, T, unit.P)
            h_out = _finite(props_out.get("H"))
            cp = _finite(props_out.get("Cp"))
            if h_in is not None and h_out is not None and cp and cp > 0:
                return (F_in * (h_in - h_out) + Q) / (n_tot * cp)
        except Exception:  # noqa: BLE001 — fall through to the constant-Cp form
            pass

    cp = _finite(getattr(unit, "cp", None))
    if cp and cp > 0:
        return (F_in * (T_in - T) + Q / cp) / n_tot
    raise RuntimeError(
        f"Tank '{getattr(unit, 'name', '?')}' cannot evaluate its energy "
        "balance: the thermodynamic provider did not return usable H/Cp and "
        "no explicit 'cp' unit parameter (J/mol/K) is set. Attach a provider "
        "(e.g. coolprop) or add \"cp\": <molar heat capacity> to the tank unit."
    )


def _series_interpolator(values, times: np.ndarray, default: float):
    """Return ``f(t)`` interpolating a (possibly constant) sample series.

    Falls back to a constant when the series is missing, scalar, or its
    length cannot be aligned with *times*.
    """
    if values is None:
        return lambda t: default
    if isinstance(values, (int, float)):
        return lambda t: float(values)
    try:
        arr = np.asarray([float(v) for v in values], dtype=float)
    except (TypeError, ValueError):
        return lambda t: default
    if arr.size == 0:
        return lambda t: default
    if arr.size == 1 or times.size < 2 or times[-1] <= times[0]:
        return lambda t: float(arr[0])
    xp = times if arr.size == times.size else np.linspace(times[0], times[-1], arr.size)
    if xp[-1] <= xp[0]:
        return lambda t: float(arr[0])
    return lambda t: float(np.interp(t, xp, arr))


@dataclass
class TankState:
    T: float
    P: float
    z: dict
    n_total: float


class Tank(BaseUnitMixin):
    """
    Simple well-mixed molar tank.
    Parameters expected in params:
      - vol: (optional) placeholder, not used explicitly here
      - inlet: name of inlet stream (must exist in flowsheet streams)
      - outlet_flow: molar flowrate out [mol/s] (constant)
      - initial_n: {comp: moles} initial contents (dict)
      - initial_T: initial temperature [K]
      - P: pressure in Pa (assumed constant; the energy balance evaluates
        both inlet and outlet enthalpy at this pressure)
      - duty: heat duty [W] applied to tank (positive heats)
      - flow_in: molar inlet flowrate [mol/s] (used if inlet stream not providing flowrate)
      - cp: (optional) explicit molar heat capacity [J/mol/K] used only when
        the thermodynamic provider cannot supply H/Cp
    """

    def __init__(self, name, **params):
        self.name = name
        self.params = params
        self.inlet = params.get("inlet")
        self.outflow = params.get("outlet_flow", 1.0)
        self.P = params.get("P", 101325)
        self.duty = params.get("duty", 0.0)
        # initial composition: dict of component -> moles
        self.initial_n = params.get("initial_n", {})
        self.initial_T = params.get("initial_T", 300.0)
        self.flow_in = params.get("flow_in", 1.0)  # mol/s if inlet doesn't specify
        # Explicit molar Cp override [J/mol/K] — provider properties always win.
        self.cp = params.get("cp")

    def _get_thermo_properties(self, z: dict, T: float, P: float) -> dict:
        """Route thermo calls through the attached provider when available."""
        provider = getattr(self, "_provider", None)
        if provider is None:
            raise RuntimeError(f"Unit '{self.name}' requires a thermodynamic provider to compute enthalpy and Cp.")
        return provider.get_thermo_properties({"z": z, "T": T, "P": P})

    def _run_impl(self, inlet):
        """
        Steady-state tank operation.
        The tank acts as a well-mixed vessel: the outlet takes the tank's
        flow rate and pressure, the composition passes through, and the
        temperature follows from an enthalpy balance at the tank pressure
        (``H(T_out, P_tank) = H(T_in, P_tank) + duty/F``).

        Args:
            inlet (dict): Inlet stream with T, P, z, flowrate

        Returns:
            dict: Outlet stream with tank's outlet flow rate and conditions
        """
        outlet = copy.deepcopy(inlet)

        # Tank outlet conditions
        outlet["flowrate"] = self.outflow
        outlet["P"] = self.P

        T_in = float(inlet.get("T", self.initial_T))
        z = inlet.get("z", {})
        F = self.outflow if self.outflow > 1e-12 else self.flow_in

        if self.duty and F > 1e-12:
            outlet["T"] = self._steady_outlet_temperature(z, T_in, float(self.duty) / F)
        else:
            outlet["T"] = T_in

        outlet["unit"] = "Tank"
        return outlet

    def _steady_outlet_temperature(self, z: dict, T_in: float, dh: float) -> float:
        """Solve ``H(T_out, P_tank) = H(T_in, P_tank) + dh`` for ``T_out``."""
        provider = getattr(self, "_provider", None)
        if provider is None:
            raise RuntimeError(
                f"Unit '{self.name}' requires a thermodynamic provider to compute enthalpy."
            )

        def H(T_val):
            return float(provider.get_thermo_properties(
                {"z": z, "T": T_val, "P": self.P}
            )["H"])

        h_target = H(T_in) + dh

        def residual(T_val):
            return H(T_val) - h_target

        lo, hi = (T_in, 2000.0) if dh > 0 else (100.0, T_in)
        try:
            return float(brentq(residual, lo, hi))
        except ValueError:
            # Provider enthalpy not bracket-solvable; fall back to a Cp-based
            # estimate from the provider Cp (or explicit cp override).
            cp = None
            try:
                cp = _finite(provider.get_thermo_properties(
                    {"z": z, "T": T_in, "P": self.P}
                ).get("Cp"))
            except Exception:  # noqa: BLE001
                pass
            if not (cp and cp > 0):
                cp = _finite(self.cp)
            if not (cp and cp > 0):
                raise RuntimeError(
                    f"Tank '{self.name}': cannot solve the steady-state energy "
                    "balance — the provider enthalpy is not bracket-solvable "
                    "and no Cp is available. Add \"cp\" to the tank unit."
                )
            return T_in + dh / cp

    def _ode_rhs(self, t, y, inlet_stream):
        """
        Computes the right-hand side of the ordinary differential equations (ODEs) for the tank unit.
        This method calculates the time derivatives of the molar amounts of species and the temperature
        based on inlet and outlet flows, enthalpy balances, and heat duty. It is used in numerical
        integration of the tank's dynamic behavior.
        Parameters
        ----------
        t : float
            Current time (in seconds). Not used in calculations but required for ODE solver interface.
        y : array-like
            State vector containing molar amounts of species (n) followed by temperature (T).
            Shape: (len(comps) + 1,), where comps are the components from inlet_stream.
        inlet_stream : dict
            Dictionary containing inlet stream properties:
            - "z": dict of mole fractions for each component.
            - "flowrate": inlet molar flowrate (mol/s), defaults to self.flow_in if not provided.
            - "T": inlet temperature (K).
        Returns
        -------
        np.ndarray
            Array of time derivatives: [dn_dt[0], dn_dt[1], ..., dn_dt[N], dT_dt],
            where N is the number of components. Shape: (len(comps) + 1,).
        Notes
        -----
        - Molar balances are computed for each species based on inlet and outlet flows.
        - The energy balance is the well-mixed enthalpy balance at the tank's own
          pressure (see :func:`tank_energy_derivative`); the inlet stream pressure
          is deliberately not used for enthalpy so pressure differences cannot
          create spurious temperature drift.
        - Assumes constant outlet flowrate (self.outflow) and tank pressure (self.P).
        """
        n = y[:-1]
        T = y[-1]
        comps = list(inlet_stream["z"].keys())
        n_tot = max(1e-12, sum(n))
        x = {comps[i]: max(0.0, n[i] / n_tot) for i in range(len(comps))}

        # inlet properties
        Fin = inlet_stream.get("flowrate", self.flow_in)  # mol/s
        z_in = inlet_stream["z"]
        T_in = inlet_stream["T"]

        # outlet flow (constant)
        Fout = self.outflow

        # species molar balances
        dn_dt = np.zeros(len(comps))
        for i, comp in enumerate(comps):
            dn_dt[i] = Fin * z_in.get(comp, 0.0) - Fout * x.get(comp, 0.0)

        # energy balance (molar basis) at the tank pressure
        dT_dt = tank_energy_derivative(self, z_in, T_in, x, T, Fin, n_tot)

        return np.concatenate([dn_dt, [dT_dt]])

    def _inlet_interpolator(self, inlet_stream_ts, times: np.ndarray):
        """Build a callable inlet snapshot ``f(t)`` from the inlet timeseries.

        ``T``, ``P`` and ``flowrate`` are linearly interpolated between
        samples; compositions are interpolated per component and
        renormalised. Missing or non-alignable series fall back to constants.
        """
        scalar_fns = {}
        for prop, default in (
            ("T", self.initial_T),
            ("P", self.P),
            ("flowrate", self.flow_in),
        ):
            scalar_fns[prop] = _series_interpolator(
                inlet_stream_ts.get(prop), times, default
            )
        z_fns = {
            comp: _series_interpolator(series, times, 0.0)
            for comp, series in (inlet_stream_ts.get("z") or {}).items()
        }

        def inlet_at(t):
            snap = {prop: fn(t) for prop, fn in scalar_fns.items()}
            z = {comp: max(0.0, fn(t)) for comp, fn in z_fns.items()}
            z_tot = sum(z.values())
            if z_tot > 0:
                z = {comp: val / z_tot for comp, val in z.items()}
            snap["z"] = z
            return snap

        return inlet_at

    def run_dynamic(self, inlet_stream_ts, t_span, t_eval, solver):
        """
        Integrate the tank ODE over the full inlet timeseries.

        The inlet conditions (T, P, flowrate, z) are interpolated at the
        solver's internal timesteps, so a downstream tank actually sees
        upstream transients instead of being frozen at the t=0 inlet values.

        Args:
            inlet_stream_ts: dict of timeseries data for the inlet stream
            t_span: (t0, tf)
            t_eval: list or array of times
            solver: solver object (currently unused but good for future)

        Returns:
            dict of timeseries data for the outlet stream
        """
        times = np.asarray(
            inlet_stream_ts.get("time") or t_eval, dtype=float
        )
        inlet_at = self._inlet_interpolator(inlet_stream_ts, times)

        comps = list((inlet_stream_ts.get("z") or {}).keys())

        # initial state vector
        n0 = [self.initial_n.get(c, 0.0) for c in comps]
        y0 = np.array(n0 + [self.initial_T])

        # integrate, sampling the (time-varying) inlet at each solver step
        sol = solve_ivp(
            lambda t, y: self._ode_rhs(t, y, inlet_at(t)),
            t_span,
            y0,
            t_eval=t_eval,
            rtol=1e-6,
            atol=1e-8,
        )

        # Reformat results into a timeseries dictionary for the outlet stream
        outlet_ts = {
            "time": sol.t.tolist(),
            "T": sol.y[-1, :].tolist(),
            "P": [self.P] * len(sol.t),
            "flowrate": [self.outflow] * len(sol.t),
        }
        # Add composition timeseries
        z_ts = {c: [] for c in comps}
        for i in range(len(sol.t)):
            n_vec = sol.y[:-1, i]
            n_tot = max(1e-12, sum(n_vec))
            for j, comp in enumerate(comps):
                z_ts[comp].append(n_vec[j] / n_tot)
        outlet_ts["z"] = z_ts

        return outlet_ts


register_unit("Tank", Tank)
