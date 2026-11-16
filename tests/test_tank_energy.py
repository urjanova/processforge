"""Tank energy balance tests.

Regression tests for the dynamic tank energy balance.  The old implementation
computed ``dT/dt = (F_in*T_in - F_out*T + duty/100) / n_tot``, which is only
correct when ``F_in == F_out``; with a recycle (``F_in > F_out``) it added a
spurious heating term that drove tank temperatures into the silent 500 K clamp.
The correct well-mixed, constant-Cp balance is::

    n_tot * Cp * dT/dt = F_in * Cp * (T_in - T) + Q
    dT/dt = (F_in * (T_in - T) + Q / Cp) / n_tot

matching the Modelica writer (``modelica/unit_equations.tank_eq``).

All fixtures use public flowsheets (``flowsheets/`` root and
``flowsheets/idaes/``) or in-test dict configs — nothing from
``flowsheets/custom/``.
"""
from __future__ import annotations

import json
import math
import sys
import types
from pathlib import Path

import pytest
from loguru import logger

from processforge.flowsheet import Flowsheet
from processforge.utils.validate_flowsheet import _validate_config_impl

TESTS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = TESTS_DIR.parent

CLOSED_LOOP = PROJECT_DIR / "flowsheets" / "closed-loop-chain.json"
CLOSED_LOOP_IDAES = PROJECT_DIR / "flowsheets" / "idaes" / "closed-loop-chain-idaes.json"

# Pump dT for the closed-loop chain: 200 kPa at 75% efficiency on 1 mol/s
# with the built-in water-property pump model -> ~0.016 K per pass.
EXPECTED_PUMP_DT = 0.044


def _load_flowsheet(path: Path) -> dict:
    """Load and validate a public flowsheet JSON file."""
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return _validate_config_impl(raw, source_name=str(path))


def _tank_energy_config(
    feed_T: float = 298.15,
    initial_T: float = 298.15,
    initial_n: float = 10.0,
    duty: float = 0.0,
    outlet_flow: float = 1.0,
) -> dict:
    """Minimal dynamic flowsheet: feed -> tank -> out (no recycle).

    A single tank with no cycle exercises the sequential dynamic path, which
    uses ``Tank.run_dynamic`` (the ``solve_ivp`` ODE with provider enthalpies).
    """
    return {
        "metadata": {"name": "tank-energy"},
        "providers": {"coolprop": {"type": "coolprop"}},
        "materials": {"Water": {"id": 1}, "Steel": {"id": 2}},
        "streams": {
            "feed": {"T": feed_T, "P": 101325.0, "flowrate": 1.0, "z": {"Water": 1.0}}
        },
        "units": {
            "tank_1": {
                "type": "Tank",
                "in": "feed",
                "out": "tank_out",
                "outlet_flow": outlet_flow,
                "initial_n": {"Water": initial_n},
                "initial_T": initial_T,
                "P": 101325,
                "duty": duty,
                "material": 1,
            }
        },
        "simulation": {"mode": "dynamic", "t0": 0.0, "tf": 5.0, "dt": 1.0},
    }


def _recycle_tank_config(
    feed_T: float = 298.15,
    initial_T: float = 298.15,
    initial_n: float = 10.0,
    duty: float = 0.0,
) -> dict:
    """Minimal dynamic flowsheet with a recycle loop through the tank.

    ``F_in = feed (1.0) + recycle (1.0) = 2.0`` while ``F_out = 1.0``, so
    ``F_in != F_out`` — the exact condition under which the old energy balance
    diverged.  Exercises ``Flowsheet._integrate_tank_step`` (Euler path).
    """
    return {
        "metadata": {"name": "tank-recycle-energy"},
        "providers": {"coolprop": {"type": "coolprop"}},
        "materials": {"Water": {"id": 1}, "Steel": {"id": 2}},
        "streams": {
            "feed": {"T": feed_T, "P": 101325.0, "flowrate": 1.0, "z": {"Water": 1.0}}
        },
        "units": {
            "tank_1": {
                "type": "Tank",
                "in": ["feed", "after_pipe"],
                "out": "tank_out",
                "outlet_flow": 1.0,
                "initial_n": {"Water": initial_n},
                "initial_T": initial_T,
                "P": 101325,
                "duty": duty,
                "material": 1,
            },
            "pipe_1": {
                "type": "Pipes",
                "in": "tank_out",
                "out": "after_pipe",
                "delta_p": 1000.0,
                "diameter": 0.05,
                "material": 2,
            },
        },
        "simulation": {"mode": "dynamic", "t0": 0.0, "tf": 5.0, "dt": 1.0},
    }


def _install_fake_idaes(monkeypatch):
    """Install a fake ``idaes`` package (mirrors test_idaes_provider.py)."""
    mock_idaes = types.ModuleType("idaes")
    monkeypatch.setitem(sys.modules, "idaes", mock_idaes)
    mock_core = types.ModuleType("idaes.core")
    mock_core.PureComponentPhaseBlock = object
    monkeypatch.setitem(sys.modules, "idaes.core", mock_core)
    for name in [
        "idaes.models",
        "idaes.models.properties",
        "idaes.models.properties.modular_properties",
    ]:
        if name not in sys.modules:
            monkeypatch.setitem(sys.modules, name, types.ModuleType(name))

    class _FakeInit:
        def get_k_value(self, comp, T, P):
            return 1.1

        def get_enthalpy(self, comp, T, P):
            return 1000.0

        def get_cp(self, comp, T, P):
            return 75.0

    monkeypatch.setitem(
        sys.modules,
        "idaes.models.properties.modular_properties",
        types.ModuleType("idaes.models.properties.modular_properties"),
    )
    sys.modules[
        "idaes.models.properties.modular_properties"
    ].ModularPropertiesInitializer = _FakeInit  # type: ignore[attr-defined]
    return mock_idaes


class _FakeTankUnit:
    """Minimal stand-in for a Tank unit, for direct helper tests."""

    def __init__(self, duty=0.0, outflow=1.0, P=101325.0, cp=75.3, fail=False):
        self.name = "fake_tank"
        self.duty = duty
        self.outflow = outflow
        self.P = P
        self._cp = cp
        self._fail = fail

    def _get_thermo_properties(self, z, T, P):
        if self._fail:
            raise RuntimeError("no provider attached")
        return {"H": 0.0, "Cp": self._cp}


# ---------------------------------------------------------------------------
# 1. Public closed-loop recycle flowsheet (CoolProp)
# ---------------------------------------------------------------------------


class TestClosedLoopRecycle:
    """``flowsheets/closed-loop-chain.json`` — dynamic recycle with a Tank."""

    @pytest.fixture
    def results(self):
        return Flowsheet(_load_flowsheet(CLOSED_LOOP)).run()

    def test_tank_holds_temperature(self, results):
        """No-duty tank stays within 0.5 K of initial_T over the whole run."""
        tank_T = results["after_tank_1"]["T"]
        assert tank_T[0] == pytest.approx(298.15)
        for t in tank_T:
            assert abs(t - 298.15) < 0.5

    def test_no_startup_dip(self, results):
        """Regression: t=0 used to plunge to 271 K because the empty recycle
        stream made F_out >> F_in in the old (wrong) energy balance."""
        assert results["after_tank_1"]["T"][0] == pytest.approx(298.15)

    def test_never_hits_temperature_clamp(self, results):
        """No stream may reach the 500 K clamp (or the 200 K floor)."""
        for name, stream in results.items():
            if not isinstance(stream, dict) or "T" not in stream:
                continue
            for t in stream["T"]:
                assert 200 < t < 500, f"{name} reached {t} K"

    def test_valve_and_pipe_outlets_match_inlets(self, results):
        """Valves and pipes are isothermal pass-throughs at every timestep."""
        pipe_pairs = [
            ("after_tank_1", "after_pipe_1"),
            ("after_valve_1", "after_pipe_2"),
            ("after_valve_2", "before_pump_1"),
            ("after_pump_1", "before_valve_3"),
            ("after_valve_3", "after_pipe_3"),
            ("after_valve_4", "recycle_to_tank"),
        ]
        valve_pairs = [
            ("after_pipe_1", "after_valve_1"),
            ("after_pipe_2", "after_valve_2"),
            ("before_valve_3", "after_valve_3"),
            ("after_pipe_3", "after_valve_4"),
        ]
        for upstream, downstream in pipe_pairs + valve_pairs:
            for t_in, t_out in zip(
                results[upstream]["T"], results[downstream]["T"]
            ):
                assert t_out == pytest.approx(t_in), (
                    f"{downstream} changed T relative to {upstream}"
                )

    def test_pump_adds_small_temperature_rise(self, results):
        """Positive control: the pump's adiabatic inefficiency still heats."""
        for t_in, t_out in zip(
            results["before_pump_1"]["T"], results["after_pump_1"]["T"]
        ):
            assert t_out - t_in == pytest.approx(EXPECTED_PUMP_DT, abs=0.005)


# ---------------------------------------------------------------------------
# 2. Same flowsheet with the IDAES provider on hydraulic units
# ---------------------------------------------------------------------------


class TestClosedLoopIdaes:
    """``flowsheets/idaes/closed-loop-chain-idaes.json`` — provider-independence.

    The tank itself has no ``provider`` key (IDAES does not model Tank), so it
    resolves to the built-in CoolProp fallback; the pipes/valves/pump use the
    (fake) IDAES provider.  The energy balance result must be identical.
    """

    def test_tank_holds_temperature_with_idaes(self, monkeypatch):
        _install_fake_idaes(monkeypatch)
        with open(CLOSED_LOOP_IDAES, "r", encoding="utf-8") as f:
            cfg = json.load(f)
        results = Flowsheet(cfg).run()

        tank_T = results["after_tank_1"]["T"]
        assert tank_T[0] == pytest.approx(298.15)
        for t in tank_T:
            assert abs(t - 298.15) < 0.5
            assert t < 400  # nowhere near the 500 K clamp


# ---------------------------------------------------------------------------
# 3. Sequential (non-recycle) dynamic path — Tank.run_dynamic ODE
# ---------------------------------------------------------------------------


class TestTankOdePath:
    """Single-tank flowsheet: exercises ``Tank._ode_rhs`` via ``run_dynamic``."""

    def test_isothermal_when_inlet_matches_initial(self):
        """T_in == T_initial, duty=0 -> T holds exactly at 298.15 K."""
        results = Flowsheet(_tank_energy_config()).run()
        for t in results["tank_out"]["T"]:
            assert t == pytest.approx(298.15)

    def test_relaxes_toward_inlet_temperature(self):
        """T_in=313 K warms the tank monotonically toward (not past) T_in.

        Analytic check: n=10 mol, F=1 mol/s, time constant n/F = 10 s, so
        T(5 s) = 298.15 + (313-298.15)*(1 - exp(-0.5)) ~= 304.0 K.
        """
        results = Flowsheet(_tank_energy_config(feed_T=313.0)).run()
        T = [float(t) for t in results["tank_out"]["T"]]
        assert all(T[i] <= T[i + 1] + 1e-9 for i in range(len(T) - 1))
        assert all(t < 313.0 for t in T)
        assert T[-1] == pytest.approx(304.0, abs=0.5)

    def test_duty_heats_at_expected_rate(self):
        """duty=753 W on n=10 mol water: steady state T_in + Q/(F*Cp) = 308.15 K.

        After 5 s (time constant 10 s): T ~ 298.15 + 10*(1-exp(-0.5)) ~ 302.1 K.
        """
        results = Flowsheet(_tank_energy_config(duty=753.0)).run()
        T = [float(t) for t in results["tank_out"]["T"]]
        rise = T[-1] - T[0]
        assert 3.0 < rise < 8.0
        assert T[-1] == pytest.approx(302.1, abs=0.5)

    def test_ode_holds_temperature_with_unbalanced_flows(self):
        """F_in != F_out must not create a spurious T term in the ODE either.

        outlet_flow (2.0) != feed (1.0): the tank drains, but with T_in == T
        and duty=0 the temperature must still hold at 298.15 K.
        """
        results = Flowsheet(
            _tank_energy_config(outlet_flow=2.0, initial_n=100.0)
        ).run()
        for t in results["tank_out"]["T"]:
            assert t == pytest.approx(298.15, abs=1e-6)


# ---------------------------------------------------------------------------
# 4. Recycle Euler path — Flowsheet._integrate_tank_step
# ---------------------------------------------------------------------------


class TestTankRecycleEulerPath:
    """Minimal recycle loop: F_in=2 vs F_out=1 — the divergence condition."""

    def test_holds_temperature_with_unbalanced_flows(self):
        """Regression: old balance gave (2*298.15 - 1*298.15)/10 ~ +30 K/s."""
        results = Flowsheet(_recycle_tank_config()).run()
        for t in results["tank_out"]["T"]:
            assert t == pytest.approx(298.15, abs=1e-9)

    def test_relaxes_toward_inlet_without_overshoot(self):
        """Warm feed: monotonic rise toward T_in, never above it, no clamp."""
        results = Flowsheet(_recycle_tank_config(feed_T=313.0)).run()
        T = [float(t) for t in results["tank_out"]["T"]]
        assert all(T[i] <= T[i + 1] + 1e-9 for i in range(len(T) - 1))
        assert all(t < 313.0 + 1e-9 for t in T)
        assert T[-1] > 299.0  # actually warming
        assert T[-1] < 400.0  # far from the clamp

    def test_duty_heats_recycle_tank(self):
        """duty=753 W, n=10 mol: ~1 K/s initial ramp, bounded relaxation."""
        results = Flowsheet(_recycle_tank_config(duty=753.0)).run()
        T = [float(t) for t in results["tank_out"]["T"]]
        assert T[-1] > T[0] + 1.0
        assert T[-1] < 400.0


# ---------------------------------------------------------------------------
# 5. Cp fallback and clamp warning behaviour
# ---------------------------------------------------------------------------


class TestTankCpHelper:
    def test_provider_cp_used(self):
        from processforge.flowsheet import Flowsheet as _F

        fs = _F.__new__(_F)  # bypass __init__; helper needs no flowsheet state
        unit = _FakeTankUnit(cp=146.0)
        assert fs._tank_cp(unit, {"Water": 1.0}, 298.15) == pytest.approx(146.0)

    def test_fallback_when_provider_raises(self):
        from processforge.flowsheet import Flowsheet as _F
        from processforge.units.tank import TANK_CP_FALLBACK

        fs = _F.__new__(_F)
        unit = _FakeTankUnit(fail=True)
        assert fs._tank_cp(unit, {"Water": 1.0}, 298.15) == TANK_CP_FALLBACK

    def test_fallback_when_provider_returns_zero_cp(self):
        """IDAES returns Cp=0.0 when its property package fails to load."""
        from processforge.flowsheet import Flowsheet as _F
        from processforge.units.tank import TANK_CP_FALLBACK

        fs = _F.__new__(_F)
        unit = _FakeTankUnit(cp=0.0)
        assert fs._tank_cp(unit, {"Water": 1.0}, 298.15) == TANK_CP_FALLBACK


class TestTankClampWarning:
    def test_out_of_bounds_temperature_is_clamped_and_warned(self):
        """A runaway energy balance clamps to 500 K and logs a warning.

        duty=1e6 W on 1 mol with Cp=75.3 gives dT/dt ~ 1.3e4 K/s — far past
        the bound.  The clamp must engage AND be visible in the logs (the old
        code clamped silently, which is how 500 K outputs went unnoticed).
        """
        records: list = []
        handler_id = logger.add(records.append, level="WARNING")
        try:
            from processforge.flowsheet import Flowsheet as _F

            fs = _F.__new__(_F)
            unit = _FakeTankUnit(duty=1_000_000.0, outflow=1.0)
            outlet, state = fs._integrate_tank_step(
                unit,
                {"T": 298.15, "P": 101325.0, "flowrate": 0.0, "z": {"Water": 1.0}},
                {"n": {"Water": 1.0}, "T": 298.15},
                dt=1.0,
                components={"Water"},
            )
        finally:
            logger.remove(handler_id)

        assert outlet["T"] == 500.0
        assert state["T"] == 500.0
        # loguru calls a plain-callable sink with the formatted message string
        assert any("clamping" in str(record) for record in records)
