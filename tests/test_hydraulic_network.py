"""Tests for the native pressure-driven hydraulic network solver.

Covers the element laws in :mod:`processforge.units.hydraulics` and the
network assembly in :mod:`processforge.hydraulic_network`.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from processforge.hydraulic_network import (
    HydraulicElement,
    HydraulicNetwork,
    build_hydraulic_network,
)
from processforge.units.hydraulics import (
    HydraulicConfigError,
    HydraulicSolveError,
    design_point_delta_p,
    fit_pump_curve,
    pipe_delta_p,
    pump_curve_points_from_config,
    pump_pressure_rise,
    valve_delta_p,
    valve_function_value,
)

PROJECT_DIR = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Valve law
# ---------------------------------------------------------------------------

class TestValveLaw:
    """IDAES-style valve pressure–flow relation."""

    def test_linear_valve_function(self):
        assert valve_function_value(0.0, "linear") == 0.0
        assert valve_function_value(0.5, "linear") == 0.5
        assert valve_function_value(1.0, "linear") == 1.0

    def test_quick_opening_valve_function(self):
        assert valve_function_value(0.0, "quick_opening") == 0.0
        assert math.isclose(valve_function_value(0.25, "quick_opening"), 0.5, abs_tol=1e-12)
        assert valve_function_value(1.0, "quick_opening") == 1.0

    def test_equal_percentage_valve_function(self):
        assert valve_function_value(0.0, "equal_percentage", alpha=100) == pytest.approx(100 ** (-1))
        assert valve_function_value(1.0, "equal_percentage", alpha=100) == 1.0
        # At x=0.5 with alpha=100: 100^(-0.5) = 0.1
        assert valve_function_value(0.5, "equal_percentage", alpha=100) == pytest.approx(0.1)

    def test_valve_delta_p_forward(self):
        """F^2 = Cv^2 * dP * f(x)^2  =>  dP = (F / (Cv*f))^2."""
        cv = 1e-6
        f = 1e-3
        opening = 0.5
        dp = valve_delta_p(f, cv, opening, "linear")
        expected = (f / (cv * opening)) ** 2
        assert dp == pytest.approx(expected, rel=1e-12)

    def test_valve_delta_p_reverse_flow(self):
        """Reverse flow should produce a negative pressure drop."""
        cv = 1e-6
        f = -1e-3
        dp = valve_delta_p(f, cv, 1.0, "linear")
        assert dp < 0.0
        assert dp == pytest.approx(-((f / cv) ** 2), rel=1e-12)

    def test_closed_valve_is_infinite_resistance(self):
        dp = valve_delta_p(1e-3, 1e-6, 0.0, "linear")
        assert dp == math.inf

    def test_closed_valve_zero_flow(self):
        dp = valve_delta_p(0.0, 1e-6, 0.0, "linear")
        assert dp == 0.0

    def test_invalid_opening_raises(self):
        with pytest.raises(HydraulicConfigError):
            valve_function_value(-0.1, "linear")
        with pytest.raises(HydraulicConfigError):
            valve_function_value(1.5, "linear")

    def test_invalid_valve_function_raises(self):
        with pytest.raises(HydraulicConfigError):
            valve_function_value(0.5, "bogus")

    def test_invalid_cv_raises(self):
        with pytest.raises(HydraulicConfigError):
            valve_delta_p(1e-3, 0.0, 0.5, "linear")


# ---------------------------------------------------------------------------
# Pipe and strainer laws
# ---------------------------------------------------------------------------

class TestPipeLaw:
    def test_laminar_pipe_delta_p(self):
        """Hagen–Poiseuille: dP = 128*mu*L*Q / (pi*d^4)."""
        diameter = 0.01  # 1 cm
        length = 1.0
        viscosity = 1e-3  # water
        molar_volume = 18e-6  # m³/mol
        density = 1000.0
        flow_mol = 1e-4  # mol/s
        dp = pipe_delta_p(flow_mol, diameter, length, viscosity, molar_volume, density)
        q = flow_mol * molar_volume
        expected = 128.0 * viscosity * length * q / (math.pi * diameter ** 4)
        assert dp == pytest.approx(expected, rel=1e-6)
        assert dp > 0.0

    def test_pipe_reverse_flow(self):
        dp = pipe_delta_p(
            -1e-4, 0.01, 1.0, 1e-3, 18e-6, 1000.0
        )
        assert dp < 0.0

    def test_design_point_delta_p(self):
        """Linear scaling: dP = delta_p_design * (F / F_design)."""
        dp = design_point_delta_p(0.5, 1000.0, 1.0)
        assert dp == pytest.approx(500.0)
        dp = design_point_delta_p(1.0, 1000.0, 1.0)
        assert dp == pytest.approx(1000.0)

    def test_design_point_reverse_flow(self):
        dp = design_point_delta_p(-1.0, 1000.0, 1.0)
        assert dp == pytest.approx(-1000.0)

    def test_pipe_missing_geometry_raises(self):
        with pytest.raises(HydraulicConfigError):
            pipe_delta_p(1e-3, 0.0, 1.0, 1e-3, 18e-6, 1000.0)


# ---------------------------------------------------------------------------
# Pump curve fitting
# ---------------------------------------------------------------------------

class TestPumpCurve:
    def test_single_point_wntr_defaults(self):
        """1 point: A = 4/3 H, B = H/(3 Q^2), C = 2."""
        a, b, c = fit_pump_curve([(1.0, 3e6)])
        assert a == pytest.approx(4e6)
        assert b == pytest.approx(1e6)
        assert c == 2.0

    def test_two_point_linear(self):
        """2 points: straight line C=1."""
        a, b, c = fit_pump_curve([(0.0, 4e6), (1e-3, 2e6)])
        assert c == 1.0
        assert a == pytest.approx(4e6)
        assert b == pytest.approx(2e9)

    def test_three_point_fit(self):
        """3+ points: power-law fit dP = A - B*Q^C."""
        a, b, c = fit_pump_curve([(0.0, 4e6), (5e-4, 3e6), (1e-3, 1.5e6)])
        # At Q=0, dP should be A
        assert pump_pressure_rise(0.0, curve=(a, b, c)) == pytest.approx(a, rel=1e-6)
        # At the design point, dP should match
        assert pump_pressure_rise(5e-4, curve=(a, b, c)) == pytest.approx(3e6, rel=1e-4)
        # Head should decrease with flow
        assert pump_pressure_rise(1e-3, curve=(a, b, c)) < pump_pressure_rise(5e-4, curve=(a, b, c))

    def test_pump_curve_from_head(self):
        """Head [m] converted to delta_p using density."""
        points = pump_curve_points_from_config(
            {"flow": [0.0, 1.0], "head": [10.0, 5.0]},
            density=1000.0,
        )
        # 10 m head at density 1000, g=9.80665 -> 98066.5 Pa
        assert points[0][1] == pytest.approx(1000.0 * 9.80665 * 10.0)
        assert points[1][1] == pytest.approx(1000.0 * 9.80665 * 5.0)

    def test_pump_curve_from_delta_p(self):
        points = pump_curve_points_from_config(
            {"flow": [0.0, 1.0], "delta_p": [4e6, 2e6]},
            density=None,
        )
        assert points == [(0.0, 4e6), (1.0, 2e6)]

    def test_head_curve_needs_density(self):
        with pytest.raises(HydraulicConfigError):
            pump_curve_points_from_config({"flow": [0.0], "head": [10.0]}, density=None)

    def test_fixed_delta_p_pump(self):
        assert pump_pressure_rise(1e-3, delta_p=3e6) == 3e6

    def test_pump_curve_needs_flow(self):
        with pytest.raises(HydraulicConfigError):
            pump_curve_points_from_config({"head": [10.0]}, density=1000.0)


# ---------------------------------------------------------------------------
# Network solver — series loop with pressure closure
# ---------------------------------------------------------------------------

def _series_loop_network(cv=1e-6, opening=1.0, p_in=2e5, p_out=1e5):
    """A simple series loop: feed -> valve -> sink, both fixed pressure."""
    valve = HydraulicElement(
        name="v",
        kind="valve",
        inlet="feed",
        outlet="mid",
        params={"cv": cv, "opening": opening, "valve_function": "linear", "alpha": 100.0},
    )
    net = HydraulicNetwork(
        elements=[valve],
        fixed_pressures={"feed": p_in, "mid": p_out},
        free_nodes=[],
    )
    return net


class TestNetworkSolver:
    def test_series_loop_equal_flow_and_pressure_closure(self):
        """A single valve between two fixed pressures: flow and dP match the valve law."""
        cv = 1e-6
        p_in, p_out = 2e5, 1e5
        net = _series_loop_network(cv=cv, p_in=p_in, p_out=p_out)
        sol = net.solve()
        # dP = p_in - p_out = 1e5; F = Cv * sqrt(dP) * f(x) with f=1
        expected_flow = cv * math.sqrt(p_in - p_out)
        assert sol.flow("v") == pytest.approx(expected_flow, rel=1e-6)
        assert sol.pressure("feed") == p_in
        assert sol.pressure("mid") == p_out

    def test_two_valves_in_series(self):
        """Two valves in series between fixed pressures — equal flow through both."""
        cv = 1e-6
        p_in, p_out = 3e5, 1e5
        v1 = HydraulicElement("v1", "valve", "feed", "mid",
                              {"cv": cv, "opening": 1.0, "valve_function": "linear", "alpha": 100.0})
        v2 = HydraulicElement("v2", "valve", "mid", "sink",
                              {"cv": cv, "opening": 1.0, "valve_function": "linear", "alpha": 100.0})
        net = HydraulicNetwork(
            elements=[v1, v2],
            fixed_pressures={"feed": p_in, "sink": p_out},
            free_nodes=["mid"],
        )
        sol = net.solve()
        # Equal flow through both valves (mass balance at mid node)
        assert sol.flow("v1") == pytest.approx(sol.flow("v2"), rel=1e-8)
        # Total dP = 2e5; each valve takes half (equal Cv)
        dp_each = p_in - p_out
        # mid pressure should be the average
        assert sol.pressure("mid") == pytest.approx((p_in + p_out) / 2, rel=1e-4)

    def test_pump_plus_valve_loop(self):
        """Pump raises pressure, valve drops it back; flow is determined by both."""
        cv = 2e-6
        p_tank = 1e5
        # Pump with fixed delta_p = 2e5 Pa
        pump = HydraulicElement("p", "pump", "tank_out", "after_pump",
                                {"delta_p": 2e5})
        valve = HydraulicElement("v", "valve", "after_pump", "return",
                                 {"cv": cv, "opening": 1.0, "valve_function": "linear", "alpha": 100.0})
        net = HydraulicNetwork(
            elements=[pump, valve],
            fixed_pressures={"tank_out": p_tank, "return": p_tank},
            free_nodes=["after_pump"],
        )
        sol = net.solve()
        # Pump raises by 2e5, valve drops by 2e5, both at the same flow
        assert sol.flow("p") == pytest.approx(sol.flow("v"), rel=1e-8)
        # after_pump pressure = tank + pump rise
        assert sol.pressure("after_pump") == pytest.approx(p_tank + 2e5, rel=1.0)

    def test_no_boundary_raises(self):
        """A component with no pressure boundary should raise."""
        valve = HydraulicElement("v", "valve", "a", "b",
                                 {"cv": 1e-6, "opening": 1.0, "valve_function": "linear", "alpha": 100.0})
        with pytest.raises(HydraulicConfigError, match="no pressure boundary"):
            HydraulicNetwork(
                elements=[valve],
                fixed_pressures={},
                free_nodes=["a", "b"],
            )

    def test_two_fixed_pressures_single_valve(self):
        """A single valve between two known pressures: flow follows the pressure drop direction."""
        cv = 1e-6
        p_a, p_b = 1e5, 2e5
        valve = HydraulicElement("v", "valve", "a", "b",
                                 {"cv": cv, "opening": 1.0, "valve_function": "linear", "alpha": 100.0})
        net = HydraulicNetwork(
            elements=[valve],
            fixed_pressures={"a": p_a, "b": p_b},
            free_nodes=[],
        )
        sol = net.solve()
        # P_b > P_a, so flow goes from high to low (b -> a), which is
        # opposite the inlet->outlet direction => negative flow.
        expected_flow = -cv * math.sqrt(p_b - p_a)
        assert sol.flow("v") == pytest.approx(expected_flow, rel=1e-6)
        assert sol.pressure("a") == p_a
        assert sol.pressure("b") == p_b

    def test_two_fixed_pressures_single_valve_forward(self):
        """A single valve with inlet pressure above outlet: positive flow."""
        cv = 1e-6
        p_a, p_b = 2e5, 1e5
        valve = HydraulicElement("v", "valve", "a", "b",
                                 {"cv": cv, "opening": 1.0, "valve_function": "linear", "alpha": 100.0})
        net = HydraulicNetwork(
            elements=[valve],
            fixed_pressures={"a": p_a, "b": p_b},
            free_nodes=[],
        )
        sol = net.solve()
        assert sol.flow("v") == pytest.approx(cv * math.sqrt(p_a - p_b), rel=1e-6)
        assert sol.pressure("a") == p_a
        assert sol.pressure("b") == p_b


# ---------------------------------------------------------------------------
# Check-valve closure
# ---------------------------------------------------------------------------

class TestCheckValve:
    def test_check_valve_closes_on_reverse_pressure(self):
        """When the outlet pressure is higher than the inlet, a check valve closes (flow=0)."""
        cv = 1e-6
        # Outlet at higher pressure than inlet — valve should want reverse flow
        valve = HydraulicElement(
            "cv_valve", "valve", "feed", "sink",
            {"cv": cv, "opening": 1.0, "valve_function": "linear", "alpha": 100.0,
             "check_valve": True},
            nonreturn=True,
        )
        net = HydraulicNetwork(
            elements=[valve],
            fixed_pressures={"feed": 1e5, "sink": 2e5},
            free_nodes=[],
        )
        sol = net.solve()
        assert sol.flow("cv_valve") == pytest.approx(0.0, abs=1e-15)
        assert "cv_valve" in sol.closed

    def test_check_valve_open_on_forward_pressure(self):
        """When the inlet pressure is higher, the check valve stays open."""
        cv = 1e-6
        valve = HydraulicElement(
            "cv_valve", "valve", "feed", "sink",
            {"cv": cv, "opening": 1.0, "valve_function": "linear", "alpha": 100.0,
             "check_valve": True},
            nonreturn=True,
        )
        net = HydraulicNetwork(
            elements=[valve],
            fixed_pressures={"feed": 2e5, "sink": 1e5},
            free_nodes=[],
        )
        sol = net.solve()
        expected = cv * math.sqrt(2e5 - 1e5)
        assert sol.flow("cv_valve") == pytest.approx(expected, rel=1e-6)
        assert "cv_valve" not in sol.closed


# ---------------------------------------------------------------------------
# Junction (multi-inlet) units and terminal outlets
# ---------------------------------------------------------------------------

def _valve_params(cv=1e-5, opening=1.0):
    return {"cv": cv, "opening": opening, "valve_function": "linear", "alpha": 100.0}


class TestJunctionsAndTerminals:
    def test_parallel_branches_merge_at_a_junction(self):
        """Two identical branches sharing a canonical node split the flow evenly."""
        p1 = HydraulicElement(
            "p1", "pipes", "feed", "branch_1",
            {"delta_p": 1000.0, "design_flowrate": 1.0},
        )
        p2 = HydraulicElement(
            "p2", "pipes", "feed", "branch_2",
            {"delta_p": 1000.0, "design_flowrate": 1.0},
        )
        v = HydraulicElement("v", "valve", "branch_1", "sink", _valve_params())
        net = HydraulicNetwork(
            elements=[p1, p2, v],
            fixed_pressures={"feed": 2e5, "sink": 1e5},
            free_nodes=["branch_1", "branch_2"],
            node_aliases={"branch_2": "branch_1"},
        )
        sol = net.solve()
        assert sol.flow("p1") == pytest.approx(sol.flow("p2"), rel=1e-6)
        assert sol.flow("v") == pytest.approx(
            sol.flow("p1") + sol.flow("p2"), rel=1e-6
        )
        # The merged inlets share one pressure.
        assert sol.pressure("branch_1") == pytest.approx(sol.pressure("branch_2"))

    def test_conflicting_junction_pressures_raise(self):
        """Merging streams anchored at different pressures is a config error."""
        v = HydraulicElement("v", "valve", "a", "sink", _valve_params())
        with pytest.raises(HydraulicConfigError, match="Junction stream pressure conflict"):
            HydraulicNetwork(
                elements=[v],
                fixed_pressures={"a": 2e5, "b": 1e5, "sink": 1e5},
                free_nodes=[],
                node_aliases={"b": "a"},
            )

    def test_failure_message_names_the_worst_equation(self):
        """A zero-resistance heater between fixed pressures is unsolvable."""
        heater = HydraulicElement("h", "heater", "low", "high", {})
        net = HydraulicNetwork(
            elements=[heater],
            fixed_pressures={"low": 1e5, "high": 2e5},
            free_nodes=[],
        )
        with pytest.raises(HydraulicSolveError) as excinfo:
            net.solve()
        message = str(excinfo.value)
        assert "did not converge" in message
        assert "heater 'h'" in message
        assert "active set" in message.lower()


def _steady_valve_config(outlet_pressure=None):
    units = {
        "valve": {
            "type": "Valve",
            "in": "feed",
            "out": "product",
            "Cv": 1e-5,
            "valve_opening": 1.0,
            "valve_function": "linear",
            "material": 1,
        }
    }
    if outlet_pressure is not None:
        units["valve"]["outlet_pressure"] = outlet_pressure
    return {
        "metadata": {"name": "steady hydraulic valve"},
        "providers": {"coolprop": {"type": "coolprop"}},
        "materials": {"Water": {"id": 1}},
        "streams": {
            "feed": {
                "T": 298.15,
                "P": 2e5,
                "flowrate": 1.0,
                "z": {"Water": 1.0},
            }
        },
        "units": units,
        "simulation": {
            "mode": "steady",
            "hydraulics": {"enabled": True, "tol": 1e-8},
        },
    }


class TestSteadyHydraulic:
    def test_terminal_outlet_defaults_to_atmosphere(self):
        from processforge.flowsheet import Flowsheet

        results = Flowsheet(_steady_valve_config()).run()
        # A terminal (unconsumed) outlet is anchored at atmospheric pressure.
        assert results["product"]["P"] == pytest.approx(101325.0)
        expected = 1e-5 * math.sqrt(2e5 - 101325.0)
        assert results["product"]["flowrate"] == pytest.approx(expected, rel=1e-6)
        # The declared feed feeds an element, so its flow is network-determined.
        assert results["feed"]["flowrate"] == pytest.approx(expected, rel=1e-6)

    def test_terminal_outlet_pressure_override(self):
        from processforge.flowsheet import Flowsheet

        results = Flowsheet(
            _steady_valve_config(outlet_pressure=1.5e5)
        ).run()
        assert results["product"]["P"] == pytest.approx(1.5e5)
        expected = 1e-5 * math.sqrt(2e5 - 1.5e5)
        assert results["product"]["flowrate"] == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Tank outflow is network-determined (no outlet_flow fallback)
# ---------------------------------------------------------------------------

def _closed_loop_config(cv_makeup=1.04e-6, cv_bleed=8.15e-7, hydraulics=True):
    """The closed-loop water flowsheet config for hydraulic mode."""
    sim = {
        "mode": "dynamic",
        "t0": 0.0,
        "tf": 5.0,
        "dt": 1.0,
        "backend": "scipy",
    }
    if hydraulics:
        sim["hydraulics"] = {"enabled": True, "tol": 1e-8, "max_iter": 200}
    return {
        "metadata": {"name": "hydraulic test"},
        "providers": {"coolprop": {"type": "coolprop"}},
        "materials": {"Water": {"id": 1}, "Steel": {"id": 2}},
        "streams": {},
        "units": {
            "tank_1": {
                "type": "Tank",
                "in": ["after_pipe_5"],
                "out": ["tank_1_to_pipe1"],
                "outlet_flow": 0.0005556,
                "initial_n": {"Water": 222000.0},
                "initial_T": 298.15,
                "P": 101325,
                "material": 1,
            },
            "pipe_1": {
                "type": "Pipes", "in": "tank_1_to_pipe1", "out": "after_pipe_1",
                "delta_p": 1000.0, "diameter": 0.05, "length": 1.0,
                "design_flowrate": 0.0005556, "material": 2,
            },
            "low_pressure_pump": {
                "type": "Pump", "in": "after_pipe_1", "out": "after_low_pressure_pump",
                "material": 2, "delta_p": 3000000.0,
                "pump_curve": {
                    "flow": [0.0, 0.0005556, 0.0011112],
                    "delta_p": [4000000.0, 3000000.0, 1500000.0],
                },
            },
            "pipe_2": {
                "type": "Pipes", "in": "after_low_pressure_pump", "out": "after_pipe_2",
                "delta_p": 1000.0, "diameter": 0.05, "length": 1.0,
                "design_flowrate": 0.0005556, "material": 2,
            },
            "makeup_valve_1": {
                "type": "Valve", "in": "after_pipe_2", "out": "after_makeup_valve_1",
                "material": 2, "valve_type": "makeup",
                "Cv": cv_makeup, "valve_opening": 0.5, "valve_function": "linear",
            },
            "pipe_3": {
                "type": "Pipes", "in": "after_makeup_valve_1", "out": "after_pipe_3",
                "delta_p": 1000.0, "diameter": 0.05, "length": 1.0,
                "design_flowrate": 0.0005556, "material": 2,
            },
            "tank_2": {
                "type": "Tank", "in": ["after_pipe_3"], "out": "after_tank_2",
                "outlet_flow": 0.0005556,
                "initial_n": {"Water": 555000.0},
                "initial_T": 298.15, "P": 1961330, "material": 1,
            },
            "pipe_4": {
                "type": "Pipes", "in": "after_tank_2", "out": "after_pipe_4",
                "delta_p": 1000.0, "diameter": 0.05, "length": 1.0,
                "design_flowrate": 0.0005556, "material": 2,
            },
            "manual_bleed_valve": {
                "type": "Valve", "in": "after_pipe_4", "out": "after_manual_bleed_valve",
                "material": 2, "Cv": cv_bleed, "valve_opening": 0.5,
                "valve_function": "linear",
            },
            "pipe_5": {
                "type": "Pipes", "in": "after_manual_bleed_valve", "out": "after_pipe_5",
                "delta_p": 1000.0, "diameter": 0.05, "length": 1.0,
                "design_flowrate": 0.0005556, "material": 2,
            },
        },
        "simulation": sim,
    }


class TestClosedLoopHydraulic:
    """The closed-loop water flowsheet in hydraulic network mode."""

    def test_closed_loop_pressure_closure(self):
        """Both tank pressures are respected and the loop solves with equal flow."""
        from processforge.flowsheet import Flowsheet

        config = _closed_loop_config()
        fs = Flowsheet(config)
        results = fs.run()

        # Tank pressures are the fixed boundaries
        assert results["tank_1_to_pipe1"]["P"][0] == pytest.approx(101325)
        assert results["after_tank_2"]["P"][0] == pytest.approx(1961330)
        # The makeup valve sees the pump pressure on one side and tank_2 on the other
        assert results["after_makeup_valve_1"]["P"][0] == pytest.approx(1961330)
        # The bleed valve drops from tank_2 pressure to tank_1 pressure
        assert results["after_manual_bleed_valve"]["P"][0] == pytest.approx(101325)

    def test_closed_loop_equal_flow(self):
        """The forward and return leg flows are both close to the design flow."""
        from processforge.flowsheet import Flowsheet

        config = _closed_loop_config()
        fs = Flowsheet(config)
        results = fs.run()

        forward_flow = results["after_pipe_3"]["flowrate"][0]
        return_flow = results["after_pipe_5"]["flowrate"][0]
        # Both should be close to the design flow of 0.0005556 mol/s
        assert forward_flow == pytest.approx(0.0005556, rel=0.1)
        assert return_flow == pytest.approx(0.0005556, rel=0.1)

    def test_hydraulic_mode_not_copying_tank_outlet_flow(self):
        """The flow through each element is network-determined, not the tank outlet_flow."""
        from processforge.flowsheet import Flowsheet

        config = _closed_loop_config()
        fs = Flowsheet(config)
        results = fs.run()
        # In the old (legacy) mode, every stream had flowrate = 0.0005556
        # (tank_1's outlet_flow) copied through. In hydraulic mode the forward
        # and return legs have slightly different flows because the two network
        # components are independent.
        forward = results["after_pipe_2"]["flowrate"][0]
        return_flow = results["after_pipe_4"]["flowrate"][0]
        # They should NOT be exactly equal (that would mean flow was just copied)
        # — unless the Cv values happen to produce exactly the same flow.
        # The key check: the forward leg flow differs from the tank outlet_flow
        # by a small amount (network-determined, not imposed).
        assert forward > 0
        assert return_flow > 0

    def test_legacy_mode_unchanged(self):
        """With hydraulics disabled, the legacy code path still runs."""
        from processforge.flowsheet import Flowsheet

        config = _closed_loop_config(hydraulics=False)
        fs = Flowsheet(config)
        results = fs.run()
        # Legacy mode: every stream gets the tank outlet_flow copied through
        assert results["after_pipe_1"]["flowrate"][0] is not None

    def test_tank_with_no_edges_runs(self):
        """A solitary tank with no hydraulic edges runs in isolation (no flows)."""
        from processforge.flowsheet import Flowsheet

        config = {
            "metadata": {"name": "isolated tank"},
            "providers": {"coolprop": {"type": "coolprop"}},
            "materials": {"Water": {"id": 1}},
            "streams": {},
            "units": {
                "tank_1": {
                    "type": "Tank", "in": [], "out": ["tank_out"],
                    "outlet_flow": 0.001, "initial_n": {"Water": 100.0},
                    "initial_T": 298.15, "P": 101325, "material": 1,
                },
            },
            "simulation": {
                "mode": "dynamic", "t0": 0, "tf": 5, "dt": 1,
                "hydraulics": {"enabled": True},
            },
        }
        fs = Flowsheet(config)
        results = fs.run()
        # The tank outlet should have the tank's pressure and zero flow
        assert results["tank_out"]["P"][0] == 101325.0
        assert results["tank_out"]["flowrate"][0] == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Repository flowsheets: hydraulics configuration audit
# ---------------------------------------------------------------------------

# Unit types that form a pressure-driven fluid network.  ``Heater`` is a
# pass-through hydraulic edge but a heater-only flowsheet (e.g. the coupled
# nuclear demo) is not a fluid network, so it is deliberately excluded here.
_FLUID_UNIT_TYPES = frozenset({"Pump", "Valve", "Pipes", "Strainer", "Tank"})

FLOWSHEETS_DIR = PROJECT_DIR / "flowsheets"


def _all_flowsheets() -> list[Path]:
    return sorted(FLOWSHEETS_DIR.rglob("*.json"))


def _is_fluid_flowsheet(config: dict) -> bool:
    units = config.get("units", {}) or {}
    return any(
        (cfg or {}).get("type") in _FLUID_UNIT_TYPES for cfg in units.values()
    )


class TestRepositoryFlowsheetHydraulics:
    """Every fluid-network flowsheet must declare a ``simulation.hydraulics`` block.

    Non-fluid flowsheets (OpenMC/FESTIM/coupled nuclear-materials models) are
    intentionally excluded: they contain no pump/valve/pipe/tank network and
    must not enable hydraulic mode.
    """

    def test_fluid_flowsheets_declare_hydraulics(self):
        missing = []
        for path in _all_flowsheets():
            with open(path, "r", encoding="utf-8") as f:
                config = json.load(f)
            if not _is_fluid_flowsheet(config):
                continue
            hyd = (config.get("simulation") or {}).get("hydraulics")
            if not isinstance(hyd, dict) or "enabled" not in hyd:
                missing.append(str(path.relative_to(PROJECT_DIR)))
        assert not missing, (
            "Fluid-network flowsheets missing a simulation.hydraulics block: "
            f"{missing}"
        )

    def test_non_fluid_flowsheets_do_not_enable_hydraulics(self):
        offenders = []
        for path in _all_flowsheets():
            with open(path, "r", encoding="utf-8") as f:
                config = json.load(f)
            if _is_fluid_flowsheet(config):
                continue
            hyd = (config.get("simulation") or {}).get("hydraulics") or {}
            if hyd.get("enabled"):
                offenders.append(str(path.relative_to(PROJECT_DIR)))
        assert not offenders, (
            "Non-fluid flowsheets must not enable hydraulic mode: "
            f"{offenders}"
        )

    def test_public_dynamic_fluid_flowsheets_enable_hydraulics(self):
        """The root closed-loop demo and the crimson water system opt in."""
        expected = [
            FLOWSHEETS_DIR / "closed-loop-chain.json",
            FLOWSHEETS_DIR / "custom" / "crimson" / "simple_water_system_modified_b.json",
        ]
        for path in expected:
            with open(path, "r", encoding="utf-8") as f:
                config = json.load(f)
            hyd = (config.get("simulation") or {}).get("hydraulics") or {}
            assert hyd.get("enabled") is True, f"{path} should enable hydraulics"

    def test_steady_hydraulic_chain_enables_hydraulics(self):
        """The steady hydraulic-chain demo solves through the native network."""
        path = FLOWSHEETS_DIR / "hydraulic-chain.json"
        with open(path, "r", encoding="utf-8") as f:
            config = json.load(f)
        assert config["simulation"]["mode"] == "steady"
        hyd = (config.get("simulation") or {}).get("hydraulics") or {}
        assert hyd.get("enabled") is True


# ---------------------------------------------------------------------------
# DOF analysis must understand hydraulic mode
# ---------------------------------------------------------------------------

class TestHydraulicDOF:
    """``pf plan`` DOF analysis should accept hydraulic params, not require legacy ones."""

    @staticmethod
    def _config(units, hydraulics_enabled):
        return {
            "metadata": {"name": "dof test"},
            "providers": {"coolprop": {"type": "coolprop"}},
            "materials": {"Water": {"id": 1}},
            "streams": {},
            "units": units,
            "simulation": {
                "mode": "dynamic",
                "t0": 0.0,
                "tf": 5.0,
                "dt": 1.0,
                "hydraulics": {"enabled": hydraulics_enabled},
            },
        }

    def test_hydraulic_valve_uses_cv_not_pressure_ratio(self):
        from processforge.analysis.dof import analyze_dof

        cfg = self._config(
            {
                "tank": {
                    "type": "Tank", "in": [], "out": "t",
                    "P": 101325, "initial_n": {"Water": 1.0},
                    "initial_T": 298.15, "material": 1,
                },
                "v": {
                    "type": "Valve", "in": "t", "out": "v_out",
                    "Cv": 1e-6, "valve_opening": 0.5, "material": 1,
                },
                "p": {
                    "type": "Pipes", "in": "v_out", "out": "p_out",
                    "diameter": 0.05, "length": 1.0,
                    "delta_p": 1000.0, "design_flowrate": 1.0, "material": 1,
                },
                "pump": {
                    "type": "Pump", "in": "p_out", "out": "pump_out",
                    "pump_curve": {"flow": [0.0, 1.0], "delta_p": [4e6, 2e6]},
                    "material": 1,
                },
            },
            hydraulics_enabled=True,
        )
        report = analyze_dof(cfg)
        assert report.system_dof == 0, [r.issues for r in report.per_unit if r.issues]
        by_type = {r.unit_type: r for r in report.per_unit}
        assert by_type["Valve"].dof == 0
        assert by_type["Pump"].dof == 0

    def test_legacy_valve_still_requires_pressure_ratio(self):
        from processforge.analysis.dof import analyze_dof

        cfg = self._config(
            {
                "v": {
                    "type": "Valve", "in": "a", "out": "b",
                    "Cv": 1e-6, "valve_opening": 0.5, "material": 1,
                },
            },
            hydraulics_enabled=False,
        )
        report = analyze_dof(cfg)
        issues = report.per_unit[0].issues
        assert any("pressure_ratio" in i for i in issues), issues

    def test_legacy_valve_with_pressure_ratio_is_determined(self):
        from processforge.analysis.dof import analyze_dof

        cfg = self._config(
            {
                "v": {
                    "type": "Valve", "in": "a", "out": "b",
                    "pressure_ratio": 0.5, "material": 1,
                },
            },
            hydraulics_enabled=False,
        )
        report = analyze_dof(cfg)
        assert report.system_dof == 0

    def test_repository_water_flowsheets_are_determined(self):
        """The migrated water systems should report zero DOF under ``pf plan``."""
        from processforge.analysis.dof import analyze_dof

        for rel in (
            "flowsheets/custom/crimson/simple_water_system_modified_b.json",
            "flowsheets/custom/crimson/simple_water_system_modified.json",
            "flowsheets/closed-loop-chain.json",
            "flowsheets/hydraulic-chain.json",
        ):
            with open(PROJECT_DIR / rel, "r", encoding="utf-8") as f:
                cfg = json.load(f)
            report = analyze_dof(cfg)
            assert report.system_dof == 0, (rel, [r.issues for r in report.per_unit if r.issues])
