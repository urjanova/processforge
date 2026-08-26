"""Tests for StateManager drift detection."""
from processforge.state import StateManager
from processforge.types import SnapshotState


def _make_state(config):
    """Helper to create a minimal SnapshotState."""
    return SnapshotState(config=config, x=[], var_names=[])


# --- detect_drift: stream changes ---


def test_detect_drift_stream_T():
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 1.0}}, "units": {}}
    new = {"streams": {"feed": {"T": 350.0, "P": 101325, "flowrate": 1.0}}, "units": {}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == ["streams.feed.T"]


def test_detect_drift_stream_P():
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 1.0}}, "units": {}}
    new = {"streams": {"feed": {"T": 298.15, "P": 200000, "flowrate": 1.0}}, "units": {}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == ["streams.feed.P"]


def test_detect_drift_stream_flowrate():
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 1.0}}, "units": {}}
    new = {"streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 2.5}}, "units": {}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == ["streams.feed.flowrate"]


def test_detect_drift_stream_composition():
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 1.0, "z": {"H2O": 0.8, "TOL": 0.2}}}, "units": {}}
    new = {"streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 1.0, "z": {"H2O": 0.5, "TOL": 0.5}}}, "units": {}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert set(drifted) == {"streams.feed.z.H2O", "streams.feed.z.TOL"}


# --- detect_drift: unit changes ---


def test_detect_drift_unit_top_level_scalar():
    """Top-level unit scalar keys (like delta_p) should be detected."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {"pump_1": {"type": "Pump", "delta_p": 100000.0, "efficiency": 0.75}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "delta_p": 200000.0, "efficiency": 0.75}}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == ["units.pump_1.delta_p"]


def test_detect_drift_unit_multiple_scalars():
    """Multiple top-level scalar changes on the same unit."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {"pump_1": {"type": "Pump", "delta_p": 100000.0, "efficiency": 0.75}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "delta_p": 200000.0, "efficiency": 0.9}}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert set(drifted) == {"units.pump_1.delta_p", "units.pump_1.efficiency"}


def test_detect_drift_unit_nested_parameters():
    """Nested 'parameters' dict changes should still be detected."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {"pump_1": {"type": "Pump", "parameters": {"delta_p": 100000.0}}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "parameters": {"delta_p": 200000.0}}}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == ["units.pump_1.parameters.delta_p"]


def test_detect_drift_skips_non_scalar_keys():
    """Keys in the skip set (type, in, out, etc.) should NOT appear in drift paths."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {"pump_1": {"type": "Pump", "in": "feed", "out": "out1", "delta_p": 100000.0}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "in": "feed2", "out": "out2", "delta_p": 100000.0}}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert "units.pump_1.in" not in drifted
    assert "units.pump_1.out" not in drifted
    assert "units.pump_1.type" not in drifted
    assert drifted == []


def test_detect_drift_no_changes():
    """Identical configs produce no drift."""
    sm = StateManager("/tmp/_test_state.pfstate")
    config = {
        "streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 1.0}},
        "units": {"pump_1": {"type": "Pump", "delta_p": 100000.0, "efficiency": 0.75}},
    }
    drifted = sm.detect_drift(config, _make_state(config))
    assert drifted == []


# --- detect_drift: combined stream + unit changes ---


def test_detect_drift_combined():
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {
        "streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 1.0}},
        "units": {"pump_1": {"type": "Pump", "delta_p": 100000.0}},
    }
    new = {
        "streams": {"feed": {"T": 350.0, "P": 101325, "flowrate": 1.0}},
        "units": {"pump_1": {"type": "Pump", "delta_p": 200000.0}},
    }
    drifted = sm.detect_drift(new, _make_state(old))
    assert set(drifted) == {"streams.feed.T", "units.pump_1.delta_p"}


# --- detect_drift: dict-based state (not SnapshotState) ---


def test_detect_drift_dict_state():
    """detect_drift should work with dict-based state, not just SnapshotState."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"config": {"streams": {}, "units": {"pump_1": {"type": "Pump", "delta_p": 100000.0}}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "delta_p": 200000.0}}}
    drifted = sm.detect_drift(new, old)
    assert drifted == ["units.pump_1.delta_p"]


# --- detect_drift: other solve-affecting sections ---


def test_detect_drift_simulation_tf():
    """Changes to the simulation section (e.g. tf) must be detected."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {}, "simulation": {"mode": "steady", "t0": 0.0, "tf": 25.0, "dt": 1.0}}
    new = {"streams": {}, "units": {}, "simulation": {"mode": "steady", "t0": 0.0, "tf": 30.0, "dt": 1.0}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == ["simulation.tf"]


def test_detect_drift_materials():
    """Changes to materials (density, nested extra dicts) must be detected."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {
        "streams": {},
        "units": {},
        "materials": {
            "salt": {"id": 3, "density": 2.2, "temperature": 300.0},
            "tungsten": {"id": 1, "extra": {"D_0": 4.1e-7, "E_D": 0.39}},
        },
    }
    new = {
        "streams": {},
        "units": {},
        "materials": {
            "salt": {"id": 3, "density": 2.5, "temperature": 300.0},
            "tungsten": {"id": 1, "extra": {"D_0": 5.0e-7, "E_D": 0.39}},
        },
    }
    drifted = sm.detect_drift(new, _make_state(old))
    assert set(drifted) == {"materials.salt.density", "materials.tungsten.extra.D_0"}


def test_detect_drift_nested_unit_config():
    """Deeply nested unit solver_config must be detected (OpenMC/FESTIM style)."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {
        "streams": {},
        "units": {
            "openmc_solver": {
                "type": "SolverUnit",
                "provider": "openmc",
                "solver_config": {"batches": 20, "particles": 20000},
                "geometry_config": {"core_radius": 72.5},
            },
            "tds_solver": {
                "type": "SolverUnit",
                "provider": "festim",
                "solver_config": {"final_time": 500, "mesh": [{"start": 0.0}]},
            },
        },
    }
    new = {
        "streams": {},
        "units": {
            "openmc_solver": {
                "type": "SolverUnit",
                "provider": "openmc",
                "solver_config": {"batches": 40, "particles": 20000},
                "geometry_config": {"core_radius": 80.0},
            },
            "tds_solver": {
                "type": "SolverUnit",
                "provider": "festim",
                "solver_config": {"final_time": 600, "mesh": [{"start": 0.0}]},
            },
        },
    }
    drifted = sm.detect_drift(new, _make_state(old))
    assert set(drifted) == {
        "units.openmc_solver.solver_config.batches",
        "units.openmc_solver.geometry_config.core_radius",
        "units.tds_solver.solver_config.final_time",
    }


def test_detect_drift_skips_metadata():
    """Cosmetic metadata changes must NOT trigger drift."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {
        "streams": {},
        "units": {"pump_1": {"type": "Pump", "delta_p": 100000.0}},
        "metadata": {"name": "Old", "version": "1.0"},
    }
    new = {
        "streams": {},
        "units": {"pump_1": {"type": "Pump", "delta_p": 100000.0}},
        "metadata": {"name": "New Name", "version": "2.0"},
    }
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == []
