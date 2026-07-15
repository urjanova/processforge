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
    """Top-level unit scalar keys (like deltaP) should be detected."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {"pump_1": {"type": "Pump", "deltaP": 100000.0, "efficiency": 0.75}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "deltaP": 200000.0, "efficiency": 0.75}}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == ["units.pump_1.deltaP"]


def test_detect_drift_unit_multiple_scalars():
    """Multiple top-level scalar changes on the same unit."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {"pump_1": {"type": "Pump", "deltaP": 100000.0, "efficiency": 0.75}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "deltaP": 200000.0, "efficiency": 0.9}}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert set(drifted) == {"units.pump_1.deltaP", "units.pump_1.efficiency"}


def test_detect_drift_unit_nested_parameters():
    """Nested 'parameters' dict changes should still be detected."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {"pump_1": {"type": "Pump", "parameters": {"deltaP": 100000.0}}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "parameters": {"deltaP": 200000.0}}}}
    drifted = sm.detect_drift(new, _make_state(old))
    assert drifted == ["units.pump_1.parameters.deltaP"]


def test_detect_drift_skips_non_scalar_keys():
    """Keys in the skip set (type, in, out, etc.) should NOT appear in drift paths."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"streams": {}, "units": {"pump_1": {"type": "Pump", "in": "feed", "out": "out1", "deltaP": 100000.0}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "in": "feed2", "out": "out2", "deltaP": 100000.0}}}
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
        "units": {"pump_1": {"type": "Pump", "deltaP": 100000.0, "efficiency": 0.75}},
    }
    drifted = sm.detect_drift(config, _make_state(config))
    assert drifted == []


# --- detect_drift: combined stream + unit changes ---


def test_detect_drift_combined():
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {
        "streams": {"feed": {"T": 298.15, "P": 101325, "flowrate": 1.0}},
        "units": {"pump_1": {"type": "Pump", "deltaP": 100000.0}},
    }
    new = {
        "streams": {"feed": {"T": 350.0, "P": 101325, "flowrate": 1.0}},
        "units": {"pump_1": {"type": "Pump", "deltaP": 200000.0}},
    }
    drifted = sm.detect_drift(new, _make_state(old))
    assert set(drifted) == {"streams.feed.T", "units.pump_1.deltaP"}


# --- detect_drift: dict-based state (not SnapshotState) ---


def test_detect_drift_dict_state():
    """detect_drift should work with dict-based state, not just SnapshotState."""
    sm = StateManager("/tmp/_test_state.pfstate")
    old = {"config": {"streams": {}, "units": {"pump_1": {"type": "Pump", "deltaP": 100000.0}}}}
    new = {"streams": {}, "units": {"pump_1": {"type": "Pump", "deltaP": 200000.0}}}
    drifted = sm.detect_drift(new, old)
    assert drifted == ["units.pump_1.deltaP"]
