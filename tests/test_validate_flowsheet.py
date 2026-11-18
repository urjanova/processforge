"""Tests for flowsheet stream-connectivity validation."""
import pytest

from processforge.utils.validate_flowsheet import check_stream_connectivity


def _closed_loop(initial_n=True):
    """A feed-less loop: tank -> pipe -> pump -> pipe -> tank."""
    tank = {
        "type": "Tank",
        "in": "after_pipe_2",
        "out": "tank_out",
    }
    if initial_n:
        tank["initial_n"] = {"Water": 1.0}
    return {
        "streams": {},
        "units": {
            "tank_1": tank,
            "pipe_1": {
                "type": "Pipes",
                "in": "tank_out",
                "out": "after_pipe_1",
            },
            "pump_1": {
                "type": "Pump",
                "in": "after_pipe_1",
                "out": "after_pump",
            },
            "pipe_2": {
                "type": "Pipes",
                "in": "after_pump",
                "out": "after_pipe_2",
            },
        },
    }


def test_feedless_dynamic_closed_loop_is_connected():
    config = _closed_loop()
    config["simulation"] = {"mode": "dynamic"}

    # Returns None on success; raises ValueError on any failure.
    assert check_stream_connectivity(config) is None


def test_feedless_dynamic_loop_without_initial_n_is_unreachable():
    config = _closed_loop(initial_n=False)
    config["simulation"] = {"mode": "dynamic"}

    with pytest.raises(ValueError, match="Unreachable units detected"):
        check_stream_connectivity(config)


def test_feedless_steady_loop_is_unreachable():
    config = _closed_loop()
    config["simulation"] = {"mode": "steady"}

    with pytest.raises(ValueError, match="Unreachable units detected"):
        check_stream_connectivity(config)


def test_disconnected_island_is_still_reported():
    config = {
        "streams": {"feed": {}},
        "simulation": {"mode": "dynamic"},
        "units": {
            "pump_1": {"type": "Pump", "in": "feed", "out": "after_pump"},
            "pipe_1": {"type": "Pipes", "in": "loop_a", "out": "loop_b"},
            "pump_2": {"type": "Pump", "in": "loop_b", "out": "loop_a"},
        },
    }

    with pytest.raises(ValueError, match="Unreachable units detected") as excinfo:
        check_stream_connectivity(config)

    message = str(excinfo.value)
    assert "pipe_1" in message
    assert "pump_2" in message
    assert "pump_1" not in message


def test_feed_based_flowsheet_still_passes():
    config = {
        "streams": {"feed": {}},
        "simulation": {"mode": "steady"},
        "units": {
            "pump_1": {"type": "Pump", "in": "feed", "out": "after_pump"},
            "pipe_1": {"type": "Pipes", "in": "after_pump", "out": "after_pipe"},
            "valve_1": {"type": "Valve", "in": "after_pipe", "out": "product"},
        },
    }

    assert check_stream_connectivity(config) is None
