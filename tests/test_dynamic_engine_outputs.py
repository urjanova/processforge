"""Dynamic flowsheets must collect ``SolverUnit`` ``EngineOutput`` records.

Steady-state runs already aggregated engine outputs; the dynamic paths
(``_run_dynamic_sequential`` / ``_run_dynamic_with_recycle``) previously
discarded them (and the sequential path raised ``KeyError`` on ``cfg["in"]``
for a unit with no inlet). These tests pin the fixed behaviour.
"""

import pytest

from processforge.flowsheet import Flowsheet
from processforge.providers.manager import ProviderMap
from processforge.types import (
    CoolPropProviderConfig,
    EngineOutput,
    OutputField,
    Quantity,
    _PROVIDER_CONFIG_REGISTRY,
)


class _FakeEngine:
    """Minimal provider returning a fixed ``EngineOutput`` once per call."""

    def __init__(self):
        self.calls = 0

    def initialize(self, cfg, flowsheet):
        pass

    def teardown(self):
        pass

    def compute_unit(self, *args, **kwargs):
        return None

    def run_simulation(self, unit_cfg, inlet=None):
        self.calls += 1
        return EngineOutput(
            engine="fake",
            sim_type="dyn",
            fields=[OutputField(name="power", quantity=Quantity(value=42.0, unit="W"))],
        )


def _install_fake_provider(monkeypatch):
    monkeypatch.setitem(_PROVIDER_CONFIG_REGISTRY, "fake", CoolPropProviderConfig)
    engine = _FakeEngine()
    pmap = ProviderMap(providers={"fake": engine}, default=None)
    monkeypatch.setattr(
        "processforge.flowsheet.build_provider_map", lambda *a, **k: pmap
    )
    monkeypatch.setattr(
        "processforge.flowsheet.teardown_providers", lambda *a, **k: None
    )
    return engine


def _sequential_config():
    return {
        "metadata": {"name": "dyn-engine"},
        "providers": {"fake": {"type": "fake"}},
        "materials": {"Water": {"id": 1}},
        "streams": {
            "feed": {
                "T": 300.0,
                "P": 101325.0,
                "flowrate": 1.0,
                "z": {"Water": 1.0},
            }
        },
        "units": {
            "engine": {
                "type": "SolverUnit",
                "provider": "fake",
                "sim_type": "dyn",
                "solver_config": {},
            }
        },
        "simulation": {"mode": "dynamic", "t0": 0.0, "tf": 2.0, "dt": 1.0},
    }


def _recycle_config():
    return {
        "metadata": {"name": "dyn-engine-recycle"},
        "providers": {"fake": {"type": "fake"}},
        "materials": {"Water": {"id": 1}, "Steel": {"id": 2}},
        "streams": {
            "feed": {
                "T": 300.0,
                "P": 101325.0,
                "flowrate": 1.0,
                "z": {"Water": 1.0},
            }
        },
        "units": {
            "tank_1": {
                "type": "Tank",
                "provider": "fake",
                "in": ["feed", "after_pipe"],
                "out": "tank_out",
                "outlet_flow": 1.0,
                "initial_n": {"Water": 1.0},
                "initial_T": 298.15,
                "P": 101325,
                "material": 1,
            },
            "pipe_1": {
                "type": "Pipes",
                "provider": "fake",
                "in": "tank_out",
                "out": "after_pipe",
                "delta_p": 1000.0,
                "diameter": 0.05,
                "material": 2,
            },
            "engine": {
                "type": "SolverUnit",
                "provider": "fake",
                "sim_type": "dyn",
                "solver_config": {},
            },
        },
        "simulation": {"mode": "dynamic", "t0": 0.0, "tf": 3.0, "dt": 1.0},
    }


def test_dynamic_sequential_collects_engine_output(monkeypatch):
    engine = _install_fake_provider(monkeypatch)

    fs = Flowsheet(_sequential_config())
    results = fs.run()

    assert isinstance(fs.engine_outputs.get("engine"), EngineOutput)
    assert fs.engine_outputs["engine"].field_value("power") == pytest.approx(42.0)
    assert isinstance(results["engine"], EngineOutput)
    # The standalone engine runs once, not once per timestep.
    assert engine.calls == 1

    manifest = fs.collect_outputs("20260101T000000Z_aaaaaa", "dynamic", "dyn-engine")
    assert "engine" in manifest.units


def test_dynamic_recycle_collects_engine_output(monkeypatch):
    engine = _install_fake_provider(monkeypatch)

    fs = Flowsheet(_recycle_config())
    results = fs.run()

    assert fs.has_recycles is True
    assert isinstance(fs.engine_outputs.get("engine"), EngineOutput)
    assert isinstance(results["engine"], EngineOutput)
    # Once, despite the timestep loop iterating over processing_order.
    assert engine.calls == 1

    manifest = fs.collect_outputs(
        "20260101T000000Z_bbbbbb", "dynamic", "dyn-engine-recycle"
    )
    assert "engine" in manifest.units
    # The engine output must not leak into the stream outputs.
    assert "engine" not in manifest.streams
