"""Tests for the multi-physics coupling resolver and flowsheet driver.

Covers:
* ``ParameterStore`` / ``CouplingResolver`` unit behaviour (reduce, unit
  conversion, deep-set/merge, error paths).
* End-to-end coupling through ``Flowsheet.run`` using two in-process fake
  providers, verifying that resolved + unit-converted values are injected and
  that acyclic and cyclic couplings converge.
"""

import pytest

from processforge.coupling import (
    CouplingError,
    CouplingResolver,
    ParameterStore,
    deep_merge,
    _deep_set,
    _reduce_value,
)
from processforge.flowsheet import Flowsheet
from processforge.providers.manager import ProviderMap
from processforge.types import EngineOutput, OutputField, Quantity


# ---------------------------------------------------------------------------
# Fake providers for the flowsheet integration tests
# ---------------------------------------------------------------------------


class _FakeA:
    """Produces a scalar ``power`` (W) and an array ``temperature`` (K)."""

    def __init__(self):
        self.last_cfg = None

    def initialize(self, cfg, flowsheet):
        pass

    def compute_unit(self, *a, **k):
        return None

    def run_simulation(self, unit_cfg, inlet=None):
        self.last_cfg = unit_cfg
        return EngineOutput(
            engine="fakeA",
            sim_type="a",
            fields=[
                OutputField(
                    name="power", quantity=Quantity(value=1000.0, unit="W")
                ),
                OutputField(
                    name="temperature",
                    quantity=Quantity(value=[800.0, 820.0, 840.0], unit="K"),
                ),
            ],
        )


class _FakeB:
    """Consumes coupled inputs; echoes them back for assertions."""

    def __init__(self):
        self.last_cfg = None

    def initialize(self, cfg, flowsheet):
        pass

    def compute_unit(self, *a, **k):
        return None

    def run_simulation(self, unit_cfg, inlet=None):
        self.last_cfg = unit_cfg
        hs = unit_cfg.solver_config.get("heat_source")
        t = unit_cfg.solver_config.get("temperature")
        return EngineOutput(
            engine="fakeB",
            sim_type="b",
            fields=[
                OutputField(name="heat_source_in", quantity=Quantity(value=float(hs), unit="W")),
                OutputField(name="temperature_in", quantity=Quantity(value=float(t), unit="K")),
                OutputField(name="result", quantity=Quantity(value=float(hs) * 2.0, unit="W")),
            ],
        )


def _make_provider_map():
    a, b = _FakeA(), _FakeB()
    pmap = ProviderMap(_providers={"fakeA": a, "fakeB": b}, _default=None)
    return pmap, a, b


# ---------------------------------------------------------------------------
# ParameterStore / CouplingResolver unit tests
# ---------------------------------------------------------------------------


def test_store_register_and_get():
    store = ParameterStore()
    store.register("u", "power", Quantity(value=5.0, unit="W"))
    assert "u.power" in store
    assert store.get("u.power").value == 5.0
    assert store.get("missing.field") is None


def test_resolver_reduce_mean_and_unit_convert():
    store = ParameterStore()
    store.register("u", "X", None)  # placeholder to keep structure
    store._quantities["u.temp"] = Quantity(value=[800.0, 820.0, 840.0], unit="K")
    store._quantities["u.power"] = Quantity(value=1000.0, unit="W")
    resolver = CouplingResolver(store)

    # mean of the array
    mean_val = resolver.resolve({"ref": "u.temp", "reduce": "mean"})
    assert mean_val == pytest.approx(820.0)

    # unit conversion W -> kW
    kw = resolver.resolve({"ref": "u.power", "as_unit": "kW"})
    assert kw == pytest.approx(1.0)

    # explicit unit passthrough
    w = resolver.resolve({"ref": "u.power", "as_unit": "W"})
    assert w == pytest.approx(1000.0)


def test_resolver_invalid_ref_raises():
    resolver = CouplingResolver(ParameterStore())
    with pytest.raises(CouplingError):
        resolver.resolve({"ref": "no_dot"})
    with pytest.raises(CouplingError):
        resolver.resolve({"ref": "missing.field"})


def test_resolver_missing_ref_raises():
    resolver = CouplingResolver(ParameterStore())
    with pytest.raises(CouplingError):
        resolver.resolve({"ref": "u.ghost"})


def test_deep_set_and_merge():
    d = {}
    _deep_set(d, "solver_config.heat_source", 3.0)
    assert d == {"solver_config": {"heat_source": 3.0}}

    base = {"solver_config": {"batches": 20, "particles": 1}}
    merged = deep_merge(base, {"solver_config": {"particles": 2}})
    assert merged["solver_config"]["batches"] == 20
    assert merged["solver_config"]["particles"] == 2
    # original untouched
    assert base["solver_config"]["particles"] == 1


def test_reduce_modes():
    assert _reduce_value([1.0, 3.0], "sum") == pytest.approx(4.0)
    assert _reduce_value([1.0, 3.0], "max") == pytest.approx(3.0)
    assert _reduce_value([1.0, 3.0], "min") == pytest.approx(1.0)
    assert _reduce_value(7.0, "mean") == pytest.approx(7.0)


# ---------------------------------------------------------------------------
# Flowsheet integration tests
# ---------------------------------------------------------------------------


def _coupled_config(cyclic=False):
    dst_inputs = {
        "solver_config.heat_source": {"ref": "src.power", "as_unit": "W"},
        "solver_config.temperature": {"ref": "src.temperature", "reduce": "mean", "as_unit": "K"},
    }
    src_inputs = {}
    if cyclic:
        # src consumes dst.result -> genuine cycle, must iterate to converge
        src_inputs = {"solver_config.feedback": {"ref": "dst.result", "as_unit": "W"}}

    return {
        "metadata": {"name": "coupling-test"},
        "providers": {"fakeA": {"type": "fakeA"}, "fakeB": {"type": "fakeB"}},
        "materials": {},
        "streams": {},
        "units": {
            "src": {
                "type": "SolverUnit",
                "provider": "fakeA",
                "sim_type": "a",
                "solver_config": {},
                "inputs": src_inputs,
            },
            "dst": {
                "type": "SolverUnit",
                "provider": "fakeB",
                "sim_type": "b",
                "solver_config": {},
                "inputs": dst_inputs,
            },
        },
        "simulation": {"mode": "steady"},
    }


def test_flowsheet_coupling_injects_resolved_values(monkeypatch):
    from processforge.types import (
        CoolPropProviderConfig,
        _PROVIDER_CONFIG_REGISTRY,
    )

    monkeypatch.setitem(_PROVIDER_CONFIG_REGISTRY, "fakeA", CoolPropProviderConfig)
    monkeypatch.setitem(_PROVIDER_CONFIG_REGISTRY, "fakeB", CoolPropProviderConfig)
    pmap, fake_a, fake_b = _make_provider_map()
    monkeypatch.setattr(
        "processforge.flowsheet.build_provider_map", lambda *a, **k: pmap
    )
    monkeypatch.setattr(
        "processforge.flowsheet.teardown_providers", lambda *a, **k: None
    )

    fs = Flowsheet(_coupled_config())
    results = fs.run()

    # The fake B must have received the unit-converted coupled values.
    sc = fake_b.last_cfg.solver_config
    assert sc["heat_source"] == pytest.approx(1000.0)  # W -> W
    assert sc["temperature"] == pytest.approx(820.0)   # mean of [800,820,840] K

    # And the resolved parameters also surface as engine outputs.
    assert fs.engine_outputs["dst"].field_value("heat_source_in") == pytest.approx(1000.0)
    assert fs.engine_outputs["dst"].field_value("temperature_in") == pytest.approx(820.0)
    assert results is not None


def test_flowsheet_cyclic_coupling_converges(monkeypatch):
    from processforge.types import (
        CoolPropProviderConfig,
        _PROVIDER_CONFIG_REGISTRY,
    )

    monkeypatch.setitem(_PROVIDER_CONFIG_REGISTRY, "fakeA", CoolPropProviderConfig)
    monkeypatch.setitem(_PROVIDER_CONFIG_REGISTRY, "fakeB", CoolPropProviderConfig)
    pmap, fake_a, fake_b = _make_provider_map()
    monkeypatch.setattr(
        "processforge.flowsheet.build_provider_map", lambda *a, **k: pmap
    )
    monkeypatch.setattr(
        "processforge.flowsheet.teardown_providers", lambda *a, **k: None
    )

    fs = Flowsheet(_coupled_config(cyclic=True))
    # Should not raise / infinite-loop; constant outputs converge quickly.
    fs.run()
    # Both units produced outputs and the coupling store resolved the cycle.
    assert "src" in fs.engine_outputs
    assert "dst" in fs.engine_outputs
