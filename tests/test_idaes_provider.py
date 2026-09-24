"""Tests for the IDAES provider — config, registry, schema and unit delegation.

IDAES is not installed in CI, so ``sys.modules['idaes']`` is mocked via
``monkeypatch`` where needed; the not-installed path is also exercised.
"""
from __future__ import annotations

import json
import sys
import types

import pytest

from processforge.providers.errors import ProviderNotAvailableError
from processforge.providers.idaes_provider import IdaesProvider
from processforge.types import FlowsheetConfig, IdaesProviderConfig, provider_config_from_dict


# ---------------------------------------------------------------------------
# Helpers — fake idaes package
# ---------------------------------------------------------------------------


def _install_fake_idaes(monkeypatch):
    mock_idaes = types.ModuleType("idaes")
    monkeypatch.setitem(sys.modules, "idaes", mock_idaes)
    # idaes.core.PureComponentPhaseBlock — imported in get_thermo_properties
    mock_core = types.ModuleType("idaes.core")
    mock_core.PureComponentPhaseBlock = object
    monkeypatch.setitem(sys.modules, "idaes.core", mock_core)
    # idaes.models.properties.modular_properties.ModularPropertiesInitializer
    for name in ["idaes.models", "idaes.models.properties", "idaes.models.properties.modular_properties"]:
        if name not in sys.modules:
            monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    class _FakeInit:
        def get_k_value(self, comp, T, P): return 1.1
        def get_enthalpy(self, comp, T, P): return 1000.0
        def get_cp(self, comp, T, P): return 75.0
    monkeypatch.setitem(
        sys.modules,
        "idaes.models.properties.modular_properties",
        types.ModuleType("idaes.models.properties.modular_properties"),
    )
    sys.modules["idaes.models.properties.modular_properties"].ModularPropertiesInitializer = _FakeInit  # type: ignore
    return mock_idaes


# ---------------------------------------------------------------------------
# 1. IdaesProviderConfig
# ---------------------------------------------------------------------------


class TestIdaesProviderConfig:
    def test_defaults(self):
        c = IdaesProviderConfig()
        assert c.type == "idaes"
        assert c.package == "ideal_pure_thermo"

    def test_from_dict(self):
        c = IdaesProviderConfig.from_dict({"package": "my_pkg"})
        assert c.package == "my_pkg"

    def test_from_dict_default(self):
        c = IdaesProviderConfig.from_dict({})
        assert c.package == "ideal_pure_thermo"

    def test_provider_config_from_dict_dispatch(self):
        c = provider_config_from_dict({"type": "idaes", "package": "foo"})
        assert isinstance(c, IdaesProviderConfig)
        assert c.package == "foo"


# ---------------------------------------------------------------------------
# 2. Registry / catalog
# ---------------------------------------------------------------------------


class TestIdaesRegistry:
    def test_catalog_entry(self):
        from processforge.providers.registry import _PROVIDER_CATALOG

        entry = _PROVIDER_CATALOG["idaes"]
        assert entry.module == "processforge.providers.idaes_provider"
        assert entry.class_name == "IdaesProvider"
        assert entry.optional_dep == "idaes"
        assert entry.docker_image is None

    def test_list_providers_includes_idaes(self):
        from processforge.providers.registry import list_providers

        types = {p["type"] for p in list_providers()}
        assert "idaes" in types

    def test_get_provider_class(self):
        from processforge.providers.registry import get_provider_class

        assert get_provider_class("idaes") is IdaesProvider

    def test_is_not_containerized(self):
        from processforge.providers.registry import is_containerized

        assert is_containerized("idaes") is False


# ---------------------------------------------------------------------------
# 3. Schema
# ---------------------------------------------------------------------------


class TestIdaesSchema:
    def test_provider_type_enum_contains_idaes(self):
        from processforge._schema import load_flowsheet_schema

        enum = load_flowsheet_schema()["properties"]["providers"]["patternProperties"]["^[a-zA-Z0-9_]+$"]["properties"]["type"]["enum"]
        assert "idaes" in enum


# ---------------------------------------------------------------------------
# 4. initialize
# ---------------------------------------------------------------------------


class TestIdaesInitialize:
    def test_raises_when_not_installed(self, monkeypatch):
        monkeypatch.delitem(sys.modules, "idaes", raising=False)
        # ensure submodules also removed
        for k in list(sys.modules):
            if k.startswith("idaes"):
                monkeypatch.delitem(sys.modules, k, raising=False)
        p = IdaesProvider()
        with pytest.raises(ProviderNotAvailableError, match="processforge\\[idaes\\]"):
            p.initialize(IdaesProviderConfig(), FlowsheetConfig())

    def test_succeeds_with_mock(self, monkeypatch):
        _install_fake_idaes(monkeypatch)
        p = IdaesProvider()
        p.initialize(IdaesProviderConfig(package="my_pkg"), FlowsheetConfig())
        assert p._package_name == "my_pkg"  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# 5. compute_unit — hydraulic units
# ---------------------------------------------------------------------------


class TestIdaesComputeUnit:
    @pytest.fixture(autouse=True)
    def _provider(self, monkeypatch):
        _install_fake_idaes(monkeypatch)
        self.p = IdaesProvider()
        self.p.initialize(IdaesProviderConfig(), FlowsheetConfig())

    def test_pump(self):
        inlet = {"T": 298.15, "P": 101325.0, "flowrate": 1.0}
        out = self.p.compute_unit("Pump", {"delta_p": 100_000, "efficiency": 0.75}, inlet)
        assert out["P"] == pytest.approx(201325.0)
        assert out["unit"] == "Pump"
        assert "power" in out

    def test_pump_uses_custom_delta_p(self):
        # regression: BaseUnitMixin must pass real delta_p, not default 1e5
        inlet = {"T": 298.15, "P": 101325.0, "flowrate": 1.0}
        out = self.p.compute_unit("Pump", {"delta_p": 200_000, "efficiency": 0.8}, inlet)
        assert out["P"] == pytest.approx(301325.0)

    def test_valve(self):
        inlet = {"T": 300.0, "P": 200_000.0}
        out = self.p.compute_unit("Valve", {"pressure_ratio": 0.5}, inlet)
        assert out["P"] == pytest.approx(100_000.0)
        assert out["unit"] == "Valve"

    def test_strainer(self):
        inlet = {"T": 300.0, "P": 100_000.0}
        out = self.p.compute_unit("Strainer", {"delta_p": 5000.0}, inlet)
        assert out["P"] == pytest.approx(95_000.0)

    def test_pipes(self):
        inlet = {"T": 300.0, "P": 10_000.0}
        out = self.p.compute_unit("Pipes", {"delta_p": 1000.0, "diameter": 0.1}, inlet)
        assert out["P"] == pytest.approx(9_000.0)
        assert out["unit"] == "Pipes"
        assert "flow" in out

    def test_heater_returns_none(self):
        assert self.p.compute_unit("Heater", {}, {"T": 300, "P": 1e5}) is None

    def test_flash_returns_none(self):
        assert self.p.compute_unit("Flash", {}, {"T": 300, "P": 1e5}) is None

    def test_unknown_returns_none(self):
        assert self.p.compute_unit("CSTR", {}, {"T": 300, "P": 1e5}) is None

    def test_valve_floor(self):
        inlet = {"T": 300.0, "P": 1500.0}
        out = self.p.compute_unit("Valve", {"pressure_ratio": 0.1}, inlet)
        assert out["P"] == pytest.approx(1000.0)


# ---------------------------------------------------------------------------
# 6. get_thermo_properties
# ---------------------------------------------------------------------------


class TestIdaesThermo:
    def test_returns_h_cp_k(self, monkeypatch):
        _install_fake_idaes(monkeypatch)
        p = IdaesProvider()
        p.initialize(IdaesProviderConfig(), FlowsheetConfig())
        props = p.get_thermo_properties({"T": 300, "P": 101325, "z": {"Water": 0.8, "Toluene": 0.2}})
        assert "H" in props and "Cp" in props and "K_values" in props
        assert props["K_values"]["Water"] == pytest.approx(1.1)


# ---------------------------------------------------------------------------
# 7. BaseUnitMixin delegation (params fix)
# ---------------------------------------------------------------------------


class TestBaseUnitMixinIdaesDelegation:
    def test_pump_delegates_with_correct_delta_p(self, monkeypatch):
        _install_fake_idaes(monkeypatch)
        from processforge.flowsheet import Flowsheet

        cfg = json.load(open("flowsheets/hydraulic-chain-idaes.json"))
        # Flowsheet does dict-based init; mock stays via sys.modules
        fs = Flowsheet(cfg)
        result = fs.run()
        # pump_2 has delta_p 200k — would be 100k if delegation used default
        assert result["after_pump_2"]["P"] == pytest.approx(293162.5)
        assert result["product"]["P"] == pytest.approx(292162.5)

    def test_units_store_params(self):
        from processforge.units.pump import Pump
        from processforge.units.valve import Valve
        from processforge.units.strainer import Strainer
        from processforge.units.pipes import Pipes

        assert Pump("x", delta_p=123).params["delta_p"] == 123
        assert Valve("x", pressure_ratio=0.3).params["pressure_ratio"] == 0.3
        assert Strainer("x", delta_p=777).params["delta_p"] == 777
        assert Pipes("x", delta_p=555, diameter=0.2).params["diameter"] == 0.2


# ---------------------------------------------------------------------------
# 8. Flowsheet integration — validation + end-to-end
# ---------------------------------------------------------------------------


class TestFlowsheetIntegration:
    def test_hydraulic_chain_idaes_validates(self):
        from processforge.utils.validate_flowsheet import validate_flowsheet_dict

        cfg = json.load(open("flowsheets/hydraulic-chain-idaes.json"))
        assert validate_flowsheet_dict(cfg, "hydraulic-chain-idaes") is not None

    def test_closed_loop_idaes_validates(self):
        from processforge.utils.validate_flowsheet import validate_flowsheet_dict

        cfg = json.load(open("flowsheets/closed-loop-chain-idaes.json"))
        assert validate_flowsheet_dict(cfg, "closed-loop-chain-idaes") is not None

    def test_closed_loop_dynamic_runs(self, monkeypatch):
        _install_fake_idaes(monkeypatch)
        from processforge.flowsheet import Flowsheet

        cfg = json.load(open("flowsheets/closed-loop-chain-idaes.json"))
        result = Flowsheet(cfg).run()
        assert "recycle_to_tank" in result
        # dynamic result has list for P
        assert isinstance(result["recycle_to_tank"]["P"], list)
        assert len(result["recycle_to_tank"]["P"]) == 21

    def test_idaes_flowsheet_fails_without_idaes(self, monkeypatch):
        # remove fake if present
        for k in list(sys.modules):
            if k.startswith("idaes"):
                monkeypatch.delitem(sys.modules, k, raising=False)
        from processforge.flowsheet import Flowsheet

        cfg = json.load(open("flowsheets/hydraulic-chain-idaes.json"))
        with pytest.raises(ProviderNotAvailableError):
            Flowsheet(cfg).run()
