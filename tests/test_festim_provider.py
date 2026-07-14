"""Tests for the FESTIM provider.

These tests validate:
- The provider registers correctly.
- The provider satisfies the AbstractProvider contract.
- Material validation catches missing D_0/E_D.
- Sim type registry works correctly.
- Build helpers produce valid FESTIM objects (when FESTIM is installed).
"""
from __future__ import annotations

import importlib

import pytest

from processforge.providers.base import AbstractProvider
from processforge.providers.registry import _PROVIDERS

# Ensure the festim provider module is imported (self-registers on import)
importlib.import_module("processforge.providers.festim_provider")


# ---------------------------------------------------------------------------
# 1. Registration
# ---------------------------------------------------------------------------


class TestFestimRegistration:
    """FESTIM provider must be registered in the global provider registry."""

    def test_festim_registered(self):
        assert "festim" in _PROVIDERS

    def test_festim_class_is_festim_provider(self):
        from processforge.providers.festim_provider import FestimProvider

        assert _PROVIDERS["festim"] is FestimProvider

    def test_festim_subclasses_abstract_provider(self):
        from processforge.providers.festim_provider import FestimProvider

        assert issubclass(FestimProvider, AbstractProvider)

    def test_all_abstract_methods_implemented(self):
        from processforge.providers.festim_provider import FestimProvider

        abstract_methods = set(AbstractProvider.__abstractmethods__)
        missing = abstract_methods - set(dir(FestimProvider))
        assert not missing, (
            f"FestimProvider is missing abstract method implementations: {sorted(missing)}"
        )


# ---------------------------------------------------------------------------
# 2. Sim type registry
# ---------------------------------------------------------------------------


class TestFestimSimTypeRegistry:
    """The FESTIM sim type registry must work correctly."""

    def test_hydrogen_transport_registered(self):
        from processforge.providers.festim_provider import get_registered_sim_types

        registry = get_registered_sim_types()
        assert "hydrogen_transport" in registry

    def test_registry_returns_copy(self):
        from processforge.providers.festim_provider import get_registered_sim_types

        r1 = get_registered_sim_types()
        r2 = get_registered_sim_types()
        assert r1 is not r2  # must be a copy
        assert r1 == r2

    def test_custom_sim_type_registration(self):
        from processforge.providers.festim_provider import (
            FestimSimStrategy,
            get_registered_sim_types,
            register_festim_sim_type,
        )

        class _DummyStrategy(FestimSimStrategy):
            def build(self, F, unit_cfg, materials_map, helpers):
                pass

        register_festim_sim_type("test_dummy", _DummyStrategy)
        try:
            registry = get_registered_sim_types()
            assert "test_dummy" in registry
        finally:
            # Cleanup — remove the test entry
            from processforge.providers.festim_provider import _SIM_TYPE_REGISTRY

            _SIM_TYPE_REGISTRY.pop("test_dummy", None)


# ---------------------------------------------------------------------------
# 3. Material validation
# ---------------------------------------------------------------------------


class TestFestimMaterialValidation:
    """validate_material should enforce D_0 and E_D in extra."""

    def test_valid_material(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import MaterialDef, UnitConfig

        mat_def = MaterialDef(id=1, extra={"D_0": 1e-7, "E_D": 0.2})
        unit_cfg = UnitConfig(type="SolverUnit", material=0, sim_type="hydrogen_transport")
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert errors == []

    def test_missing_D_0(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import MaterialDef, UnitConfig

        mat_def = MaterialDef(id=1, extra={"E_D": 0.2})
        unit_cfg = UnitConfig(type="SolverUnit", material=0, sim_type="hydrogen_transport")
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 1
        assert "D_0" in errors[0]

    def test_missing_E_D(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import MaterialDef, UnitConfig

        mat_def = MaterialDef(id=1, extra={"D_0": 1e-7})
        unit_cfg = UnitConfig(type="SolverUnit", material=0, sim_type="hydrogen_transport")
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 1
        assert "E_D" in errors[0]

    def test_missing_both(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import MaterialDef, UnitConfig

        mat_def = MaterialDef(id=1, extra={})
        unit_cfg = UnitConfig(type="SolverUnit", material=0, sim_type="hydrogen_transport")
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 2

    def test_missing_extra_entirely(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import MaterialDef, UnitConfig

        mat_def = MaterialDef(id=1)
        unit_cfg = UnitConfig(type="SolverUnit", material=0, sim_type="hydrogen_transport")
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 2


# ---------------------------------------------------------------------------
# 4. ProviderConfig
# ---------------------------------------------------------------------------


class TestFestimProviderConfig:
    """FestimProviderConfig should parse from dict correctly."""

    def test_default_config(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig()
        assert cfg.type == "festim"
        assert cfg.output_dir == "outputs/festim"

    def test_from_dict_defaults(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict({"type": "festim"})
        assert cfg.output_dir == "outputs/festim"

    def test_from_dict_custom_output_dir(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict({"type": "festim", "output_dir": "/tmp/festim"})
        assert cfg.output_dir == "/tmp/festim"

    def test_in_registry(self):
        from processforge.types import _PROVIDER_CONFIG_REGISTRY

        assert "festim" in _PROVIDER_CONFIG_REGISTRY

    def test_in_union(self):
        from processforge.types import provider_config_from_dict

        cfg = provider_config_from_dict({"type": "festim"})
        from processforge.types import FestimProviderConfig

        assert isinstance(cfg, FestimProviderConfig)
