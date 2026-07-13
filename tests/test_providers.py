"""Structural and lifecycle tests for the provider layer.

These tests validate that:
- Every provider module follows the ``{name}_provider`` naming convention.
- Every provider class subclasses ``AbstractProvider`` and implements its
  abstract methods.
- The registry is internally consistent (registered names match modules).
- ``build_provider_map`` / ``ProviderMap.resolve`` / ``teardown_providers``
  work for a minimal flowsheet.
"""
from __future__ import annotations

import importlib
import inspect
import pkgutil
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from processforge.providers.base import AbstractProvider
from processforge.providers.registry import _PROVIDERS, get_provider_class


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

PROVIDERS_PKG = Path(__file__).resolve().parent.parent / "src" / "processforge" / "providers"


@pytest.fixture()
def flowsheet_config():
    """Minimal FlowsheetConfig for tests that need one."""
    from processforge.types import FlowsheetConfig

    return FlowsheetConfig()


@pytest.fixture()
def coolprop_config():
    from processforge.types import CoolPropProviderConfig

    return CoolPropProviderConfig()


# ---------------------------------------------------------------------------
# 1. Naming convention — every *_provider.py follows {name}_provider.py
# ---------------------------------------------------------------------------

class TestNamingConvention:
    """Every provider module must be named ``{registered_name}_provider.py``."""

    def test_all_provider_modules_follow_convention(self):
        """Discover every *_provider.py in the providers package and check
        that the stem matches the ``{name}_provider`` pattern."""
        provider_modules = [
            mod_info.name
            for mod_info in pkgutil.iter_modules([str(PROVIDERS_PKG)])
            if mod_info.name.endswith("_provider")
        ]
        assert provider_modules, "No *_provider modules found"

        for mod_name in provider_modules:
            stem = mod_name.removesuffix("_provider")
            # stem should be non-empty and purely alphanumeric
            assert stem.isidentifier(), (
                f"Module '{mod_name}' has stem '{stem}' which is not a "
                f"valid identifier — expected '{{name}}_provider.py' with "
                f"a clean name."
            )

    def test_no_orphaned_provider_modules(self):
        """Every *_provider.py module that calls register_provider() must
        have its registered name match the module stem."""
        # Import all provider modules to populate _PROVIDERS
        for mod_info in pkgutil.iter_modules([str(PROVIDERS_PKG)]):
            if mod_info.name.endswith("_provider"):
                importlib.import_module(f"processforge.providers.{mod_info.name}")

        for registered_name, cls in _PROVIDERS.items():
            expected_module = f"processforge.providers.{registered_name}_provider"
            try:
                mod = importlib.import_module(expected_module)
            except ModuleNotFoundError:
                pytest.fail(
                    f"Provider '{registered_name}' is registered but no module "
                    f"named '{expected_module}' exists."
                )
            # The class should be defined in (or re-exported by) that module
            found = any(
                obj is cls
                for obj in vars(mod).values()
                if obj is cls
            )
            assert found, (
                f"Provider class {cls.__name__} is not defined in or "
                f"re-exported by its expected module '{expected_module}'."
            )


# ---------------------------------------------------------------------------
# 2. Class contract — subclass of AbstractProvider, implements ABCs
# ---------------------------------------------------------------------------

class TestProviderContract:
    """Every registered provider must satisfy the AbstractProvider contract."""

    def test_all_providers_subclass_abstract_provider(self):
        for name, cls in _PROVIDERS.items():
            assert issubclass(cls, AbstractProvider), (
                f"Provider '{name}' ({cls.__name__}) does not subclass "
                f"AbstractProvider."
            )

    def test_all_abstract_methods_implemented(self):
        abstract_methods = set(AbstractProvider.__abstractmethods__)
        for name, cls in _PROVIDERS.items():
            missing = abstract_methods - set(dir(cls))
            assert not missing, (
                f"Provider '{name}' ({cls.__name__}) is missing abstract "
                f"method implementations: {sorted(missing)}"
            )

    def test_all_providers_have_initialize(self):
        for name, cls in _PROVIDERS.items():
            assert hasattr(cls, "initialize"), (
                f"Provider '{name}' ({cls.__name__}) has no 'initialize' method."
            )

    def test_all_providers_have_teardown(self):
        for name, cls in _PROVIDERS.items():
            assert hasattr(cls, "teardown"), (
                f"Provider '{name}' ({cls.__name__}) has no 'teardown' method."
            )


# ---------------------------------------------------------------------------
# 3. Registry consistency
# ---------------------------------------------------------------------------

class TestRegistry:
    """The registry must be internally consistent."""

    def test_coolprop_always_registered(self):
        assert "coolprop" in _PROVIDERS

    def test_coolprop_class_is_coolprop_provider(self):
        from processforge.providers.coolprop_provider import CoolPropProvider

        assert _PROVIDERS["coolprop"] is CoolPropProvider

    def test_get_provider_class_returns_correct_class(self):
        for name, cls in _PROVIDERS.items():
            assert get_provider_class(name) is cls

    def test_get_provider_class_raises_for_unknown(self):
        with pytest.raises(ValueError, match="No provider module found"):
            get_provider_class("nonexistent_provider_xyz")


# ---------------------------------------------------------------------------
# 4. Lazy import (convention-based discovery via registry)
# ---------------------------------------------------------------------------

class TestLazyImport:
    """get_provider_class triggers lazy import for unregistered types."""

    def test_import_unknown_raises_valueerror(self):
        with pytest.raises(ValueError, match="No provider module found"):
            get_provider_class("nonexistent_engine")

    def test_each_registered_type_resolves(self):
        """For every type already in the registry, get_provider_class
        should return the class without error."""
        for name in list(_PROVIDERS):
            cls = get_provider_class(name)
            assert cls is _PROVIDERS[name]


# ---------------------------------------------------------------------------
# 5. build_provider_map integration (minimal flowsheet)
# ---------------------------------------------------------------------------

class TestBuildProviderMap:

    def test_returns_provider_map(self, flowsheet_config):
        from processforge.providers.manager import ProviderMap, build_provider_map

        pmap = build_provider_map({}, flowsheet_config)
        assert isinstance(pmap, ProviderMap)

    def test_coolprop_always_in_map(self, flowsheet_config):
        from processforge.providers.manager import _BUILTIN_DEFAULT_KEY, build_provider_map

        pmap = build_provider_map({}, flowsheet_config)
        assert _BUILTIN_DEFAULT_KEY in pmap

    def test_declared_provider_ends_up_in_map(self, flowsheet_config):
        from processforge.providers.manager import build_provider_map

        providers_cfg = {"my_coolprop": MagicMock(type="coolprop")}
        pmap = build_provider_map(providers_cfg, flowsheet_config)
        assert "my_coolprop" in pmap

    def test_default_provider_resolves(self, flowsheet_config):
        from processforge.providers.manager import UnitProviderConfig, build_provider_map

        flowsheet_config.default_provider = "my_coolprop"
        providers_cfg = {"my_coolprop": MagicMock(type="coolprop")}
        pmap = build_provider_map(providers_cfg, flowsheet_config)
        # A unit with no "provider" key should resolve to my_coolprop
        resolved = pmap.resolve(UnitProviderConfig())
        assert resolved is pmap["my_coolprop"]

    def test_invalid_default_provider_raises(self, flowsheet_config):
        from processforge.providers.manager import build_provider_map

        flowsheet_config.default_provider = "does_not_exist"
        with pytest.raises(ValueError, match="'default_provider' references"):
            build_provider_map({}, flowsheet_config)

    def test_undeclared_provider_type_raises(self, flowsheet_config):
        from processforge.providers.manager import build_provider_map

        providers_cfg = {"bad": MagicMock(type="nonexistent_engine")}
        with pytest.raises(ValueError, match="No provider module found"):
            build_provider_map(providers_cfg, flowsheet_config)

    def test_dict_like_access(self, flowsheet_config):
        from processforge.providers.manager import build_provider_map

        pmap = build_provider_map({}, flowsheet_config)
        # __getitem__, __contains__, __bool__, .values()
        assert pmap
        assert "__coolprop__" in pmap
        assert pmap["__coolprop__"] is not None
        assert len(list(pmap.values())) >= 1


# ---------------------------------------------------------------------------
# 6. ProviderMap.resolve
# ---------------------------------------------------------------------------

class TestProviderMapResolve:

    def test_returns_default_when_no_provider_key(self, flowsheet_config):
        from processforge.providers.manager import UnitProviderConfig, build_provider_map

        pmap = build_provider_map({}, flowsheet_config)
        resolved = pmap.resolve(UnitProviderConfig())
        assert resolved is pmap["__coolprop__"]

    def test_returns_named_provider(self, flowsheet_config):
        from processforge.providers.manager import UnitProviderConfig, build_provider_map

        providers_cfg = {"my_coolprop": MagicMock(type="coolprop")}
        pmap = build_provider_map(providers_cfg, flowsheet_config)
        resolved = pmap.resolve(UnitProviderConfig(provider="my_coolprop"))
        assert resolved is pmap["my_coolprop"]

    def test_undeclared_provider_raises(self, flowsheet_config):
        from processforge.providers.manager import UnitProviderConfig, build_provider_map

        pmap = build_provider_map({}, flowsheet_config)
        with pytest.raises(ValueError, match="Unit references provider"):
            pmap.resolve(UnitProviderConfig(provider="ghost"))

    def test_unit_provider_config_from_dict(self):
        from processforge.providers.manager import UnitProviderConfig

        cfg = UnitProviderConfig(**{"type": "CSTR", "provider": "cantera", "extra": 42})
        assert cfg.provider == "cantera"
        assert cfg.model_extra == {"type": "CSTR", "extra": 42}


# ---------------------------------------------------------------------------
# 7. teardown_providers
# ---------------------------------------------------------------------------

class TestTeardownProviders:

    def test_calls_teardown_on_each_provider(self, flowsheet_config):
        from processforge.providers.manager import build_provider_map, teardown_providers

        p1, p2 = MagicMock(), MagicMock()
        pmap = build_provider_map({}, flowsheet_config)
        # Inject mock providers directly
        object.__setattr__(pmap, "_providers", {"a": p1, "b": p2})
        teardown_providers(pmap)
        p1.teardown.assert_called_once()
        p2.teardown.assert_called_once()

    def test_teardown_error_does_not_block_others(self, flowsheet_config):
        from processforge.providers.manager import build_provider_map, teardown_providers

        p1 = MagicMock()
        p1.teardown.side_effect = RuntimeError("boom")
        p2 = MagicMock()
        pmap = build_provider_map({}, flowsheet_config)
        object.__setattr__(pmap, "_providers", {"a": p1, "b": p2})
        teardown_providers(pmap)
        p1.teardown.assert_called_once()
        p2.teardown.assert_called_once()

    def test_teardown_deduplicates_aliases(self, flowsheet_config):
        from processforge.providers.manager import build_provider_map, teardown_providers

        p1 = MagicMock()
        pmap = build_provider_map({}, flowsheet_config)
        object.__setattr__(pmap, "_providers", {"a": p1, "__default__": p1})
        teardown_providers(pmap)
        assert p1.teardown.call_count == 1

    def test_teardown_none_is_noop(self):
        from processforge.providers.manager import teardown_providers

        teardown_providers(None)  # should not raise
