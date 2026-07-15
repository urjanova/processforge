"""Tests for the FESTIM provider (HTTP client architecture).

These tests validate:
- The provider registers correctly.
- The provider satisfies the AbstractProvider contract.
- Material validation catches missing D_0/E_D.
- FestimProviderConfig parses url from dict correctly.
- The HTTP client serialises requests and deserialises responses.
"""
from __future__ import annotations

import importlib
import json
from unittest.mock import MagicMock, patch

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
# 2. Material validation
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
# 3. ProviderConfig
# ---------------------------------------------------------------------------


class TestFestimProviderConfig:
    """FestimProviderConfig should parse from dict correctly."""

    def test_default_config(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig()
        assert cfg.type == "festim"
        assert cfg.url is None
        assert cfg.output_dir == "outputs/festim"

    def test_from_dict_defaults(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict({"type": "festim"})
        assert cfg.url is None
        assert cfg.output_dir == "outputs/festim"

    def test_from_dict_with_url(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict({
            "type": "festim",
            "url": "http://localhost:9002",
        })
        assert cfg.url == "http://localhost:9002"

    def test_from_dict_cloud_url(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict({
            "type": "festim",
            "url": "https://festim-production.up.railway.app",
        })
        assert cfg.url == "https://festim-production.up.railway.app"

    def test_from_dict_custom_output_dir(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict({
            "type": "festim",
            "url": "http://localhost:9002",
            "output_dir": "/tmp/festim",
        })
        assert cfg.output_dir == "/tmp/festim"

    def test_from_dict_with_docker_image(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict({
            "type": "festim",
            "docker_image": "my-org/custom-festim:v2",
        })
        assert cfg.docker_image == "my-org/custom-festim:v2"

    def test_docker_image_defaults_to_none(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig()
        assert cfg.docker_image is None

    def test_in_registry(self):
        from processforge.types import _PROVIDER_CONFIG_REGISTRY

        assert "festim" in _PROVIDER_CONFIG_REGISTRY

    def test_in_union(self):
        from processforge.types import provider_config_from_dict

        cfg = provider_config_from_dict({"type": "festim"})
        from processforge.types import FestimProviderConfig

        assert isinstance(cfg, FestimProviderConfig)


# ---------------------------------------------------------------------------
# 4. Catalog metadata
# ---------------------------------------------------------------------------


class TestFestimCatalog:
    """Festim catalog entry must reflect Docker service architecture."""

    def test_has_docker_image(self):
        from processforge.providers.registry import _PROVIDER_CATALOG

        info = _PROVIDER_CATALOG["festim"]
        assert info["docker_image"] == "ghcr.io/urjanova/processforge-festim:latest"

    def test_has_default_port(self):
        from processforge.providers.registry import _PROVIDER_CATALOG

        info = _PROVIDER_CATALOG["festim"]
        assert info["default_port"] == 9002

    def test_is_containerized(self):
        from processforge.providers.registry import is_containerized

        assert is_containerized("festim") is True

    def test_no_optional_dep(self):
        from processforge.providers.registry import _PROVIDER_CATALOG

        info = _PROVIDER_CATALOG["festim"]
        assert info["optional_dep"] is None


# ---------------------------------------------------------------------------
# 5. Sim-type registry
# ---------------------------------------------------------------------------


class TestFestimSimTypeRegistry:
    """FESTIM sim_type registry must be importable and pre-seeded."""

    def test_get_registered_sim_types_callable(self):
        from processforge.providers.festim_provider import get_registered_sim_types

        result = get_registered_sim_types()
        assert isinstance(result, dict)

    def test_hydrogen_transport_registered(self):
        from processforge.providers.festim_provider import get_registered_sim_types

        registry = get_registered_sim_types()
        assert "hydrogen_transport" in registry
        assert "description" in registry["hydrogen_transport"]

    def test_register_new_sim_type(self):
        from processforge.providers.festim_provider import (
            _SIM_TYPE_REGISTRY,
            get_registered_sim_types,
            register_festim_sim_type,
        )

        register_festim_sim_type("custom_thing", description="a test type")
        try:
            registry = get_registered_sim_types()
            assert "custom_thing" in registry
        finally:
            del _SIM_TYPE_REGISTRY["custom_thing"]

    def test_returns_copy(self):
        from processforge.providers.festim_provider import get_registered_sim_types

        r1 = get_registered_sim_types()
        r2 = get_registered_sim_types()
        assert r1 == r2
        assert r1 is not r2


# ---------------------------------------------------------------------------
# 6. HTTP client behaviour
# ---------------------------------------------------------------------------


class TestFestimHTTPClient:
    """FestimProvider should communicate with the Docker service via HTTP."""

    def _make_provider(self):
        from processforge.providers.festim_provider import FestimProvider

        return FestimProvider()

    def _make_config(self, url="http://localhost:9002"):
        from processforge.types import FestimProviderConfig

        return FestimProviderConfig(url=url)

    def _make_flowsheet_config(self):
        from processforge.types import FlowsheetConfig, MaterialDef

        return FlowsheetConfig(
            materials={
                "steel": MaterialDef(id=1, extra={"D_0": 1e-7, "E_D": 0.2}),
            }
        )

    @patch("processforge.providers.festim_provider.urllib.request.urlopen")
    def test_initialize_health_check(self, mock_urlopen):
        health_resp = MagicMock()
        health_resp.read.return_value = json.dumps({
            "status": "ready",
            "provider_type": "festim",
        }).encode()
        health_resp.__enter__ = MagicMock(return_value=health_resp)
        health_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = health_resp

        provider = self._make_provider()
        provider.initialize(self._make_config(), self._make_flowsheet_config())

        assert provider._initialized is True
        assert provider._url == "http://localhost:9002"
        mock_urlopen.assert_called_once()

    @patch("processforge.providers.festim_provider.urllib.request.urlopen")
    def test_initialize_failure_raises(self, mock_urlopen):
        import urllib.error

        mock_urlopen.side_effect = urllib.error.URLError("connection refused")

        provider = self._make_provider()
        with pytest.raises(RuntimeError, match="unreachable"):
            provider.initialize(self._make_config(), self._make_flowsheet_config())

    @patch("processforge.providers.festim_provider.urllib.request.urlopen")
    def test_run_simulation_posts_json(self, mock_urlopen):
        from processforge.types import UnitConfig

        # Mock health check
        health_resp = MagicMock()
        health_resp.read.return_value = json.dumps({"status": "ready"}).encode()
        health_resp.__enter__ = MagicMock(return_value=health_resp)
        health_resp.__exit__ = MagicMock(return_value=False)

        # Mock run response
        run_resp = MagicMock()
        run_resp.read.return_value = json.dumps({
            "status": "completed",
            "sim_type": "hydrogen_transport",
            "scalars": {"hydrogen_flux": 42.0},
            "metadata": {"run_dir": "/data/festim"},
        }).encode()
        run_resp.__enter__ = MagicMock(return_value=run_resp)
        run_resp.__exit__ = MagicMock(return_value=False)

        mock_urlopen.side_effect = [health_resp, run_resp]

        provider = self._make_provider()
        provider.initialize(self._make_config(), self._make_flowsheet_config())

        unit_cfg = UnitConfig(
            type="SolverUnit",
            sim_type="hydrogen_transport",
            material=1,
            solver_config={"atol": 1e-6, "rtol": 1e-4},
        )
        result = provider.run_simulation(unit_cfg, {})

        assert result.status == "completed"
        assert result.sim_type == "hydrogen_transport"
        assert result.scalars["hydrogen_flux"] == 42.0

        # Verify the POST request was made with correct JSON body
        call_args = mock_urlopen.call_args_list
        post_call = call_args[1]  # second call is the POST /run
        req = post_call[0][0]
        assert req.full_url == "http://localhost:9002/run"
        assert req.method == "POST"

        body = json.loads(req.data.decode())
        assert body["unit_config"]["sim_type"] == "hydrogen_transport"
        assert "steel" in body["materials"]
        assert body["materials"]["steel"]["id"] == 1

    def test_run_simulation_before_init_raises(self):
        provider = self._make_provider()
        from processforge.types import UnitConfig

        unit_cfg = UnitConfig(type="SolverUnit", sim_type="hydrogen_transport")
        with pytest.raises(RuntimeError, match="not been initialized"):
            provider.run_simulation(unit_cfg, {})

    def test_teardown_resets_state(self):
        from processforge.providers.festim_provider import FestimProvider

        provider = FestimProvider()
        provider._initialized = True
        provider.teardown()
        assert provider._initialized is False

    @patch("processforge.providers.festim_provider.urllib.request.urlopen")
    def test_default_url_derived_from_port(self, mock_urlopen):
        from processforge.types import FestimProviderConfig

        health_resp = MagicMock()
        health_resp.read.return_value = json.dumps({"status": "ready"}).encode()
        health_resp.__enter__ = MagicMock(return_value=health_resp)
        health_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = health_resp

        provider = self._make_provider()
        # No url in config — should default to http://localhost:9002
        provider.initialize(FestimProviderConfig(), self._make_flowsheet_config())
        assert provider._url == "http://localhost:9002"

    def test_compute_unit_returns_none(self):
        provider = self._make_provider()
        assert provider.compute_unit("SolverUnit", {}, {}) is None

    def test_get_thermo_properties_raises(self):
        provider = self._make_provider()
        with pytest.raises(NotImplementedError):
            provider.get_thermo_properties({})
