"""Tests for pf init and pf validate commands."""
from __future__ import annotations

import json
import os
from pathlib import Path
from unittest.mock import patch

import pytest

from processforge.lock import read_lock, write_lock, LOCK_VERSION
from processforge.compose import generate_compose, COMPOSE_FILENAME
from processforge.providers.registry import (
    _PROVIDER_CATALOG,
    get_provider_docker_image,
    get_provider_default_port,
    is_containerized,
    list_providers,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

TESTS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = TESTS_DIR.parent


@pytest.fixture
def tmp_pf_dir(tmp_path):
    """Create a temporary .processforge/ directory."""
    pf_dir = tmp_path / ".processforge"
    pf_dir.mkdir()
    return pf_dir


@pytest.fixture
def coolprop_flowsheet(tmp_path):
    """Write a minimal flowsheet with only coolprop."""
    flowsheet = {
        "providers": {"coolprop": {"type": "coolprop"}},
        "streams": {},
        "units": {},
        "simulation": {"mode": "steady"},
    }
    path = tmp_path / "coolprop_only.json"
    with open(path, "w") as f:
        json.dump(flowsheet, f)
    return path


@pytest.fixture
def openmc_flowsheet(tmp_path):
    """Write a minimal flowsheet with openmc (containerized)."""
    flowsheet = {
        "providers": {
            "openmc": {
                "type": "openmc",
                "url": "http://localhost:9001",
                "output_dir": "openmc/run",
            }
        },
        "streams": {},
        "units": {},
        "simulation": {"mode": "steady"},
    }
    path = tmp_path / "openmc_flowsheet.json"
    with open(path, "w") as f:
        json.dump(flowsheet, f)
    return path


@pytest.fixture
def mixed_flowsheet(tmp_path):
    """Write a flowsheet with both pip and containerized providers."""
    flowsheet = {
        "providers": {
            "coolprop": {"type": "coolprop"},
            "openmc": {
                "type": "openmc",
                "url": "http://localhost:9001",
            },
        },
        "streams": {},
        "units": {},
        "simulation": {"mode": "steady"},
    }
    path = tmp_path / "mixed.json"
    with open(path, "w") as f:
        json.dump(flowsheet, f)
    return path


# ---------------------------------------------------------------------------
# Registry tests
# ---------------------------------------------------------------------------


class TestProviderCatalog:
    def test_coolprop_no_docker(self):
        assert _PROVIDER_CATALOG["coolprop"]["docker_image"] is None
        assert _PROVIDER_CATALOG["coolprop"]["default_port"] is None

    def test_openmc_has_docker(self):
        assert _PROVIDER_CATALOG["openmc"]["docker_image"] is not None
        assert _PROVIDER_CATALOG["openmc"]["default_port"] == 9001

    def test_festim_in_catalog(self):
        assert "festim" in _PROVIDER_CATALOG
        assert _PROVIDER_CATALOG["festim"]["docker_image"] is not None
        assert _PROVIDER_CATALOG["festim"]["default_port"] == 9002

    def test_get_provider_docker_image(self):
        assert get_provider_docker_image("openmc") is not None
        assert get_provider_docker_image("coolprop") is None

    def test_get_provider_docker_image_unknown(self):
        with pytest.raises(ValueError, match="Unknown provider type"):
            get_provider_docker_image("nonexistent")

    def test_get_provider_default_port(self):
        assert get_provider_default_port("openmc") == 9001
        assert get_provider_default_port("coolprop") is None

    def test_is_containerized(self):
        assert is_containerized("openmc") is True
        assert is_containerized("festim") is True
        assert is_containerized("coolprop") is False
        assert is_containerized("cantera") is False

    def test_list_providers_includes_docker_fields(self):
        providers = list_providers()
        openmc = next(p for p in providers if p["type"] == "openmc")
        assert "docker_image" in openmc
        assert "default_port" in openmc
        assert openmc["docker_image"] is not None
        assert openmc["default_port"] == 9001


# ---------------------------------------------------------------------------
# Lock file tests
# ---------------------------------------------------------------------------


class TestLockFile:
    def test_read_lock_missing(self, tmp_pf_dir):
        assert read_lock(str(tmp_pf_dir)) is None

    def test_write_and_read_lock(self, tmp_pf_dir):
        providers = {
            "openmc": {"docker_image": "ghcr.io/test:latest", "url": "http://localhost:9001"},
            "coolprop": {"docker_image": None, "url": None},
        }
        write_lock(str(tmp_pf_dir), "test.json", providers, "0.2.34")
        lock = read_lock(str(tmp_pf_dir))

        assert lock is not None
        assert lock["version"] == LOCK_VERSION
        assert lock["processforge_version"] == "0.2.34"
        assert lock["flowsheet"] == "test.json"
        assert "openmc" in lock["providers"]
        assert "coolprop" in lock["providers"]
        assert lock["providers"]["openmc"]["url"] == "http://localhost:9001"

    def test_lock_creates_directory(self, tmp_path):
        pf_dir = tmp_path / ".processforge"
        providers = {"coolprop": {"docker_image": None, "url": None}}
        write_lock(str(pf_dir), "test.json", providers, "0.2.34")
        assert read_lock(str(pf_dir)) is not None


# ---------------------------------------------------------------------------
# Compose generation tests
# ---------------------------------------------------------------------------


class TestComposeGeneration:
    def test_generate_compose_with_docker_providers(self, tmp_pf_dir):
        docker_providers = {
            "openmc": {
                "docker_image": "ghcr.io/urjanova/processforge-openmc:latest",
                "port": 9001,
            },
        }
        generate_compose(str(tmp_pf_dir), docker_providers)
        compose_path = tmp_pf_dir / COMPOSE_FILENAME
        assert compose_path.exists()

        content = compose_path.read_text()
        assert "services:" in content
        assert "openmc:" in content
        assert "ghcr.io/urjanova/processforge-openmc:latest" in content
        assert "9001:9001" in content

    def test_generate_compose_empty_when_no_docker(self, tmp_pf_dir):
        generate_compose(str(tmp_pf_dir), {})
        compose_path = tmp_pf_dir / COMPOSE_FILENAME
        assert not compose_path.exists()

    def test_generate_compose_multiple_providers(self, tmp_pf_dir):
        docker_providers = {
            "openmc": {"docker_image": "ghcr.io/test/openmc:latest", "port": 9001},
            "festim": {"docker_image": "ghcr.io/test/festim:latest", "port": 9002},
        }
        generate_compose(str(tmp_pf_dir), docker_providers)
        content = (tmp_pf_dir / COMPOSE_FILENAME).read_text()
        assert "openmc:" in content
        assert "festim:" in content

    def test_generate_compose_custom_docker_image(self, tmp_pf_dir):
        docker_providers = {
            "openmc": {
                "docker_image": "my-org/custom-openmc:v2",
                "port": 9001,
            },
        }
        generate_compose(str(tmp_pf_dir), docker_providers)
        content = (tmp_pf_dir / COMPOSE_FILENAME).read_text()
        assert "my-org/custom-openmc:v2" in content
        assert "ghcr.io/urjanova/processforge-openmc" not in content


# ---------------------------------------------------------------------------
# Extract providers tests
# ---------------------------------------------------------------------------


class TestExtractProviders:
    def test_extract_from_coolprop_flowsheet(self, coolprop_flowsheet):
        from processforge.cli.common import extract_providers

        providers = extract_providers(str(coolprop_flowsheet))
        assert "coolprop" in providers
        assert providers["coolprop"]["type"] == "coolprop"

    def test_extract_from_openmc_flowsheet(self, openmc_flowsheet):
        from processforge.cli.common import extract_providers

        providers = extract_providers(str(openmc_flowsheet))
        assert "openmc" in providers
        assert providers["openmc"]["type"] == "openmc"
        assert providers["openmc"]["url"] == "http://localhost:9001"

    def test_extract_from_mixed_flowsheet(self, mixed_flowsheet):
        from processforge.cli.common import extract_providers

        providers = extract_providers(str(mixed_flowsheet))
        assert "coolprop" in providers
        assert "openmc" in providers

    def test_extract_missing_file(self, tmp_path):
        from processforge.cli.common import extract_providers

        with pytest.raises(SystemExit):
            extract_providers(str(tmp_path / "nonexistent.json"))
