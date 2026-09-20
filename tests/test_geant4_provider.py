"""Tests for the Geant4 provider (in-process architecture).

These tests validate:
- The provider registers correctly.
- The provider satisfies the AbstractProvider contract.
- Geant4ProviderConfig parses from dict correctly.
- The Geant4 schema validates a reactor_core geometry configuration.
- The strategy registry and dispatch work correctly.
- A Geant4 flowsheet (flowsheets/geant4/msre_shielding.json) parses correctly.
"""

from __future__ import annotations

import importlib
import importlib.util
import json
from pathlib import Path

import pytest

from processforge.providers.base import AbstractProvider
from processforge.providers.registry import _PROVIDERS
from processforge.types import (
    FlowsheetConfig,
    Geant4ProviderConfig,
    MaterialDef,
    UnitConfig,
)

# Ensure the geant4 provider module is imported (self-registers on import)
importlib.import_module("processforge.providers.geant4")


def test_geant4_provider_is_registered():
    assert "geant4" in _PROVIDERS
    assert issubclass(_PROVIDERS["geant4"], AbstractProvider)


def test_geant4_provider_config_from_dict():
    cfg = Geant4ProviderConfig.from_dict(
        {
            "type": "geant4",
            "url": "http://localhost:9003",
            "output_dir": "outputs/geant4",
        }
    )
    assert cfg.type == "geant4"
    assert cfg.url == "http://localhost:9003"
    assert cfg.output_dir == "outputs/geant4"


def test_flowsheet_config_parses_geant4_provider():
    cfg = Geant4ProviderConfig.from_dict({"type": "geant4"})
    assert cfg.output_dir == "outputs/geant4"


def test_geant4_schema_reactor_core_geometry():
    from processforge.schemas.geant4.geant4_model import (
        Geant4Setting,
        ReactorCoreGeometryConfig,
    )

    context = {"materials": {"salt", "graphite", "inconel", "helium", "inor"}}
    geo = ReactorCoreGeometryConfig.model_validate(
        {
            "type": "reactor_core",
            "core_radius": 72.5,
            "core_height": 160.0,
            "reflector_thickness": 50.0,
            "vessel_thickness": 5.0,
            "gap_thickness": 3.0,
            "structure_thickness": 10.0,
            "core_material": "salt",
            "reflector_material": "graphite",
            "vessel_material": "inconel",
            "gap_material": "helium",
            "structure_material": "inor",
            "source": {
                "particle": "neutron",
                "energy_MeV": 2.0,
                "distribution": "isotropic",
                "position": {"x": 0.0, "y": 0.0, "z": 0.0},
            },
        },
        context=context,
    )
    assert geo.core_radius == 72.5
    assert geo.source.energy_MeV == 2.0


def test_geant4_schema_rejects_unknown_material():
    from processforge.schemas.geant4.geant4_model import ReactorCoreGeometryConfig

    context = {"materials": {"salt"}}
    with pytest.raises(ValueError):
        ReactorCoreGeometryConfig.model_validate(
            {
                "type": "reactor_core",
                "core_radius": 72.5,
                "core_height": 160.0,
                "core_material": "salt",
                "reflector_material": "graphite",
            },
            context=context,
        )


def test_geant4_setting_defaults():
    from processforge.schemas.geant4.geant4_model import Geant4Setting

    cfg = Geant4Setting()
    assert cfg.physics_list == "FTFP_BERT"
    assert cfg.events == 1000


def test_geant4_strategy_registry_contains_shielding_attenuation():
    from processforge.providers.geant4.strategies import get_registered_sim_types

    assert "shielding_attenuation" in get_registered_sim_types()


def test_msre_shielding_flowsheet_parses():
    path = Path(__file__).parents[1] / "flowsheets" / "geant4" / "msre_shielding.json"
    raw = json.loads(path.read_text())
    cfg = FlowsheetConfig.from_dict(raw)

    assert cfg.providers["geant4"].type == "geant4"
    assert cfg.units["geant4_solver"].sim_type == "shielding_attenuation"
    assert cfg.units["geant4_solver"].solver_config["events"] == 10000
    assert cfg.materials["salt"].density == 2.2


@pytest.mark.skipif(
    importlib.util.find_spec("geant4") is None,
    reason="geant4 (geant4-pybind) is not installed",
)
def test_geant4_provider_initializes_with_valid_config():
    from processforge.providers.geant4.provider import Geant4Provider

    provider = Geant4Provider()
    cfg = Geant4ProviderConfig.from_dict(
        {"type": "geant4", "output_dir": "outputs/geant4"}
    )
    flowsheet_cfg = FlowsheetConfig.from_dict(
        {
            "providers": {"geant4": cfg.model_dump()},
            "materials": {
                "water": {
                    "id": 1,
                    "density": 1.0,
                    "density_units": "g/cm3",
                    "temperature": 300.0,
                    "elements": [
                        {"element": "H", "percent": 0.111, "percent_type": "wo"}
                    ],
                    "extra": {
                        "elements": [
                            {"element": "O", "percent": 0.889, "percent_type": "wo"}
                        ]
                    },
                }
            },
        }
    )

    provider.initialize(cfg, flowsheet_cfg)
    assert provider.is_initialized


def test_geant4_provider_validate_material_flags_missing_density():
    from processforge.providers.geant4.provider import Geant4Provider

    mat = MaterialDef.from_dict(
        {
            "id": 1,
            "density_units": "g/cm3",
        }
    )
    errors = Geant4Provider.validate_material("mat1", mat, None)
    assert any("missing 'density'" in e for e in errors)


def test_geant4_provider_validate_material_requires_density_units():
    from processforge.providers.geant4.provider import Geant4Provider

    mat = MaterialDef.from_dict(
        {
            "id": 1,
            "density": 1.0,
        }
    )
    errors = Geant4Provider.validate_material("mat1", mat, None)
    assert any("missing 'density_units'" in e for e in errors)
