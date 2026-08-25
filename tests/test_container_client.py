"""Tests for the CLI-side container provider HTTP client serialization."""

from processforge.providers.container_client import ContainerProviderClient
from processforge.types import UnitConfig


def test_serialize_unit_config_includes_geometry_config():
    """geometry_config must survive the CLI -> container body serialization.

    Regression test for the bug where geometry_config (a known UnitConfig
    field) was dropped by the hardcoded serialization list, so the container
    received an empty geometry_config and failed validation.
    """
    geometry_config = {
        "type": "reactor_core",
        "core_radius": 72.5,
        "core_height": 160.0,
        "core_material": "salt",
    }
    uc = UnitConfig(
        type="SolverUnit",
        provider="openmc",
        material=3,
        sim_type="eigenvalue_reactor",
        solver_config={"batches": 20},
        geometry_config=geometry_config,
    )

    body = ContainerProviderClient._serialize_unit_config(uc)

    assert body["geometry_config"] == geometry_config
    assert body["solver_config"] == {"batches": 20}
    assert body["sim_type"] == "eigenvalue_reactor"
