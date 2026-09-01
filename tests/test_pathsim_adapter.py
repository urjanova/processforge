"""Tests for the PathSim adapter."""
from __future__ import annotations

import json
import os
import tempfile

import pytest

from processforge.fmu import build_fmu
from processforge.pathsim import ProcessForgeFMU

# Re-use the fake flowsheet fixture from the FMU tests.
from tests.test_fmu_solverunit import fake_flowsheet_path  # noqa: F401


@pytest.fixture
def fmu_and_manifest(fake_flowsheet_path: str) -> tuple[str, str]:
    pytest.importorskip("pythonfmu")
    pytest.importorskip("fmpy")
    pytest.importorskip("pathsim")

    with tempfile.TemporaryDirectory() as output_dir:
        fmu_path = build_fmu(
            fake_flowsheet_path, output_dir=output_dir, backend="scipy"
        )
        manifest_path = fmu_path.replace(".fmu", "_fmu_interface.json")
        yield fmu_path, manifest_path


def test_adapter_port_mapping(fmu_and_manifest: tuple[str, str]) -> None:
    """The adapter exposes logical port names backed by the manifest."""
    fmu_path, manifest_path = fmu_and_manifest

    plant = ProcessForgeFMU(fmu_path, manifest_path=manifest_path, dt=1.0)

    assert "fake_solver.solver_config.temperature" in plant.inputs
    assert "fake_solver.k_eff" in plant.outputs
    assert "fake_solver.power" in plant.outputs

    # Ports should be PathSim register entries (backend-specific objects).
    assert plant["fake_solver.solver_config.temperature"] is not None
    assert plant["fake_solver.k_eff"] is not None


def test_adapter_run_in_pathsim(fmu_and_manifest: tuple[str, str]) -> None:
    """Run a minimal PathSim simulation using the adapter."""
    pytest.importorskip("pathsim")
    fmu_path, manifest_path = fmu_and_manifest

    import pathsim as ps

    plant = ProcessForgeFMU(fmu_path, manifest_path=manifest_path, dt=1.0)

    # Drive the SolverUnit input with a constant block.
    setpoint = ps.blocks.Constant(value=950.0)
    sim = ps.Simulation(
        blocks=[setpoint, plant.block],
        connections=[
            ps.Connection(setpoint, plant["fake_solver.solver_config.temperature"]),
        ],
    )
    sim.run(5.0)

    # The output should reflect the fake provider result.
    assert abs(plant.value("fake_solver.k_eff") - 1.05) < 1e-6
