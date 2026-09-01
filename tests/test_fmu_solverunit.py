"""Tests for SolverUnit-aware FMU export."""
from __future__ import annotations

import json
import os
import tempfile

import pytest

from processforge.fmu import build_fmu
from processforge.providers.base import AbstractProvider
from processforge.providers.registry import register_provider
from processforge.quantity import Quantity
from processforge.types import EngineOutput, OutputField, OutputProvenance


class FakeProvider(AbstractProvider):
    """Minimal provider that returns a fixed EngineOutput for testing."""

    def __init__(self) -> None:
        self._initialized = False

    def initialize(self, provider_config, flowsheet_config) -> None:
        self._initialized = True

    def get_thermo_properties(self, stream: dict) -> dict:
        raise NotImplementedError

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        return None

    def teardown(self) -> None:
        self._initialized = False

    def run_simulation(self, unit_config, inlet: dict) -> EngineOutput:
        return EngineOutput(
            engine="fake",
            sim_type=unit_config.sim_type,
            status="completed",
            fields=[
                OutputField(
                    name="k_eff",
                    quantity=Quantity(value=1.05, unit=""),
                    kind="scalar",
                ),
                OutputField(
                    name="power",
                    quantity=Quantity(value=1_000_000.0, unit="W"),
                    kind="scalar",
                ),
            ],
            provenance=OutputProvenance(),
        )


# Register the fake provider under the "modelica" type so the JSON schema is
# satisfied, without pulling in any real external solver dependencies.
register_provider("modelica", FakeProvider)


_FAKE_FLOW_SHEET = {
    "metadata": {"name": "Fake SolverUnit Test", "version": "1.0"},
    "providers": {"modelica_provider": {"type": "modelica"}},
    "materials": {
        "salt": {
            "id": 1,
            "density": 2.2,
            "density_units": "g/cm3",
            "temperature": 900.0,
        }
    },
    "streams": {
        "feed": {"T": 900.0, "P": 101325.0, "flowrate": 1.0, "z": {"salt": 1.0}}
    },
    "units": {
        "fake_solver": {
            "type": "SolverUnit",
            "provider": "modelica_provider",
            "sim_type": "test_sim",
            "material": 1,
            "solver_config": {"temperature": 900.0},
            "fmu": {
                "inputs": ["solver_config.temperature"],
                "outputs": ["k_eff", "power"],
            },
        }
    },
    "simulation": {"mode": "steady"},
}


_FAKE_DYNAMIC_FLOW_SHEET = {
    "metadata": {"name": "Tank SolverUnit Test", "version": "1.0"},
    "providers": {"modelica_provider": {"type": "modelica"}},
    "materials": {
        "salt": {
            "id": 1,
            "density": 2.2,
            "density_units": "g/cm3",
            "temperature": 900.0,
        }
    },
    "streams": {
        "feed": {"T": 900.0, "P": 101325.0, "flowrate": 1.0, "z": {"salt": 1.0}}
    },
    "units": {
        "tank": {
            "type": "Tank",
            "in": "feed",
            "out": "tank_out",
            "material": 1,
            "volume": 10.0,
            "initial_level": 0.5,
            "outlet_flow": 0.5,
            "initial_n": {"salt": 10.0},
            "initial_T": 900.0,
            "P": 101325.0,
        },
        "fake_solver": {
            "type": "SolverUnit",
            "provider": "modelica_provider",
            "sim_type": "test_sim",
            "material": 1,
            "solver_config": {"temperature": 900.0},
            "fmu": {
                "inputs": ["solver_config.temperature"],
                "outputs": ["k_eff", "power"],
            },
        },
    },
    "simulation": {"mode": "dynamic", "t0": 0, "tf": 10, "dt": 1},
}


@pytest.fixture
def fake_flowsheet_path() -> str:
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "fake_solverunit.json")
        with open(path, "w", encoding="utf-8") as f:
            json.dump(_FAKE_FLOW_SHEET, f)
        yield path


@pytest.fixture
def fake_dynamic_flowsheet_path() -> str:
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "fake_dynamic_solverunit.json")
        with open(path, "w", encoding="utf-8") as f:
            json.dump(_FAKE_DYNAMIC_FLOW_SHEET, f)
        yield path


def test_build_solverunit_fmu(fake_flowsheet_path: str) -> None:
    """Export a SolverUnit flowsheet and verify the FMU + manifest are produced."""
    with tempfile.TemporaryDirectory() as output_dir:
        fmu_path = build_fmu(
            fake_flowsheet_path, output_dir=output_dir, backend="scipy"
        )
        assert os.path.isfile(fmu_path)

        manifest_path = fmu_path.replace(".fmu", "_fmu_interface.json")
        assert os.path.isfile(manifest_path)

        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        assert manifest["version"] == "1.0"
        assert manifest["flowsheet_name"] == "Fake SolverUnit Test"

        variables = {v["logical_name"]: v for v in manifest["variables"]}
        assert "fake_solver.solver_config.temperature" in variables
        assert variables["fake_solver.solver_config.temperature"]["causality"] == "input"
        assert "fake_solver.k_eff" in variables
        assert variables["fake_solver.k_eff"]["causality"] == "output"
        assert "fake_solver.power" in variables


def test_fmu_manifest_matches_slave_variables(fake_flowsheet_path: str) -> None:
    """The manifest attr_names must match the variables declared by pythonfmu."""
    pytest.importorskip("pythonfmu")
    pytest.importorskip("fmpy")

    with tempfile.TemporaryDirectory() as output_dir:
        fmu_path = build_fmu(
            fake_flowsheet_path, output_dir=output_dir, backend="scipy"
        )
        manifest_path = fmu_path.replace(".fmu", "_fmu_interface.json")
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        from fmpy import read_model_description

        md = read_model_description(fmu_path)
        fmu_vars = {sv.name for sv in md.modelVariables}
        for var in manifest["variables"]:
            assert var["attr_name"] in fmu_vars, f"missing {var['attr_name']}"


def test_solverunit_fmu_simulation(fake_flowsheet_path: str) -> None:
    """Run the exported FMU and check that SolverUnit outputs are populated."""
    pytest.importorskip("pythonfmu")
    pytest.importorskip("fmpy")

    with tempfile.TemporaryDirectory() as output_dir:
        fmu_path = build_fmu(
            fake_flowsheet_path, output_dir=output_dir, backend="scipy"
        )

        from fmpy import simulate_fmu

        result = simulate_fmu(
            fmu_path,
            stop_time=1.0,
            output_interval=1.0,
            output=["out_fake_solver_k_eff", "out_fake_solver_power"],
        )

        # result is a structured array; last row is final time.
        final = result[-1]
        assert abs(final["out_fake_solver_k_eff"] - 1.05) < 1e-9
        assert abs(final["out_fake_solver_power"] - 1_000_000.0) < 1.0


def test_build_solverunit_dynamic_fmu(fake_dynamic_flowsheet_path: str) -> None:
    """Export a Tank + SolverUnit flowsheet and verify the manifest contents."""
    with tempfile.TemporaryDirectory() as output_dir:
        fmu_path = build_fmu(
            fake_dynamic_flowsheet_path, output_dir=output_dir, backend="scipy"
        )
        assert os.path.isfile(fmu_path)

        manifest_path = fmu_path.replace(".fmu", "_fmu_interface.json")
        assert os.path.isfile(manifest_path)

        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        variables = {v["logical_name"]: v for v in manifest["variables"]}
        # SolverUnit ports
        assert "fake_solver.solver_config.temperature" in variables
        assert "fake_solver.k_eff" in variables
        # Tank state outputs are exposed in solverunit_dynamic mode
        assert "state_tank_T" in variables
        assert "state_tank_n_salt" in variables
        # Stream outputs
        assert "out_tank_out_T" in variables


def test_solverunit_dynamic_fmu_simulation(fake_dynamic_flowsheet_path: str) -> None:
    """Run the exported dynamic FMU and check Tank and SolverUnit outputs."""
    pytest.importorskip("pythonfmu")
    pytest.importorskip("fmpy")

    with tempfile.TemporaryDirectory() as output_dir:
        fmu_path = build_fmu(
            fake_dynamic_flowsheet_path, output_dir=output_dir, backend="scipy"
        )

        from fmpy import simulate_fmu

        result = simulate_fmu(
            fmu_path,
            stop_time=5.0,
            output_interval=1.0,
            output=[
                "out_tank_out_T",
                "out_tank_out_flowrate",
                "state_tank_n_salt",
                "out_fake_solver_k_eff",
                "out_fake_solver_power",
            ],
        )

        final = result[-1]
        # SolverUnit outputs are populated each step
        assert abs(final["out_fake_solver_k_eff"] - 1.05) < 1e-9
        assert abs(final["out_fake_solver_power"] - 1_000_000.0) < 1.0
        # Tank holdup changes over the co-simulation
        assert final["state_tank_n_salt"] != pytest.approx(10.0)
        # Outlet flow rate tracks the constant tank outlet_flow parameter
        assert abs(final["out_tank_out_flowrate"] - 0.5) < 1e-9
