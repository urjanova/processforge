"""Tests for ``save_results_zarr`` handling of EngineOutput / StreamOutput."""
import json
import os
import tempfile
import warnings

import zarr

import numpy as np

from processforge.result import _convert_value, save_results_zarr, summarize_zarr_store
from processforge.types import (
    EngineOutput,
    OutputArtifact,
    OutputField,
    Quantity,
)


def _make_openmc_eo(remote: bool) -> EngineOutput:
    art = OutputArtifact(
        name="statepoint",
        kind="h5",
        local_path=None if remote else "/data/openmc/msre_run/statepoint.h5",
        remote_uris=["s3://bucket/prefix/statepoint.h5"] if remote else [],
        source="remote" if remote else "local",
    )
    return EngineOutput(
        engine="openmc",
        sim_type="eigenvalue_reactor",
        status="completed",
        fields=[
            OutputField(
                name="k_eff",
                quantity=Quantity(value=1.05, unit=""),
                kind="scalar",
                source="keff",
            )
        ],
        artifacts=[art],
    )


def test_save_results_zarr_engine_output_remote_artifact():
    eo = _make_openmc_eo(remote=True)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "results.zarr")
        save_results_zarr({"openmc_solver": eo}, path)

        store = zarr.storage.LocalStore(path)
        root = zarr.open_group(store=store, mode="r")

        assert "openmc_solver" in root
        g = root["openmc_solver"]
        assert g.attrs["engine"] == "openmc"
        assert "k_eff" in g.array_keys()
        assert g["k_eff"].attrs["source"] == "keff"

        assert "artifacts" in g
        assert "statepoint" in g["artifacts"].attrs
        dumped = json.loads(g["artifacts"].attrs["statepoint"])
        assert dumped["remote_uris"] == ["s3://bucket/prefix/statepoint.h5"]


def test_save_results_zarr_engine_output_local_artifact():
    eo = _make_openmc_eo(remote=False)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "results.zarr")
        save_results_zarr({"openmc_solver": eo}, path)

        store = zarr.storage.LocalStore(path)
        root = zarr.open_group(store=store, mode="r")
        g = root["openmc_solver"]
        dumped = json.loads(g["artifacts"].attrs["statepoint"])
        assert dumped["local_path"] == "/data/openmc/msre_run/statepoint.h5"
        assert dumped["source"] == "local"


def test_save_results_zarr_string_fields_use_stable_utf8():
    """String stream fields (e.g. 'unit') must not use unstable FixedLengthUTF32."""
    results = {
        "after_pipe_1": {
            "T": 313.0,
            "P": 101325.0,
            "flowrate": 1.0,
            "z": {"Water": 1.0},
            "unit": "Pipes",
        },
        "after_low_pressure_pump": {
            "T": 313.0,
            "P": 102325.0,
            "flowrate": 1.0,
            "z": {"Water": 1.0},
            "unit": "Pump",
        },
    }
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "results.zarr")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            save_results_zarr(results, path)

        unstable = [
            w for w in caught
            if "FixedLengthUTF32" in str(w.message)
            and w.category.__name__ == "UnstableSpecificationWarning"
        ]
        assert not unstable

        store = zarr.storage.LocalStore(path)
        root = zarr.open_group(store=store, mode="r")

        assert root["after_pipe_1"]["unit"][:] == "Pipes"
        assert root["after_low_pressure_pump"]["unit"][:] == "Pump"

        schema_path = path + ".schema.json"
        with open(schema_path, encoding="utf-8") as f:
            schema = json.load(f)
        assert schema["streams"]["after_pipe_1"]["dtypes"]["unit"] == "str"
        assert schema["streams"]["after_low_pressure_pump"]["dtypes"]["unit"] == "str"


def test_convert_value_empty_array():
    assert _convert_value(np.array([])) == ""


def test_convert_value_scalar_array():
    assert _convert_value(np.array([3.14])) == 3.14
    assert _convert_value(np.array(2.0)) == 2.0


def test_convert_value_multi_element_array():
    arr = np.array([1.0, 2.0, 3.0])
    assert _convert_value(arr) == [1.0, 2.0, 3.0]


def test_convert_value_nested_arrays():
    arr = np.array([[1.0, 2.0], [3.0, 4.0]])
    assert _convert_value(arr) == [[1.0, 2.0], [3.0, 4.0]]


def test_summarize_zarr_store_handles_multi_element_stream():
    results = {
        "after_drain_valve_1": {
            "T": np.array([313.0, 313.0, 313.0]),
            "P": np.array([50162.5, 50162.5, 50162.5]),
            "flowrate": np.array([1.0, 1.0, 1.0]),
            "z": {"Water": np.array([1.0, 1.0, 1.0])},
        }
    }
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "results.zarr")
        save_results_zarr(results, path)
        summary = summarize_zarr_store(path)
        assert summary["present"]
        stream = summary["streams"]["after_drain_valve_1"]
        assert stream["fields"]["P"]["value"] == [50162.5, 50162.5, 50162.5]
        assert stream["fields"]["T"]["value"] == [313.0, 313.0, 313.0]
