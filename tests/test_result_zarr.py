"""Tests for ``save_results_zarr`` handling of EngineOutput / StreamOutput."""
import json
import os
import tempfile

import zarr

from processforge.result import save_results_zarr
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
