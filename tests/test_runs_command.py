"""Tests for per-run Zarr outputs and the ``pf runs`` listing command."""
import json
import os

import zarr
from typer.testing import CliRunner

import typer

from processforge.cli.runs import runs
from processforge.result import relink_latest_results, save_results_zarr
from processforge.types import (
    EngineOutput,
    OutputArtifact,
    OutputField,
    Quantity,
    RunManifest,
)


def _eo(run_id_suffix: str, keff: float) -> EngineOutput:
    return EngineOutput(
        engine="openmc",
        sim_type="eigenvalue_reactor",
        status="completed",
        fields=[
            OutputField(
                name="k_eff",
                quantity=Quantity(value=keff, unit=""),
                kind="scalar",
                source="keff",
            )
        ],
        artifacts=[
            OutputArtifact(
                name="statepoint",
                kind="h5",
                local_path=f"/data/openmc/msre_run/{run_id_suffix}/statepoint.h5",
                remote_uris=[],
                source="local",
            )
        ],
    )


def test_relink_latest_results_symlinks_to_latest(tmp_path):
    d = str(tmp_path)
    rid_a = "20260101T000000Z_aaaaaa"
    rid_b = "20260102T000000Z_bbbbbb"

    save_results_zarr({"openmc_solver": _eo("a", 1.04)},
                      os.path.join(d, "results", rid_a, "results.zarr"), None)
    relink_latest_results(d, rid_a)
    link = os.path.join(d, "results.zarr")
    assert os.path.islink(link)
    assert os.path.realpath(link) == os.path.realpath(
        os.path.join(d, "results", rid_a, "results.zarr"))

    # A second run re-points the symlink without clobbering run A's zarr.
    save_results_zarr({"openmc_solver": _eo("b", 1.06)},
                      os.path.join(d, "results", rid_b, "results.zarr"), None)
    relink_latest_results(d, rid_b)
    assert os.path.realpath(link) == os.path.realpath(
        os.path.join(d, "results", rid_b, "results.zarr"))
    assert os.path.isdir(os.path.join(d, "results", rid_a, "results.zarr"))

    root = zarr.open_group(link, mode="r")
    assert abs(root["openmc_solver"]["k_eff"][0] - 1.06) < 1e-9


def test_pf_runs_lists_runs_with_disk_status(tmp_path, monkeypatch):
    monkeypatch.setenv("PROCESSFORGE_OUTPUT_DIR", str(tmp_path))
    base = "test_flowsheet"
    archive_path = tmp_path / f"{base}.pfarchive"
    runs_dir = archive_path / "runs"
    results_dir = archive_path / "results"
    runs_dir.mkdir(parents=True)

    rid_a = "20260101T000000Z_aaaaaa"
    rid_b = "20260102T000000Z_bbbbbb"

    # Run A has its zarr on disk; run B does not (deleted).
    save_results_zarr({"openmc_solver": _eo("a", 1.04)},
                      str(results_dir / rid_a / "results.zarr"), None)
    for rid, keff in ((rid_a, 1.04), (rid_b, 1.06)):
        manifest = RunManifest(
            run_id=rid, timestamp=rid.split("_")[0], mode="steady",
            units={"openmc_solver": _eo(rid[-6:], keff)},
        )
        (runs_dir / f"{rid}.json").write_text(manifest.model_dump_json())
    (archive_path / "latest_run").write_text(rid_b)

    app = typer.Typer()
    app.command()(runs)
    result = CliRunner().invoke(app, [f"/tmp/{base}.json"])
    assert result.exit_code == 0, result.output
    assert rid_a in result.output
    assert rid_b in result.output
    assert "*" in result.output  # latest marker on run B
    assert "LATEST" in result.output
    assert "RUN ID" in result.output


def test_pf_runs_shows_manifest_and_artifacts(tmp_path, monkeypatch):
    monkeypatch.setenv("PROCESSFORGE_OUTPUT_DIR", str(tmp_path))
    base = "test_flowsheet"
    archive_path = tmp_path / f"{base}.pfarchive"
    runs_dir = archive_path / "runs"
    runs_dir.mkdir(parents=True)

    rid = "20260101T000000Z_aaaaaa"
    manifest = RunManifest(
        run_id=rid, timestamp=rid.split("_")[0], mode="steady",
        units={"openmc_solver": _eo(rid[-6:], 1.04)},
    )
    (runs_dir / f"{rid}.json").write_text(manifest.model_dump_json())

    app = typer.Typer()
    app.command()(runs)
    result = CliRunner().invoke(app, [f"/tmp/{base}.json", rid])
    assert result.exit_code == 0, result.output
    assert "k_eff" in result.output
    assert "Artifacts on disk:" in result.output
