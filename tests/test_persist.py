"""Tests for ``cli/persist.py`` — the shared run-persistence helper."""
from __future__ import annotations

from unittest.mock import MagicMock

from processforge.cli.persist import flowsheet_hash, make_run_id, persist_run
from processforge.persistence.archive import ProcessStateArchive
from processforge.types import RunManifest


def _make_fs(run_id: str) -> MagicMock:
    fs = MagicMock()
    fs.collect_outputs.return_value = RunManifest(run_id=run_id)
    return fs


def test_make_run_id_format_and_unique():
    rid = make_run_id()
    assert len(rid) == len("20260101T000000Z_") + 6
    assert "_" in rid
    assert make_run_id() != make_run_id()


def test_flowsheet_hash_is_deterministic_and_short():
    cfg = {"a": 1, "b": [1, 2]}
    h1 = flowsheet_hash(cfg)
    h2 = flowsheet_hash(dict(cfg))
    assert h1 == h2
    assert len(h1) == 16


def test_persist_run_saves_manifest_and_zarr(tmp_path):
    archive = ProcessStateArchive(str(tmp_path / "fs.pfarchive"))
    run_id = make_run_id()
    fs = _make_fs(run_id)
    results = {"s1": {"T": [300.0, 310.0], "P": [101325.0]}}
    config = {"simulation": {"mode": "steady"}}

    returned = persist_run(
        archive, fs, run_id, results, {"mode": "steady", "backend": "scipy"},
        config, "fs",
    )

    assert returned == run_id
    # Manifest was persisted and is loadable.
    loaded = archive.load_run(run_id)
    assert loaded is not None
    assert loaded.run_id == run_id
    # Zarr copy written beside the run.
    import os

    assert os.path.exists(
        os.path.join(archive.path, "results", run_id, "results.zarr")
    )
    # collect_outputs was called with the right args.
    fs.collect_outputs.assert_called_once_with(
        run_id, "steady", "fs", provenance={"mode": "steady", "backend": "scipy"}
    )


def test_persist_run_stamps_snapshot_id(tmp_path):
    archive = ProcessStateArchive(str(tmp_path / "fs.pfarchive"))
    run_id = make_run_id()
    fs = _make_fs(run_id)
    results = {"s1": {"T": [300.0]}}

    persist_run(
        archive, fs, run_id, results, {},
        {"simulation": {"mode": "steady"}}, "fs", snapshot_id="snap-123",
    )

    loaded = archive.load_run(run_id)
    assert loaded.snapshot_id == "snap-123"
