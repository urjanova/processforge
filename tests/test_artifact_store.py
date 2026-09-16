"""Tests for ArtifactStore S3 upload helpers."""
from __future__ import annotations

import sys
import types

import pytest

from processforge.persistence.artifact_store import ArtifactStore


@pytest.fixture
def fake_s3fs(monkeypatch):
    """Install a fake s3fs module that records put() calls."""
    calls = []

    class FakeFS:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def put(self, local_path: str, uri: str, **kwargs):
            calls.append((local_path, uri))

    fake_s3fs_mod = types.ModuleType("s3fs")
    fake_s3fs_mod.S3FileSystem = FakeFS
    monkeypatch.setitem(sys.modules, "s3fs", fake_s3fs_mod)
    return calls


@pytest.fixture
def clean_s3_env(monkeypatch):
    for var in ("S3_BUCKET", "S3_PREFIX", "S3_ACCESS_KEY", "S3_SECRET_KEY"):
        monkeypatch.delenv(var, raising=False)


def test_persist_archive_uploads_all_files(tmp_path, monkeypatch, fake_s3fs, clean_s3_env):
    monkeypatch.setenv("S3_BUCKET", "test-bucket")
    monkeypatch.setenv("S3_PREFIX", "pf")

    archive = tmp_path / "example.pfarchive"
    runs_dir = archive / "runs"
    runs_dir.mkdir(parents=True)
    (runs_dir / "run_1.json").write_text('{"run_id": "run_1"}')
    snapshots_dir = archive / "snapshots" / "0001_2026-01-01T00:00:00Z"
    snapshots_dir.mkdir(parents=True)
    (snapshots_dir / "x").write_bytes(b"state")
    (archive / "latest_run").write_text("run_1")

    store = ArtifactStore()
    uris = store.persist_archive(str(archive), run_id="run_1", flowsheet_hash="abc123")

    assert len(uris) == 3
    assert all(uri.startswith("s3://test-bucket/pf/abc123/run_1/archive/") for uri in uris)
    assert any(uri.endswith("runs/run_1.json") for uri in uris)
    assert any("snapshots/0001_2026-01-01T00:00:00Z/x" in uri for uri in uris)
    assert any(uri.endswith("latest_run") for uri in uris)

    # put was called for every file.
    assert len(fake_s3fs) == 3


def test_persist_archive_no_bucket_is_noop(tmp_path, monkeypatch, fake_s3fs, clean_s3_env):
    archive = tmp_path / "example.pfarchive"
    archive.mkdir()
    (archive / "file.txt").write_text("data")

    store = ArtifactStore()
    uris = store.persist_archive(str(archive), run_id="run_1", flowsheet_hash="abc123")

    assert uris == []
    assert fake_s3fs == []


def test_upload_file_specific_key(tmp_path, monkeypatch, fake_s3fs, clean_s3_env):
    monkeypatch.setenv("S3_BUCKET", "test-bucket")
    local_file = tmp_path / "report.json"
    local_file.write_text('{"ok": true}')

    store = ArtifactStore()
    uri = store.upload_file(str(local_file), remote_key="pf/hash/run/report.json")

    assert uri == "s3://test-bucket/pf/hash/run/report.json"
    assert len(fake_s3fs) == 1
    assert fake_s3fs[0][1] == uri
