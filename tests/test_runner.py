"""Tests for the high-level processforge runner wrapper."""
from __future__ import annotations

import json
import sys
import types
from pathlib import Path

import pytest

from processforge.runner import (
    FlowsheetValidationError,
    ProviderUnavailableError,
    RunResult,
    apply_flowsheet,
    run_flowsheet,
)

TESTS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = TESTS_DIR.parent
HYDRAULIC_CHAIN = PROJECT_DIR / "flowsheets" / "hydraulic-chain.json"


@pytest.fixture
def hydraulic_flowsheet_dict():
    with open(HYDRAULIC_CHAIN, "r", encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture
def fake_s3fs(monkeypatch):
    """Install a fake s3fs module that records put() calls."""
    calls = []

    class FakeFS:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def put(self, local_path: str, uri: str, **kwargs):
            calls.append((local_path, uri))

        def ls(self, path: str):
            return []

    fake_s3fs_mod = types.ModuleType("s3fs")
    fake_s3fs_mod.S3FileSystem = FakeFS
    monkeypatch.setitem(sys.modules, "s3fs", fake_s3fs_mod)
    return calls


def test_run_flowsheet_hydraulic_chain(tmp_path, monkeypatch):
    monkeypatch.setenv("PROCESSFORGE_OUTPUT_DIR", str(tmp_path))
    result = run_flowsheet(str(HYDRAULIC_CHAIN))

    assert isinstance(result, RunResult)
    assert result.status == "converged"
    assert result.converged is True
    assert result.mode == "steady"
    assert result.run_id
    assert result.flowsheet_hash
    assert result.archive_path == str(tmp_path / "hydraulic-chain.pfarchive")
    assert Path(result.manifest_path).exists()
    assert result.remote_uris == []
    assert result.archive_uris == []


def test_run_flowsheet_with_dict_config(tmp_path, monkeypatch, hydraulic_flowsheet_dict):
    monkeypatch.setenv("PROCESSFORGE_OUTPUT_DIR", str(tmp_path))
    result = run_flowsheet(hydraulic_flowsheet_dict)

    assert result.status == "converged"
    assert result.converged is True
    assert Path(result.archive_path).exists()


def test_run_flowsheet_uploads_archive_to_s3(
    tmp_path, monkeypatch, hydraulic_flowsheet_dict, fake_s3fs
):
    monkeypatch.setenv("PROCESSFORGE_OUTPUT_DIR", str(tmp_path))
    monkeypatch.setenv("S3_BUCKET", "my-bucket")

    result = run_flowsheet(hydraulic_flowsheet_dict)

    assert result.status == "converged"
    assert result.archive_uris
    assert result.remote_uris == result.archive_uris
    assert all(uri.startswith("s3://my-bucket/processforge/") for uri in result.archive_uris)
    assert any("/archive/runs/" in uri for uri in result.archive_uris)
    assert any(uri.endswith(f"runs/{result.run_id}.json") for uri in result.archive_uris)

    # Verify the fake s3fs received the manifest upload.
    manifest_uploads = [
        (local, uri)
        for local, uri in fake_s3fs
        if uri.endswith(f"runs/{result.run_id}.json")
    ]
    assert len(manifest_uploads) == 1
    assert Path(manifest_uploads[0][0]).exists()


def test_apply_flowsheet_up_to_date(tmp_path, monkeypatch, hydraulic_flowsheet_dict):
    monkeypatch.setenv("PROCESSFORGE_OUTPUT_DIR", str(tmp_path))

    # First apply creates a snapshot.
    first = apply_flowsheet(hydraulic_flowsheet_dict)
    assert first.status == "converged"
    assert first.snapshot_id

    # Second apply with the same config detects no drift.
    second = apply_flowsheet(hydraulic_flowsheet_dict)
    assert second.status == "up_to_date"
    assert second.converged is True
    assert second.snapshot_id == first.snapshot_id


def test_apply_flowsheet_missing_file():
    with pytest.raises(FlowsheetValidationError):
        apply_flowsheet("/no/such/file.json")


def test_provider_unavailable_error(tmp_path, monkeypatch):
    import processforge.runner as runner_mod

    flowsheet = {
        "providers": {"openmc": {"type": "openmc", "url": "http://localhost:59999"}},
        "streams": {},
        "units": {},
        "simulation": {"mode": "steady"},
    }
    monkeypatch.setenv("PROCESSFORGE_OUTPUT_DIR", str(tmp_path))
    # Speed up the health-check retry loop.
    monkeypatch.setattr(runner_mod, "HEALTH_MAX_ATTEMPTS", 2)
    monkeypatch.setattr(runner_mod, "HEALTH_RETRY_DELAY", 0)

    with pytest.raises(ProviderUnavailableError):
        run_flowsheet(flowsheet)


def test_run_flowsheet_invalid_config():
    with pytest.raises(FlowsheetValidationError):
        run_flowsheet({"providers": {}, "streams": {}, "units": {}})
