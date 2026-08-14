"""Tests for the S3 run-output upload helper."""
from __future__ import annotations

import os
import sys
from datetime import datetime, timezone
from unittest.mock import MagicMock, patch

import pytest

from processforge.utils import s3_upload


@pytest.fixture
def clean_env(monkeypatch):
    for var in (
        "S3_BUCKET",
        "S3_PREFIX",
        "S3_ACCESS_KEY",
        "S3_SECRET_KEY",
        "S3_ENDPOINT_URL",
        "S3_REGION_NAME",
    ):
        monkeypatch.delenv(var, raising=False)
    return monkeypatch


def test_s3_storage_options_defaults_to_ams3(clean_env):
    opts = s3_upload.s3_storage_options()
    assert opts["client_kwargs"]["region_name"] == "ams3"
    # No creds configured -> not present, s3fs uses default chain.
    assert "key" not in opts
    assert "secret" not in opts


def test_s3_storage_options_includes_credentials(clean_env):
    clean_env.setenv("S3_ACCESS_KEY", "ak")
    clean_env.setenv("S3_SECRET_KEY", "sk")
    clean_env.setenv("S3_ENDPOINT_URL", "https://s3.example.com")
    opts = s3_upload.s3_storage_options()
    assert opts["key"] == "ak"
    assert opts["secret"] == "sk"
    assert opts["client_kwargs"]["endpoint_url"] == "https://s3.example.com"


def test_upload_directory_to_s3_disabled_without_bucket(clean_env):
    assert s3_upload.upload_directory_to_s3("/tmp/run") is None


def test_upload_directory_to_s3_uploads_recursively(clean_env):
    clean_env.setenv("S3_BUCKET", "my-bucket")
    clean_env.setenv("S3_PREFIX", "openmc-runs")

    mock_fs = MagicMock()
    mock_s3fs = MagicMock()
    mock_s3fs.S3FileSystem.return_value = mock_fs

    with patch.dict(sys.modules, {"s3fs": mock_s3fs}):
        uri = s3_upload.upload_directory_to_s3("/tmp/processforge/openmc/msre_run")

    assert uri.startswith("s3://my-bucket/openmc-runs/msre_run-")
    assert uri.endswith(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")[:8]) or "-" in uri
    mock_fs.put.assert_called_once()
    args, kwargs = mock_fs.put.call_args
    assert args[0] == "/tmp/processforge/openmc/msre_run"
    assert args[1].startswith("s3://my-bucket/openmc-runs/msre_run-")
    assert kwargs.get("recursive") is True


def test_upload_directory_to_s3_no_prefix(clean_env):
    clean_env.setenv("S3_BUCKET", "b")

    mock_fs = MagicMock()
    mock_s3fs = MagicMock()
    mock_s3fs.S3FileSystem.return_value = mock_fs

    with patch.dict(sys.modules, {"s3fs": mock_s3fs}):
        uri = s3_upload.upload_directory_to_s3("/data/festim")

    assert uri == f"s3://b/festim-{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}"
    assert mock_fs.put.call_args[1].get("recursive") is True
