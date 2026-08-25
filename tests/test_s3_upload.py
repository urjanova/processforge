"""Tests for the S3 storage-options helper."""
from __future__ import annotations

import os

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


def test_validate_s3_no_bucket_is_noop(clean_env):
    # With S3_BUCKET unset, validate_s3 must not raise (uploads disabled).
    assert s3_upload.validate_s3() is None


def test_validate_s3_unreachable_raises(clean_env, monkeypatch):
    import sys
    import types

    clean_env.setenv("S3_BUCKET", "my-bucket")

    class FakeFS:
        def ls(self, bucket):  # noqa: ANN001, D401
            raise ConnectionError("cannot reach bucket")

    fake_s3fs = types.ModuleType("s3fs")
    fake_s3fs.S3FileSystem = lambda **kwargs: FakeFS()  # noqa: ARG005
    monkeypatch.setitem(sys.modules, "s3fs", fake_s3fs)

    with pytest.raises(RuntimeError, match="unreachable"):
        s3_upload.validate_s3()

