from __future__ import annotations

import os
from typing import Literal

from pydantic import BaseModel, Field


def _get_default_bucket() -> str:
    return os.environ.get("S3_BUCKET_NAME", "")


class RunRequest(BaseModel):
    flowsheet: dict
    s3_bucket: str = Field(default_factory=_get_default_bucket)
    s3_prefix: str | None = None


class JobStatus(BaseModel):
    job_id: str
    status: Literal["queued", "running", "complete", "failed"]
    s3_urls: dict | None = None
    error: str | None = None
    started_at: str
    completed_at: str | None = None
