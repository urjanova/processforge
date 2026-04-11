from __future__ import annotations

from typing import Literal

from pydantic import BaseModel


class RunRequest(BaseModel):
    flowsheet: dict
    s3_bucket: str
    s3_prefix: str | None = None


class JobStatus(BaseModel):
    job_id: str
    status: Literal["queued", "running", "complete", "failed"]
    s3_urls: dict | None = None
    error: str | None = None
    started_at: str
    completed_at: str | None = None
