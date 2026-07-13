from __future__ import annotations

import os
import threading
from typing import Literal

from pydantic import BaseModel, Field


def _get_default_bucket() -> str:
    return os.environ.get("S3_BUCKET_NAME", "")


def _utcnow_iso() -> str:
    import datetime
    return datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z")


class UnitDef(BaseModel):
    type: str
    params: dict = Field(default_factory=dict)


class StreamDef(BaseModel):
    sources: list[str] = Field(default_factory=list)
    params: dict = Field(default_factory=dict)


class SimulationConfig(BaseModel):
    mode: Literal["steady", "dynamic"] = "steady"
    backend: str | None = None


class FlowsheetConfig(BaseModel):
    units: dict[str, UnitDef]
    streams: dict[str, StreamDef]
    simulation: SimulationConfig = Field(default_factory=SimulationConfig)


class RunRequest(BaseModel):
    flowsheet: FlowsheetConfig
    s3_bucket: str = Field(default_factory=_get_default_bucket)
    s3_prefix: str | None = None


class RunResponse(BaseModel):
    job_id: str
    status: Literal["queued", "running", "complete", "failed"]
    status_url: str


class JobStatus(BaseModel):
    job_id: str
    status: Literal["queued", "running", "complete", "failed"]
    s3_urls: dict | None = None
    error: str | None = None
    started_at: str
    completed_at: str | None = None


class JobStore:
    """Thread-safe in-memory job store."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._jobs: dict[str, JobStatus] = {}

    def create(self, job_id: str) -> JobStatus:
        entry = JobStatus(job_id=job_id, status="queued", started_at=_utcnow_iso())
        with self._lock:
            self._jobs[job_id] = entry
        return entry

    def get(self, job_id: str) -> JobStatus | None:
        with self._lock:
            return self._jobs.get(job_id)

    def update(self, job_id: str, **kwargs) -> None:
        with self._lock:
            entry = self._jobs.get(job_id)
            if entry is not None:
                for k, v in kwargs.items():
                    setattr(entry, k, v)
