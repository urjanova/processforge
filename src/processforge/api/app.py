from __future__ import annotations

import asyncio
import uuid
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException
from loguru import logger

from processforge.api.models import JobStatus, JobStore, RunRequest, RunResponse
from processforge.api.runner import run_job

try:
    from processforge import __version__ as _version
except Exception:
    _version = "unknown"

job_store = JobStore()


def _utcnow_iso() -> str:
    import datetime
    return datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z")


@asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncIterator[None]:
    import os
    required_env_vars = [
        "S3_ACCESS_KEY",
        "S3_SECRET_KEY",
        "S3_ENDPOINT_URL",
        "S3_BUCKET_NAME",
    ]
    missing = [v for v in required_env_vars if not os.environ.get(v)]
    if missing:
        raise RuntimeError(f"Missing required S3 environment variables: {', '.join(missing)}")
    yield


app = FastAPI(title="ProcessForge API", version=_version, lifespan=lifespan)


@app.middleware("http")
async def log_requests(request, call_next):
    logger.debug(f"{request.method} {request.url.path}")
    response = await call_next(request)
    logger.debug(f"{request.method} {request.url.path} -> {response.status_code}")
    return response


@app.get("/api/v1/health")
def health():
    return {"status": "ok", "version": _version}


@app.post("/api/v1/run", status_code=202, response_model=RunResponse)
async def run_flowsheet(request: RunRequest):
    job_id = str(uuid.uuid4())
    prefix = request.s3_prefix or f"jobs/{job_id}"

    job_store.create(job_id)

    asyncio.create_task(
        run_job(
            job_id,
            request.flowsheet.model_dump(),
            request.s3_bucket,
            prefix,
            job_store,
        )
    )

    return RunResponse(
        job_id=job_id,
        status="queued",
        status_url=f"/api/v1/jobs/{job_id}",
    )


@app.get("/api/v1/jobs/{job_id}", response_model=JobStatus)
def get_job(job_id: str):
    job = job_store.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return job
