from __future__ import annotations

import asyncio
import datetime
import uuid

from fastapi import BackgroundTasks, FastAPI, HTTPException

from processforge.api.models import JobStatus, RunRequest
from processforge.api.runner import run_job_sync

try:
    from processforge import __version__ as _version
except Exception:
    _version = "unknown"

app = FastAPI(title="ProcessForge API", version=_version)

_jobs: dict[str, dict] = {}


@app.on_event("startup")
async def startup_event():
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


async def _run_job_in_thread(
    job_id: str,
    flowsheet: dict,
    s3_bucket: str,
    s3_prefix: str,
) -> None:
    await asyncio.to_thread(run_job_sync, job_id, flowsheet, s3_bucket, s3_prefix, _jobs)


@app.get("/api/v1/health")
def health():
    return {"status": "ok", "version": _version}


@app.post("/api/v1/run", status_code=202)
async def run_flowsheet(request: RunRequest, background_tasks: BackgroundTasks):
    job_id = str(uuid.uuid4())
    prefix = request.s3_prefix or f"jobs/{job_id}"

    _jobs[job_id] = {
        "job_id": job_id,
        "status": "queued",
        "started_at": datetime.datetime.utcnow().isoformat() + "Z",
        "completed_at": None,
        "s3_urls": None,
        "error": None,
    }

    background_tasks.add_task(
        _run_job_in_thread,
        job_id,
        request.flowsheet,
        request.s3_bucket,
        prefix,
    )

    return {
        "job_id": job_id,
        "status": "queued",
        "status_url": f"/api/v1/jobs/{job_id}",
    }


@app.get("/api/v1/jobs/{job_id}", response_model=JobStatus)
def get_job(job_id: str):
    job = _jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return job
