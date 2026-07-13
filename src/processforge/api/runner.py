from __future__ import annotations

import asyncio
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

import boto3
from loguru import logger

from processforge.api.models import JobStore
from processforge.provenance import build_dynamic_x0, build_run_info
from processforge.result import save_results_zarr
from processforge.utils.validate_flowsheet import validate_flowsheet_dict


def _upload_outputs_dir(s3_client, bucket: str, prefix: str) -> list[str]:
    """Walk the outputs/ directory and upload every file to S3.

    Returns a list of s3:// URIs for all uploaded files.
    """
    uploaded: list[str] = []
    outputs = Path("outputs")
    if not outputs.exists():
        return uploaded
    for path in sorted(outputs.rglob("*")):
        if path.is_file():
            s3_key = f"{prefix}/outputs/{path.relative_to(outputs)}"
            s3_client.upload_file(str(path), bucket, s3_key)
            uri = f"s3://{bucket}/{s3_key}"
            uploaded.append(uri)
            logger.info("Uploaded {path} -> {uri}", path=path, uri=uri)
    return uploaded


def _build_s3_client():
    return boto3.client(
        "s3",
        aws_access_key_id=os.environ.get("S3_ACCESS_KEY"),
        aws_secret_access_key=os.environ.get("S3_SECRET_KEY"),
        endpoint_url=os.environ.get("S3_ENDPOINT_URL"),
        region_name=os.environ.get("S3_REGION_NAME", "ams3"),
    )


def _run_simulation(config: dict):
    mode = config.get("simulation", {}).get("mode", "steady")
    if mode == "dynamic":
        from processforge.flowsheet import Flowsheet

        fs = Flowsheet(config)
        results = fs.run()
        x0, var_names = build_dynamic_x0(config)
        run_info = build_run_info(config, x0=x0, var_names=var_names)
    else:
        from processforge.eo.flowsheet import EOFlowsheet

        fs = EOFlowsheet(config, backend=None)
        results = fs.run()
        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
    return results, run_info


async def run_job(
    job_id: str,
    flowsheet: dict,
    s3_bucket: str,
    s3_prefix: str,
    job_store: JobStore,
) -> None:
    """Run a simulation job and upload results to S3.

    Runs synchronously via asyncio.to_thread so the FastAPI event loop
    remains free to handle status-poll requests.
    """
    log = logger.bind(job_id=job_id)
    await asyncio.to_thread(
        _run_job_sync, job_id, flowsheet, s3_bucket, s3_prefix, job_store, log
    )


def _run_job_sync(
    job_id: str,
    flowsheet: dict,
    s3_bucket: str,
    s3_prefix: str,
    job_store: JobStore,
    log: Any,
) -> None:
    import datetime

    log_file = tempfile.NamedTemporaryFile(delete=False, suffix=".log")
    log_file.close()

    logger_id = logger.add(log_file.name)

    job_store.update(
        job_id,
        status="running",
        started_at=datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
    )

    try:
        validate_flowsheet_dict(flowsheet)
        results, run_info = _run_simulation(flowsheet)

        s3_client = _build_s3_client()

        with tempfile.TemporaryDirectory() as tmpdir:
            local_zarr_path = os.path.join(tmpdir, "results.zarr")
            save_results_zarr(results, local_zarr_path, run_info=run_info)

            zip_base_name = os.path.join(tmpdir, "results")
            shutil.make_archive(zip_base_name, "zip", root_dir=tmpdir, base_dir="results.zarr")
            zip_path = f"{zip_base_name}.zip"

            zip_s3_key = f"{s3_prefix}/results.zip"
            s3_client.upload_file(zip_path, s3_bucket, zip_s3_key)
            zarr_uri = f"s3://{s3_bucket}/{zip_s3_key}"

        log.info("Zipped Zarr written to {zarr_uri}", zarr_uri=zarr_uri)
        output_urls = _upload_outputs_dir(s3_client, s3_bucket, s3_prefix)

        logger.remove(logger_id)

        log_s3_key = f"{s3_prefix}/run.log"
        s3_client.upload_file(log_file.name, s3_bucket, log_s3_key)
        log_uri = f"s3://{s3_bucket}/{log_s3_key}"
        os.unlink(log_file.name)

        job_store.update(
            job_id,
            status="complete",
            completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
            s3_urls={"zarr": zarr_uri, "outputs": output_urls, "log": log_uri},
        )
    except Exception as exc:
        log.exception("Job failed: {exc}", exc=exc)
        logger.remove(logger_id)

        job_store.update(
            job_id,
            status="failed",
            completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
            error=str(exc),
        )

        try:
            s3_client = _build_s3_client()
            log_s3_key = f"{s3_prefix}/run.log"
            s3_client.upload_file(log_file.name, s3_bucket, log_s3_key)
            log_uri = f"s3://{s3_bucket}/{log_s3_key}"
            job = job_store.get(job_id)
            if job is not None:
                s3_urls = job.s3_urls or {}
                s3_urls["log"] = log_uri
                job_store.update(job_id, s3_urls=s3_urls)
        except Exception as e:
            log.error("Failed to upload log file: {e}", e=e)
        finally:
            if os.path.exists(log_file.name):
                os.unlink(log_file.name)
