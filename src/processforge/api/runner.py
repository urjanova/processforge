from __future__ import annotations

import datetime
from pathlib import Path

import boto3
from loguru import logger

from processforge.provenance import build_dynamic_x0, build_run_info
from processforge.result import save_results_zarr_s3
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
            logger.info(f"Uploaded {path} → {uri}")
    return uploaded


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


def run_job_sync(
    job_id: str,
    flowsheet: dict,
    s3_bucket: str,
    s3_prefix: str,
    jobs_store: dict,
) -> None:
    """Run a simulation job synchronously and upload results to S3.

    This is called from an async background task via asyncio.to_thread so
    the FastAPI event loop remains free to handle status-poll requests.
    """
    jobs_store[job_id]["status"] = "running"
    jobs_store[job_id]["started_at"] = datetime.datetime.utcnow().isoformat() + "Z"

    try:
        validate_flowsheet_dict(flowsheet)
        results, run_info = _run_simulation(flowsheet)

        zarr_uri = f"s3://{s3_bucket}/{s3_prefix}/results.zarr"
        save_results_zarr_s3(results, zarr_uri, run_info=run_info)
        logger.info(f"[{job_id}] Zarr written to {zarr_uri}")

        import os
        s3_client = boto3.client(
            "s3",
            aws_access_key_id=os.environ.get("S3_ACCESS_KEY"),
            aws_secret_access_key=os.environ.get("S3_SECRET_KEY"),
            endpoint_url=os.environ.get("S3_ENDPOINT_URL"),
            region_name=os.environ.get("S3_REGION_NAME", "ams3"),
        )
        output_urls = _upload_outputs_dir(s3_client, s3_bucket, s3_prefix)

        jobs_store[job_id].update(
            {
                "status": "complete",
                "completed_at": datetime.datetime.utcnow().isoformat() + "Z",
                "s3_urls": {
                    "zarr": zarr_uri,
                    "outputs": output_urls,
                },
            }
        )
    except Exception as exc:
        logger.exception(f"[{job_id}] Job failed: {exc}")
        jobs_store[job_id].update(
            {
                "status": "failed",
                "completed_at": datetime.datetime.utcnow().isoformat() + "Z",
                "error": str(exc),
            }
        )
