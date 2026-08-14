"""S3 upload helper for provider containers.

Uploads a provider's run-output directory (``run_dir``) to S3 so ephemeral
container-local artifacts (``statepoint.h5``, mesh tallies, ``*.h5m``) survive
after the container stops. The container opts in by setting ``S3_BUCKET`` in its
environment; credential env vars reuse the convention from
:func:`processforge.result.save_results_zarr_s3`.
"""

from __future__ import annotations

import os
from datetime import datetime, timezone


def s3_storage_options() -> dict:
    """Build s3fs storage options from environment variables.

    Returns an empty dict when no credentials are configured so s3fs falls
    back to its default credential chain (e.g. instance/role-based auth).
    """
    opts: dict = {}
    key = os.environ.get("S3_ACCESS_KEY")
    secret = os.environ.get("S3_SECRET_KEY")
    if key:
        opts["key"] = key
    if secret:
        opts["secret"] = secret
    client_kwargs = {
        k: v
        for k, v in {
            "endpoint_url": os.environ.get("S3_ENDPOINT_URL"),
            "region_name": os.environ.get("S3_REGION_NAME", "ams3"),
        }.items()
        if v
    }
    if client_kwargs:
        opts["client_kwargs"] = client_kwargs
    return opts


def upload_directory_to_s3(run_dir: str) -> str | None:
    """Upload *run_dir* recursively to S3.

    Reads the destination from the container environment:

    * ``S3_BUCKET``  — required; enables the upload. When unset, returns ``None``
      and the caller keeps its local artifacts.
    * ``S3_PREFIX``  — optional key prefix (e.g. ``openmc-runs``).

    The object key is ``<S3_PREFIX>/<run_dir_basename>-<utc_stamp>/`` so
    repeated runs do not clobber each other.

    Returns the uploaded ``s3://`` URI, or ``None`` when disabled.
    """
    bucket = os.environ.get("S3_BUCKET")
    if not bucket:
        return None

    prefix = os.environ.get("S3_PREFIX", "").strip("/")
    tag = f"{os.path.basename(run_dir.rstrip('/'))}-{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}"
    s3_dir = f"{prefix + '/' if prefix else ''}{tag}"

    import s3fs

    s3fs.S3FileSystem(**s3_storage_options()).put(
        run_dir, f"s3://{bucket}/{s3_dir}", recursive=True
    )
    return f"s3://{bucket}/{s3_dir}"
