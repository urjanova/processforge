"""S3 storage-options helper.

Builds s3fs storage options (credentials, endpoint, region) from the
container environment for use by :class:`processforge.persistence.artifact_store.ArtifactStore`.
"""

from __future__ import annotations

import os


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
