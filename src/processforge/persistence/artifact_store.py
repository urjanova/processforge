"""Content-addressed artifact storage for simulation outputs.

Replaces the old ``upload_directory_to_s3`` whole-directory dump.  Providers
declare their output files as :class:`~processforge.types.OutputArtifact`
records (with ``local_path``); the :class:`ArtifactStore` computes a content
hash, uploads each to object storage under a deterministic, reproducible key,
and fills in ``remote_uris``.  The CLI and the container server both use this
same store so artifact handling is identical everywhere.
"""
from __future__ import annotations

import hashlib
import os
from datetime import datetime, timezone
from typing import Optional

from loguru import logger

from ..types import OutputArtifact
from ..utils.s3_upload import s3_storage_options


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _artifact_kind(path: str) -> str:
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    return {
        "h5": "h5",
        "xdmf": "xdmf",
        "csv": "csv",
        "vtk": "vtk",
        "vtu": "vtk",
        "png": "png",
        "zarr": "zarr",
        "nc": "nc",
    }.get(ext, ext or "data")


class ArtifactStore:
    """Register and upload simulation artifact files to object storage."""

    def __init__(self, bucket: Optional[str] = None, prefix: str = "processforge"):
        self.bucket = bucket or os.environ.get("S3_BUCKET")
        self.prefix = prefix

    # ------------------------------------------------------------------
    # Registration
    # ------------------------------------------------------------------
    def register_file(
        self,
        path: str,
        *,
        name: Optional[str] = None,
        kind: Optional[str] = None,
        run_id: str = "",
        flowsheet_hash: str = "",
        unit: str = "",
        content_type: str = "",
    ) -> OutputArtifact:
        """Build an :class:`OutputArtifact` with upload to S3 (if configured)."""
        if not os.path.exists(path):
            raise FileNotFoundError(f"Artifact file not found: {path}")
        checksum = _sha256_file(path)
        size = os.path.getsize(path)
        artifact = OutputArtifact(
            name=name or os.path.basename(path),
            kind=kind or _artifact_kind(path),
            local_path=os.path.abspath(path),
            size_bytes=size,
            checksum=checksum,
            content_type=content_type,
            artifact_id=checksum[:16],
            source="local",
        )
        return self.upload(artifact, run_id=run_id, flowsheet_hash=flowsheet_hash, unit=unit)

    # ------------------------------------------------------------------
    # Upload
    # ------------------------------------------------------------------
    def upload(
        self,
        artifact: OutputArtifact,
        *,
        run_id: str = "",
        flowsheet_hash: str = "",
        unit: str = "",
    ) -> OutputArtifact:
        """Upload *artifact* to S3 (if configured) and return it with URIs.

        No-op (returns the artifact unchanged) when ``S3_BUCKET`` is unset.
        """
        if not self.bucket or not artifact.local_path:
            return artifact

        try:
            import s3fs

            ext = os.path.splitext(artifact.local_path)[1]
            key_parts = [self.prefix, flowsheet_hash, run_id, unit,
                         f"{artifact.name}-{artifact.artifact_id}{ext}"]
            key = "/".join(p for p in key_parts if p)
            uri = f"s3://{self.bucket}/{key}"

            fs = s3fs.S3FileSystem(**s3_storage_options())
            fs.put(artifact.local_path, uri)
            if uri not in artifact.remote_uris:
                artifact.remote_uris.append(uri)
            logger.info(f"Uploaded artifact '{artifact.name}' to {uri}")
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"S3 upload of '{artifact.name}' failed: {exc}")
        return artifact

    # ------------------------------------------------------------------
    # Container-side convenience: persist all artifacts in engine outputs
    # ------------------------------------------------------------------
    def persist_outputs(self, outputs, *, run_id: str = "", flowsheet_hash: str = "", unit: str = ""):
        """Upload every artifact referenced by *outputs*.

        Mutates the artifacts in place, filling ``remote_uris``.  ``outputs``
        may be a single object or an iterable of them.
        """
        items = outputs if isinstance(outputs, (list, tuple)) else [outputs]
        for out in items:
            for art in getattr(out, "artifacts", []):
                if art.local_path and (self.bucket or os.environ.get("S3_BUCKET")):
                    self.upload(art, run_id=run_id, flowsheet_hash=flowsheet_hash, unit=unit)
