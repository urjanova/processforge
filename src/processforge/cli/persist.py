"""Shared run-persistence helpers for the ``pf run`` / ``pf apply`` CLIs.

Both commands build a :class:`~processforge.types.RunManifest`, save it to the
:class:`~processforge.persistence.archive.ProcessStateArchive`, and write a Zarr
copy of the standardized outputs. That logic was previously duplicated (with
slightly different bugs in each copy); it now lives here in one place.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
from typing import TYPE_CHECKING, Any

from loguru import logger

if TYPE_CHECKING:
    from ..persistence.archive import ProcessStateArchive
    from ..types import RunInfo

_RUN_ID_TS_FMT = "%Y%m%dT%H%M%SZ"


def make_run_id() -> str:
    """Build a unique, sortable run id: ``<UTC timestamp>_<6 hex bytes>``."""
    return (
        datetime.datetime.now(datetime.timezone.utc).strftime(_RUN_ID_TS_FMT)
        + "_"
        + os.urandom(3).hex()
    )


def flowsheet_hash(config: dict) -> str:
    """Short (16-hex) sha256 of a flowsheet config for run/index correlation."""
    return hashlib.sha256(
        json.dumps(config, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()[:16]


def persist_run(
    archive: "ProcessStateArchive",
    fs: Any,
    run_id: str,
    results: dict,
    run_info: "RunInfo | dict",
    config: dict,
    base_name: str,
    snapshot_id: str | None = None,
) -> str:
    """Persist a run's manifest + Zarr outputs into *archive*.

    Args:
        archive: The :class:`ProcessStateArchive` for this flowsheet.
        fs: The flowsheet object that produced *results* (must expose
            ``collect_outputs``). Either ``EOFlowsheet`` or ``Flowsheet``.
        run_id: Unique run id (use :func:`make_run_id`).
        results: The solve's stream/unit result dict (the *converged* one).
        run_info: Provenance from
            :func:`processforge.provenance.build_run_info`.
        config: Validated flowsheet config (for the simulation mode).
        base_name: Flowsheet basename used for archive paths.
        snapshot_id: Optional snapshot id to stamp on the manifest (so a run
            can be linked back to the state it was derived from).

    Returns:
        The run id that was persisted.
    """
    mode = config.get("simulation", {}).get("mode", "steady")
    run_info_dump = (
        run_info.model_dump() if hasattr(run_info, "model_dump") else run_info
    )
    manifest = fs.collect_outputs(
        run_id, mode, base_name, provenance=run_info_dump
    )
    if snapshot_id is not None:
        manifest.snapshot_id = snapshot_id

    stream_results = {k: v for k, v in results.items() if not hasattr(v, "fields")}
    archive.save_run(manifest, stream_results=stream_results)

    # Always persist a Zarr copy of the standardized outputs (fields + artifacts)
    # inside the archive, mirroring `pf run`.
    try:
        from ..result import relink_latest_results, save_results_zarr

        run_results_dir = os.path.join(archive.path, "results", run_id)
        save_results_zarr(
            results,
            os.path.join(run_results_dir, "results.zarr"),
            run_info,
        )
        relink_latest_results(archive.path, run_id)
    except Exception as e:
        logger.warning(f"Failed to write results.zarr: {type(e).__name__}: {e}")

    logger.info(f"Run saved    : {run_id}")
    return run_id
