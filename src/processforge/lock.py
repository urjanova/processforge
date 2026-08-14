"""Lock file read/write for .processforge/<env-key>/lock.json.

Each flowsheet gets its own environment directory under ``.processforge/``,
keyed by a stable hash of the flowsheet's absolute path. This lets multiple
flowsheets keep separate provider environments (lock + compose) without
clobbering each other. The lock file is provenance metadata only.
"""
from __future__ import annotations

import hashlib
import json
import os
from typing import Any

LOCK_VERSION = 1
LOCK_FILENAME = "lock.json"


def flowsheet_env_key(flowsheet_path: str) -> str:
    """Stable key for a flowsheet's environment dir.

    Derived from the sha1 of the flowsheet's real (absolute) path so it is
    collision-proof across folders and stable across re-inits.
    """
    abs_path = os.path.realpath(flowsheet_path)
    return hashlib.sha1(abs_path.encode("utf-8")).hexdigest()[:12]


def flowsheet_env_dir(pf_dir: str, flowsheet_path: str) -> str:
    """Directory holding a flowsheet's lock.json and docker-compose.yml."""
    return os.path.join(pf_dir, flowsheet_env_key(flowsheet_path))


def _lock_path(pf_dir: str, flowsheet: str | None = None) -> str:
    if flowsheet is not None:
        return os.path.join(flowsheet_env_dir(pf_dir, flowsheet), LOCK_FILENAME)
    return os.path.join(pf_dir, LOCK_FILENAME)


def read_lock(pf_dir: str, flowsheet: str | None = None) -> dict[str, Any] | None:
    """Read lock.json.

    When *flowsheet* is given, reads that flowsheet's env dir. Otherwise falls
    back to the legacy root ``.processforge/lock.json`` (used during migration).
    Returns None if not found.
    """
    path = _lock_path(pf_dir, flowsheet)
    if not os.path.exists(path):
        if flowsheet is None and os.path.exists(os.path.join(pf_dir, LOCK_FILENAME)):
            path = os.path.join(pf_dir, LOCK_FILENAME)
        else:
            return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_lock(
    pf_dir: str,
    flowsheet: str,
    providers: dict[str, dict],
    processforge_version: str,
) -> None:
    """Write lock.json into the flowsheet's environment directory."""
    lock = {
        "version": LOCK_VERSION,
        "processforge_version": processforge_version,
        "flowsheet": flowsheet,
        "providers": providers,
    }
    env_dir = flowsheet_env_dir(pf_dir, flowsheet)
    os.makedirs(env_dir, exist_ok=True)
    path = os.path.join(env_dir, LOCK_FILENAME)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(lock, f, indent=2)
