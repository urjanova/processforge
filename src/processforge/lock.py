"""Lock file read/write for .processforge/lock.json."""
from __future__ import annotations

import json
import os
from typing import Any

LOCK_VERSION = 1
LOCK_FILENAME = "lock.json"


def _lock_path(pf_dir: str) -> str:
    return os.path.join(pf_dir, LOCK_FILENAME)


def read_lock(pf_dir: str) -> dict[str, Any] | None:
    """Read lock.json from .processforge/. Returns None if not found."""
    path = _lock_path(pf_dir)
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_lock(
    pf_dir: str,
    flowsheet: str,
    providers: dict[str, dict],
    processforge_version: str,
) -> None:
    """Write lock.json to .processforge/."""
    lock = {
        "version": LOCK_VERSION,
        "processforge_version": processforge_version,
        "flowsheet": flowsheet,
        "providers": providers,
    }
    os.makedirs(pf_dir, exist_ok=True)
    path = _lock_path(pf_dir)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(lock, f, indent=2)
