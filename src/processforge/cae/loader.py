"""Filesystem helpers for locating and parsing CAE cases and the index."""

from __future__ import annotations

import os
from typing import Optional

import yaml

from .types import CAECase, CAEIndex, CAERequirement, CAEArgument, CAEClaim, CAEEvidence
from .types import CAECaseMeta, CAEEvidenceRef

CAE_DIR_ENV = "PROCESSFORGE_CAE_DIR"
INDEX_FILENAME = "index.yaml"
CASES_DIRNAME = "cases"
CASE_SUFFIX = ".ca.yaml"


class CAEError(Exception):
    """Raised when a CAE case, index, or evidence cannot be resolved."""


def find_cae_root(start: str | None = None) -> Optional[str]:
    """Locate the CAE repository directory.

    Resolution order:

    1. ``$PROCESSFORGE_CAE_DIR`` when set.
    2. A ``cae/`` directory containing ``index.yaml``, walking up from *start*
       (or the current working directory).
    3. *start* itself when it already looks like a CAE directory.
    """
    env_root = os.environ.get(CAE_DIR_ENV)
    if env_root:
        env_root = os.path.abspath(os.path.expanduser(env_root))
        if os.path.isfile(os.path.join(env_root, INDEX_FILENAME)):
            return env_root
        return None

    current = os.path.abspath(os.path.expanduser(start or os.getcwd()))
    if os.path.isfile(os.path.join(current, INDEX_FILENAME)):
        return current

    for directory in _iter_ancestors(current):
        candidate = os.path.join(directory, "cae")
        if os.path.isfile(os.path.join(candidate, INDEX_FILENAME)):
            return candidate
    return None


def _iter_ancestors(path: str):
    current = path
    while True:
        yield current
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent


def load_index(cae_root: str) -> CAEIndex:
    """Read and validate ``<cae_root>/index.yaml``."""
    path = os.path.join(cae_root, INDEX_FILENAME)
    try:
        with open(path, "r", encoding="utf-8") as handle:
            raw = yaml.safe_load(handle) or {}
    except FileNotFoundError as exc:
        raise CAEError(f"CAE index not found at '{path}'.") from exc
    except yaml.YAMLError as exc:
        raise CAEError(f"Failed to parse CAE index '{path}': {exc}") from exc

    try:
        return CAEIndex.model_validate(raw)
    except Exception as exc:  # pydantic ValidationError
        raise CAEError(f"Invalid CAE index '{path}': {exc}") from exc


def load_case(path: str) -> CAECase:
    """Read and validate a single ``*.ca.yaml`` case file."""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            raw = yaml.safe_load(handle) or {}
    except FileNotFoundError as exc:
        raise CAEError(f"CAE case not found at '{path}'.") from exc
    except yaml.YAMLError as exc:
        raise CAEError(f"Failed to parse CAE case '{path}': {exc}") from exc

    try:
        return CAECase.model_validate(raw)
    except Exception as exc:  # pydantic ValidationError
        raise CAEError(f"Invalid CAE case '{path}': {exc}") from exc


def iter_case_paths(cases_dir: str) -> list[str]:
    """Return every case file in *cases_dir*, sorted for determinism."""
    if not os.path.isdir(cases_dir):
        return []
    return sorted(
        os.path.join(cases_dir, name)
        for name in os.listdir(cases_dir)
        if name.endswith(CASE_SUFFIX)
    )


def version_key(version: str) -> tuple:
    """Build a sortable key for a semantic version, tolerating ``v`` prefixes."""
    text = version.strip().lstrip("vV")
    parts: list[object] = []
    for piece in text.split("."):
        try:
            parts.append((0, int(piece)))
        except ValueError:
            parts.append((1, piece))
    return tuple(parts)


def select_version(versions: list[str], requested: str | None) -> Optional[str]:
    """Pick *requested* when present, otherwise the highest version."""
    if not versions:
        return None
    if requested:
        return requested if requested in versions else None
    return max(versions, key=version_key)


__all__ = [
    "CAE_DIR_ENV",
    "CAEError",
    "find_cae_root",
    "iter_case_paths",
    "load_case",
    "load_index",
    "select_version",
    "version_key",
    # Re-exported for convenience so callers need only import from the loader.
    "CAEArgument",
    "CAECase",
    "CAECaseMeta",
    "CAEClaim",
    "CAEEvidence",
    "CAEEvidenceRef",
    "CAEIndex",
    "CAERequirement",
]
