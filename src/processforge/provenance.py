"""Provenance and metadata utilities for Processforge runs."""
from __future__ import annotations

import platform
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np

from .types import RunInfo

_KEY_PACKAGES = [
    "numpy",
    "scipy",
    "zarr",
    "coolprop",
    "loguru",
    "pandas",
    "openpyxl",
    "matplotlib",
    "jsonschema",
    "graphviz",
]


def _git_hash() -> str:
    """Return the current git commit hash, or 'unknown' if not in a git repo."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return "unknown"


def _package_versions(packages: list[str]) -> dict[str, str]:
    """Return installed versions for each package (or 'unknown')."""
    from importlib.metadata import PackageNotFoundError, version

    versions: dict[str, str] = {}
    for pkg in packages:
        try:
            versions[pkg] = version(pkg)
        except PackageNotFoundError:
            versions[pkg] = "unknown"
    return versions


def build_run_info(
    config: dict,
    x0: np.ndarray | None = None,
    var_names: list[str] | None = None,
) -> RunInfo:
    """Build a run_info metadata dict for provenance tracking.

    Args:
        config:    The validated flowsheet configuration dict.
        x0:        Initial guess vector (EO mode) or flattened feed-stream
                   initial conditions (dynamic mode).
        var_names: Human-readable label for each element of x0 (e.g.
                   ``"feed/T"``, ``"feed/P"``, ``"product/z_H2O"``).

    Returns:
        A :class:`RunInfo` instance suitable for passing to
        ``save_results_zarr`` as *run_info*.
    """
    from . import __version__ as _pf_version

    sim_cfg = config.get("simulation", {})
    x0_list = None
    if x0 is not None:
        x0_list = np.asarray(x0, dtype=float).tolist()

    return RunInfo(
        git_hash=_git_hash(),
        timestamp=datetime.now(timezone.utc).isoformat(),
        python_version=sys.version,
        platform=platform.platform(),
        processforge_version=_pf_version,
        mode=sim_cfg.get("mode", "steady"),
        backend=sim_cfg.get("backend", "scipy"),
        pkg_versions=_package_versions(_KEY_PACKAGES),
        x0=x0_list,
        var_names=list(var_names) if var_names is not None else None,
    )


def build_dynamic_x0(config: dict) -> tuple[np.ndarray, list[str]]:
    """Flatten feed-stream initial conditions into a reproducibility vector.

    For dynamic simulations there is no single ``x0`` passed to a global
    solver; instead the starting state is fully determined by the feed streams
    in the config.  This helper packs those values into a flat numpy array so
    the same schema used for EO runs can be applied here.

    Returns:
        (x0_array, var_names) where ``x0_array[i]`` corresponds to
        ``var_names[i]`` (e.g. ``"feed/T"``, ``"feed/z_H2O"``).
    """
    values: list[float] = []
    names: list[str] = []

    for stream_name, stream in config.get("streams", {}).items():
        for prop in ("T", "P", "flowrate"):
            val = stream.get(prop, 0.0)
            values.append(float(val))
            names.append(f"{stream_name}/{prop}")
        for comp, frac in sorted(stream.get("z", {}).items()):
            values.append(float(frac))
            names.append(f"{stream_name}/z_{comp}")

    return np.asarray(values, dtype=float), names
