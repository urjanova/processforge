"""Shared CLI helpers — provider checks, state loading, metadata, divergence reports."""

from __future__ import annotations

import datetime
import hashlib
import importlib
import json
import os
import urllib.error
import urllib.request
from typing import TYPE_CHECKING, Any

from loguru import logger

from .. import __version__ as _pf_version
from ..utils.validate_flowsheet import validate_flowsheet
from ..state import StateManager

if TYPE_CHECKING:
    from ..state import SnapshotState


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

def output_root() -> str:
    """Root directory for run outputs (zarr, provider artifacts).

    Defaults to ``outputs`` for local runs; the Docker image sets
    ``PROCESSFORGE_OUTPUT_DIR=/data`` so outputs land on the mounted volume.
    """
    return os.environ.get("PROCESSFORGE_OUTPUT_DIR", "outputs")


def require_existing_file(path: str, label: str = "Flowsheet file") -> None:
    """Fail fast with a consistent message when an input path is missing."""
    if not os.path.exists(path):
        logger.error(f"{label} '{path}' not found.")
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# Flowsheet validation
# ---------------------------------------------------------------------------

def validate_runtime_flowsheet(path: str) -> dict:
    """Validate flowsheet config and attach source path for runtime providers."""
    try:
        config = validate_flowsheet(path)
    except Exception as e:
        logger.error(f"Failed to validate flowsheet file '{path}': {e}")
        logger.debug("Validation traceback:", exc_info=True)
        raise SystemExit(1)

    # Added after schema validation so additionalProperties:false doesn't reject it.
    config["_config_path"] = path
    return config


def validate_snapshot_config(state: "SnapshotState", base_name: str) -> None:
    """Warn if a loaded snapshot config has missing required fields."""
    from ..types import FlowsheetConfig

    snap_cfg = state.config
    try:
        FlowsheetConfig.from_dict(snap_cfg)
    except (TypeError, KeyError) as exc:
        logger.warning(
            f"Snapshot '{state.snapshot_id}' for '{base_name}' has an "
            f"incompatible config schema: {exc}. "
            f"Delete 'outputs/{base_name}.pfstate' and re-run to create a "
            f"fresh snapshot."
        )


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

def build_run_metadata(config: dict, solver_tol: float, solver_max_iter: int, backend: str) -> dict:
    """Build run metadata for checkpoint storage."""
    config_bytes = json.dumps(config, sort_keys=True).encode("utf-8")
    flowsheet_hash = hashlib.sha256(config_bytes).hexdigest()[:16]
    return {
        "version": _pf_version,
        "flowsheet_hash": flowsheet_hash,
        "solver_settings": {
            "tol": solver_tol,
            "max_iter": solver_max_iter,
            "backend": backend,
        },
    }


def extract_providers(flowsheet_path: str) -> dict:
    """Read flowsheet JSON and return the raw providers dict."""
    try:
        with open(flowsheet_path, "r", encoding="utf-8") as f:
            config = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.error(f"Failed to read flowsheet '{flowsheet_path}': {e}")
        raise SystemExit(1)
    providers = config.get("providers", {})
    if not providers:
        logger.warning(f"No providers declared in '{flowsheet_path}'.")
    return providers


# ---------------------------------------------------------------------------
# State management
# ---------------------------------------------------------------------------

def load_state_manager(outputs_dir: str, base_name: str) -> tuple[StateManager, "SnapshotState | None"]:
    """Create a StateManager and load the snapshot (if it exists)."""
    state_path = os.path.join(outputs_dir, f"{base_name}.pfstate")
    sm = StateManager(state_path)
    state = sm.load_state()
    return sm, state


# ---------------------------------------------------------------------------
# Provider checking
# ---------------------------------------------------------------------------

def _resolve_provider_url(cfg: dict, ptype: str) -> str:
    """Build the provider URL, falling back to the default port."""
    from ..providers.registry import get_provider_default_port

    url = cfg.get("url")
    if not url:
        port = get_provider_default_port(ptype) or 9000
        url = f"http://localhost:{port}"
    return url


def check_providers(
    config: dict,
    flowsheet_path: str,
    *,
    fail_fast: bool = True,
    verbose: bool = False,
) -> None:
    """Verify all declared providers are reachable or importable.

    Parameters
    ----------
    fail_fast:
        If ``True`` (default), raise ``SystemExit`` on the first failure.
        If ``False``, accumulate errors and raise once at the end.
    verbose:
        If ``True``, log per-provider success messages.
    """
    from ..providers.registry import is_containerized, _PROVIDER_CATALOG

    providers = config.get("providers", {})
    errors: list[str] = []

    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = _resolve_provider_url(cfg, ptype)
            try:
                urllib.request.urlopen(f"{url}/health", timeout=5)
                if verbose:
                    logger.info(f"Provider '{name}' — reachable at {url}")
            except (urllib.error.URLError, OSError, TimeoutError) as exc:
                msg = (
                    f"Provider '{name}' unreachable at {url}. "
                    f"Run: pf init {flowsheet_path}"
                )
                if fail_fast:
                    logger.error(msg)
                    raise SystemExit(1)
                logger.error(msg)
                errors.append(msg)
        else:
            catalog = _PROVIDER_CATALOG.get(ptype, {})
            module = catalog.get("module", "")
            try:
                importlib.util.find_spec(module)
                if verbose:
                    logger.info(f"Provider '{name}' — importable (pip)")
            except (ModuleNotFoundError, ValueError):
                dep = catalog.get("optional_dep")
                hint = f"pip install 'processforge[{dep}]'" if dep else "built-in"
                msg = (
                    f"Provider '{name}' not installed. "
                    f"Run: pf init {flowsheet_path}  (install with: {hint})"
                )
                if fail_fast:
                    logger.error(msg)
                    raise SystemExit(1)
                logger.error(msg)
                errors.append(msg)

    if errors:
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# Divergence reports
# ---------------------------------------------------------------------------

def log_residual_breakdown(fs: Any) -> list[dict]:
    """Log the top residual violators from a failed solve and return the breakdown."""
    breakdown = getattr(fs, "residual_breakdown", [])
    if breakdown:
        logger.error("Top residual violators:")
        for entry in breakdown:
            logger.error(
                f"  [{entry['index']:4d}] {entry['var_name']:<35s}  |F| = {entry['residual']:.4e}"
            )
    else:
        logger.debug("No residual breakdown available from solver.")
    return breakdown


def build_divergence_report(
    *,
    drifted_params: list[str],
    solver_stats: dict,
    x_last: Any,
    var_names: list[str],
    breakdown: list[dict],
) -> dict:
    """Construct the divergence report dict for failed solves."""
    x_last_list = x_last.tolist() if hasattr(x_last, "tolist") else (x_last if isinstance(x_last, list) else [])
    var_names_list = var_names if var_names else []

    if not x_last_list:
        logger.debug("x_last unavailable — divergence report will have empty x_last field.")
    if not var_names_list:
        logger.debug("var_names unavailable — divergence report will have empty var_names field.")

    return {
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "drifted_params": drifted_params,
        "final_norm": solver_stats.get("final_norm"),
        "solver_stats": solver_stats,
        "x_last": x_last_list,
        "var_names": var_names_list,
        "top_violators": breakdown,
    }


def write_divergence_report(outputs_dir: str, base_name: str, divergence: dict) -> None:
    """Write the divergence report JSON and log its path."""
    div_path = os.path.join(outputs_dir, f"{base_name}_divergence.json")
    with open(div_path, "w", encoding="utf-8") as f:
        json.dump(divergence, f, indent=2)
    logger.error(f"Divergence report written to {div_path}")


# ---------------------------------------------------------------------------
# State saving (with error handling)
# ---------------------------------------------------------------------------

def save_snapshot(
    sm: StateManager,
    config: dict,
    x_converged: Any,
    var_names: list[str],
    *,
    metadata: dict,
    parent_snapshot_id: str | None,
    label: str = "snapshot",
) -> str:
    """Save a state snapshot, logging a clear message on failure."""
    try:
        return sm.save_state(
            config, x_converged, var_names,
            metadata=metadata,
            parent_snapshot_id=parent_snapshot_id,
        )
    except Exception as exc:
        logger.error(f"Failed to save {label}: {exc}")
        logger.warning("Results were computed but state was not persisted.")
        raise SystemExit(1)
