"""High-level, library-safe Python wrapper for running processforge flowsheets.

This module is the primary integration point for web applications (FastAPI,
Django, Celery, …) and programmatic callers.  It exposes :func:`run_flowsheet`
and :func:`apply_flowsheet`, which perform the full simulation lifecycle —
validate, check providers, solve, persist the pfarchive, and upload all
processforge outputs to S3 — and return structured results instead of calling
``SystemExit``.

Provider artifacts produced by containerized backends (e.g. OpenMC h5 files)
continue to be uploaded from the container side when ``S3_BUCKET`` is
configured there.  This wrapper additionally uploads the processforge archive
(manifest, Zarr results, snapshots, registry files) from the host/FastAPI
process so the whole output footprint is durable in S3.
"""
from __future__ import annotations

import importlib.util
import json
import os
import time
import urllib.error
import urllib.request
from typing import Any, Literal

from jsonschema import ValidationError as JsonschemaValidationError
from loguru import logger
from pydantic import BaseModel, Field

from .cli.common import (
    build_divergence_report,
    build_run_metadata,
    display_backend,
    flowsheet_basename,
    load_state_manager,
    log_residual_breakdown,
    output_root,
    write_divergence_report,
)
from .cli.persist import flowsheet_hash, make_run_id, persist_run
from .eo import EOFlowsheet
from .eo.solver import EOSolver, solve_with_homotopy
from .flowsheet import Flowsheet
from .persistence.archive import ProcessStateArchive
from .persistence.artifact_store import ArtifactStore
from .provenance import build_dynamic_x0, build_run_info
from .providers.manager import teardown_providers
from .providers.registry import (
    _PROVIDER_CATALOG,
    get_provider_default_port,
    is_containerized,
)
from .types import EngineOutput, RunManifest
from .utils.s3_upload import validate_s3
from .utils.validate_flowsheet import _validate_config_impl

HEALTH_MAX_ATTEMPTS = 6
HEALTH_RETRY_DELAY = 5  # seconds


# ---------------------------------------------------------------------------
# Typed exceptions
# ---------------------------------------------------------------------------
class ProcessforgeRunError(Exception):
    """Base exception for high-level runner failures."""


class FlowsheetValidationError(ProcessforgeRunError):
    """Raised when a flowsheet fails schema or semantic validation."""


class ProviderUnavailableError(ProcessforgeRunError):
    """Raised when a declared provider is unreachable or not installed."""


class ProviderRunError(ProcessforgeRunError):
    """Raised when a provider unit reports status=failed."""


class ConvergenceError(ProcessforgeRunError):
    """Raised when a steady-state solve (or homotopy fallback) fails to converge."""

    def __init__(self, message: str, divergence_report: dict | None = None):
        super().__init__(message)
        self.divergence_report = divergence_report or {}


class StatePersistenceError(ProcessforgeRunError):
    """Raised when persisting state/manifest/Zarr outputs fails."""


# ---------------------------------------------------------------------------
# Result models
# ---------------------------------------------------------------------------
class RunResult(BaseModel):
    """Structured result returned by :func:`run_flowsheet`."""

    run_id: str
    flowsheet_hash: str
    mode: str
    converged: bool
    status: str  # converged | not_converged | up_to_date
    archive_path: str
    manifest_path: str
    remote_uris: list[str] = Field(default_factory=list)
    artifact_uris: list[str] = Field(
        default_factory=list,
        description="S3 URIs of provider-produced artifacts (populated by the container side).",
    )
    archive_uris: list[str] = Field(
        default_factory=list,
        description="S3 URIs of processforge archive files uploaded from this process.",
    )
    elapsed_s: float
    solver_stats: dict = Field(default_factory=dict)
    snapshot_id: str | None = None
    divergence_report: dict | None = None
    backend: str = ""


class ApplyResult(RunResult):
    """Structured result returned by :func:`apply_flowsheet`."""


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
def _resolve_config(flowsheet: str | dict) -> tuple[dict, str]:
    """Load and validate a flowsheet from a path or dict.

    Returns ``(config, source_name)``.  Raises :class:`FlowsheetValidationError`
    on any validation failure.  The input dict is deep-copied so validation
    mutations (e.g. expanding ``material_mix`` into ``z``) do not leak back to
    the caller.
    """
    if isinstance(flowsheet, dict):
        source_name = flowsheet.get("metadata", {}).get("name", "<dict>")
        raw = json.loads(json.dumps(flowsheet))
        try:
            config = _validate_config_impl(raw, source_name=source_name)
        except (JsonschemaValidationError, ValueError) as exc:
            raise FlowsheetValidationError(str(exc)) from exc
        return config, source_name

    if not os.path.exists(flowsheet):
        raise FlowsheetValidationError(f"Flowsheet file not found: {flowsheet}")

    with open(flowsheet, "r", encoding="utf-8") as f:
        raw = json.load(f)
    try:
        config = _validate_config_impl(raw, source_name=flowsheet)
    except (JsonschemaValidationError, ValueError) as exc:
        raise FlowsheetValidationError(str(exc)) from exc
    # Runtime providers need the original source path.
    config["_config_path"] = flowsheet
    return config, flowsheet


def _resolve_outputs_dir(outputs_dir: str | None) -> str:
    return outputs_dir or output_root()


def _resolve_s3_bucket(bucket: str | None) -> str | None:
    return bucket or os.environ.get("S3_BUCKET")


def _resolve_s3_prefix(prefix: str | None) -> str:
    return prefix or os.environ.get("S3_PREFIX", "processforge")


def _ping_provider_health(url: str, timeout: int = 5) -> tuple[bool, dict | str]:
    try:
        with urllib.request.urlopen(f"{url.rstrip('/')}/health", timeout=timeout) as resp:
            return True, json.loads(resp.read().decode())
    except (urllib.error.URLError, OSError, TimeoutError, ValueError) as exc:
        return False, str(exc)


def _resolve_provider_url(cfg: dict, ptype: str) -> str:
    """Return the provider service URL, deriving a default port if needed."""
    url = cfg.get("url")
    if not url:
        port = get_provider_default_port(ptype) or 9000
        url = f"http://localhost:{port}"
    return url


def _check_providers(config: dict, source_name: str) -> None:
    """Verify all declared providers are reachable or importable.

    Mirrors ``cli.common.check_providers`` but raises
    :class:`ProviderUnavailableError` instead of ``SystemExit``.
    """
    providers = config.get("providers", {})
    errors: list[str] = []

    if providers:
        logger.info("=== Provider / Container Health ===")

    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = _resolve_provider_url(cfg, ptype)
            ok, info = False, "not probed"
            for attempt in range(1, HEALTH_MAX_ATTEMPTS + 1):
                ok, info = _ping_provider_health(url, timeout=5)
                if ok:
                    break
                if attempt < HEALTH_MAX_ATTEMPTS:
                    logger.info(
                        f"  Waiting for '{name}' container to become healthy at {url}… "
                        f"(attempt {attempt}/{HEALTH_MAX_ATTEMPTS})"
                    )
                    time.sleep(HEALTH_RETRY_DELAY)

            if ok:
                payload = info if isinstance(info, dict) else {}
                status = payload.get("status", "?")
                provider_type = payload.get("provider_type", "?")
                logger.info(
                    f"  [OK] {name} [{ptype}] {url} — "
                    f"status={status} provider_type={provider_type}"
                )
            else:
                msg = (
                    f"  [ERR] {name} [{ptype}] {url} — unreachable after "
                    f"{HEALTH_MAX_ATTEMPTS} attempt(s): {info}."
                )
                logger.error(msg)
                errors.append(f"Provider '{name}' unreachable at {url}: {info}")
        else:
            catalog = _PROVIDER_CATALOG.get(ptype)
            module = catalog.module if catalog else ""
            try:
                importlib.util.find_spec(module)
                logger.info(f"  [OK] {name} [{ptype}] (pip — importable)")
            except (ModuleNotFoundError, ValueError):
                dep = catalog.optional_dep if catalog else None
                hint = f"pip install 'processforge[{dep}]'" if dep else "built-in"
                msg = f"  [WARN] {name} [{ptype}] — not installed. (install with: {hint})"
                logger.warning(msg)
                errors.append(f"Provider '{name}' not installed ({hint})")

    if errors:
        raise ProviderUnavailableError("; ".join(errors))


def _check_for_failed_units(fs: Any) -> None:
    """Fail loudly if any SolverUnit returned ``status='failed'``."""
    failed = {
        name: out
        for name, out in getattr(fs, "engine_outputs", {}).items()
        if isinstance(out, EngineOutput) and out.status == "failed"
    }
    if not failed:
        return

    messages: list[str] = []
    for name, out in failed.items():
        err = out.error
        category = getattr(err, "category", "unknown") if err else "unknown"
        message = getattr(err, "message", "") if err else ""
        hint = getattr(err, "hint", "") if err else ""
        msg = f"Unit '{name}' simulation FAILED [{category}]: {message}"
        if hint:
            msg += f" (hint: {hint})"
        messages.append(msg)
        logger.error(msg)
    raise ProviderRunError("; ".join(messages))


def _log_active_providers(config: dict) -> None:
    """Log a one-line summary of active providers for visibility."""
    providers = config.get("providers", {})
    if not providers:
        return

    containers: list[str] = []
    pips: list[str] = []
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = _resolve_provider_url(cfg, ptype)
            containers.append(f"{name} [{ptype}] @ {url}")
        else:
            pips.append(f"{name} [{ptype}]")

    if containers:
        logger.info("=== Active providers (containerized) ===")
        for entry in containers:
            logger.info(f"  → {entry}")
    if pips:
        logger.info("=== Active providers (local) ===")
        for entry in pips:
            logger.info(f"  → {entry}")


def _artifact_uris_from_manifest(manifest: RunManifest | None) -> list[str]:
    """Collect S3 URIs already populated on provider artifacts."""
    if manifest is None:
        return []
    uris: list[str] = []
    for out in manifest.units.values():
        for art in out.artifacts:
            uris.extend(art.remote_uris)
    for out in manifest.streams.values():
        for art in out.artifacts:
            uris.extend(art.remote_uris)
    return uris


def _upload_pfarchive(
    archive_path: str,
    run_id: str,
    flowsheet_hash: str,
    bucket: str | None,
    prefix: str,
) -> list[str]:
    """Upload the entire pfarchive to S3 (no-op when bucket is unset)."""
    if not bucket:
        return []
    store = ArtifactStore(bucket=bucket, prefix=prefix)
    uris = store.persist_archive(
        archive_path, run_id=run_id, flowsheet_hash=flowsheet_hash
    )
    logger.info(f"Uploaded pfarchive to {len(uris)} S3 object(s)")
    return uris


def _save_snapshot_safe(
    archive: ProcessStateArchive,
    config: dict,
    x_converged: Any,
    var_names: list[str],
    *,
    metadata: dict,
    parent_snapshot_id: str | None = None,
) -> str:
    """Save a snapshot and convert failures to :class:`StatePersistenceError`."""
    try:
        return archive.save_snapshot(
            config,
            x_converged,
            var_names,
            metadata=metadata,
            parent_snapshot_id=parent_snapshot_id,
        )
    except Exception as exc:
        raise StatePersistenceError(f"Failed to save snapshot: {exc}") from exc


def _build_run_result(
    *,
    fs: Any,
    config: dict,
    archive: ProcessStateArchive,
    run_id: str,
    flowsheet_hash: str,
    mode: str,
    converged: bool,
    status: str,
    elapsed_s: float,
    snapshot_id: str | None,
    s3_bucket: str | None,
    s3_prefix: str,
) -> RunResult:
    """Assemble the final :class:`RunResult` and upload the archive to S3."""
    manifest = archive.load_run(run_id)
    artifact_uris = _artifact_uris_from_manifest(manifest)
    archive_uris = _upload_pfarchive(
        archive.path, run_id, flowsheet_hash, s3_bucket, s3_prefix
    )
    backend_label = display_backend(config, getattr(fs, "backend", "dynamic"))
    logger.info(f"Backend      : {backend_label}")
    logger.info(f"Run saved    : {os.path.join(archive.path, 'runs', run_id + '.json')}")
    return RunResult(
        run_id=run_id,
        flowsheet_hash=flowsheet_hash,
        mode=mode,
        converged=converged,
        status=status,
        archive_path=archive.path,
        manifest_path=os.path.join(archive.path, "runs", f"{run_id}.json"),
        remote_uris=sorted(set(artifact_uris + archive_uris)),
        artifact_uris=sorted(set(artifact_uris)),
        archive_uris=sorted(set(archive_uris)),
        elapsed_s=elapsed_s,
        solver_stats=getattr(fs, "solver_stats", {}) or {},
        snapshot_id=snapshot_id,
        backend=backend_label,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def run_flowsheet(
    flowsheet: str | dict,
    *,
    outputs_dir: str | None = None,
    backend: Literal["scipy", "pyomo", "casadi"] | None = None,
    export_images: bool = False,
    s3_bucket: str | None = None,
    s3_prefix: str | None = None,
) -> RunResult:
    """Run a process simulation and persist all outputs (locally + S3).

    Args:
        flowsheet: Path to a flowsheet JSON file or an already-loaded dict.
        outputs_dir: Directory for local pfarchive output. Defaults to
            ``PROCESSFORGE_OUTPUT_DIR`` or ``outputs/``.
        backend: Optional EO solver backend override (``scipy``/``pyomo``/``casadi``).
        export_images: If ``True``, generate PNG plots of the results.
        s3_bucket: S3 bucket for archive upload. Defaults to ``S3_BUCKET`` env var.
        s3_prefix: S3 key prefix. Defaults to ``S3_PREFIX`` env var or
            ``"processforge"``.

    Returns:
        :class:`RunResult` with run metadata, local paths, and S3 URIs.

    Raises:
        FlowsheetValidationError: On invalid flowsheet input.
        ProviderUnavailableError: On unreachable/missing providers.
        ProviderRunError: When a provider unit fails.
        StatePersistenceError: When local persistence fails.
    """
    s3_bucket = _resolve_s3_bucket(s3_bucket)
    s3_prefix = _resolve_s3_prefix(s3_prefix)
    if s3_bucket:
        validate_s3()

    config, source_name = _resolve_config(flowsheet)
    _check_providers(config, source_name)
    _log_active_providers(config)

    base_name = (
        flowsheet_basename(source_name)
        if isinstance(flowsheet, str)
        else config.get("metadata", {}).get("name", "flowsheet")
    )
    outputs_dir = _resolve_outputs_dir(outputs_dir)
    os.makedirs(outputs_dir, exist_ok=True)

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    is_dynamic = mode == "dynamic"

    run_id = make_run_id()
    fhash = flowsheet_hash(config)

    t0 = time.perf_counter()

    if is_dynamic:
        archive_path = os.path.join(outputs_dir, f"{base_name}.pfarchive")
        archive = ProcessStateArchive(archive_path)
        state = archive.load_snapshot()
        if state is not None:
            stream_inits = archive.state_to_stream_dicts(state)
            streams_cfg = config.get("streams", {})
            for s_name, s_vals in stream_inits.items():
                if s_name in streams_cfg:
                    streams_cfg[s_name].update(s_vals)
                else:
                    logger.debug(
                        f"Stream '{s_name}' from snapshot not in current flowsheet — skipped."
                    )
            logger.info("Using .pfarchive converged state as dynamic t=0.")
        else:
            logger.debug(
                "No .pfarchive found — starting dynamic run from flowsheet defaults."
            )

        fs = Flowsheet(config)
        fs.set_run_context(run_id, fhash)
        logger.info("=== Dynamic Results ===")
        results = fs.run()

        _check_for_failed_units(fs)

        if hasattr(fs, "converged"):
            if fs.converged:
                logger.info("Dynamic simulation converged.")
            else:
                logger.warning(
                    "Dynamic simulation did NOT converge. Results may be unreliable."
                )

        try:
            x0, var_names = build_dynamic_x0(config)
        except Exception as exc:
            raise StatePersistenceError(
                f"Failed to build initial state vector: {exc}"
            ) from exc

        run_info = build_run_info(config, x0=x0, var_names=var_names)
        snap_x, snap_vn, snap_backend, snap_ok = x0, var_names, "dynamic", True
    else:
        fs = EOFlowsheet(config, backend=backend)
        fs.set_run_context(run_id, fhash)
        logger.info("=== Steady-State EO Results ===")
        results = fs.run()

        _check_for_failed_units(fs)

        if hasattr(fs, "converged"):
            if fs.converged:
                logger.info("Steady-state simulation converged.")
            else:
                logger.warning(
                    "Steady-state simulation did NOT converge. Results may be unreliable."
                )

        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
        snap_x, snap_vn = fs.x_converged, fs.var_names
        snap_backend = getattr(fs, "backend", "scipy")
        snap_ok = bool(getattr(fs, "converged", False))

    archive_path = os.path.join(outputs_dir, f"{base_name}.pfarchive")
    archive = ProcessStateArchive(archive_path)
    persist_run(archive, fs, run_id, results, run_info, config, base_name)

    snapshot_id: str | None = None
    if (is_dynamic or snap_ok) and snap_x is not None:
        meta = build_run_metadata(config, 1e-6, 50, snap_backend)
        snapshot_id = _save_snapshot_safe(
            archive,
            config,
            snap_x,
            snap_vn,
            metadata=meta,
            parent_snapshot_id=None,
        )

    result = _build_run_result(
        fs=fs,
        config=config,
        archive=archive,
        run_id=run_id,
        flowsheet_hash=fhash,
        mode=mode,
        converged=bool(getattr(fs, "converged", is_dynamic)),
        status="converged" if (is_dynamic or fs.converged) else "not_converged",
        elapsed_s=time.perf_counter() - t0,
        snapshot_id=snapshot_id,
        s3_bucket=s3_bucket,
        s3_prefix=s3_prefix,
    )

    if export_images:
        try:
            from .result import plot_results, plot_timeseries

            plot_results(results, fname=f"{base_name}_results.png")
            plot_timeseries(results, fname=f"{base_name}_timeseries.png")
            logger.info(
                f"Plots saved: {base_name}_results.png, {base_name}_timeseries.png"
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"Failed to generate plots: {type(exc).__name__}: {exc}")

    return result


def apply_flowsheet(
    flowsheet: str | dict,
    *,
    outputs_dir: str | None = None,
    backend: Literal["scipy", "pyomo", "casadi"] | None = None,
    tolerance: float = 1e-6,
    max_iter: int = 50,
    skip_homotopy: bool = False,
    s3_bucket: str | None = None,
    s3_prefix: str | None = None,
) -> ApplyResult:
    """Apply a steady-state flowsheet with warm-start / homotopy fallback.

    Args mirror :func:`run_flowsheet` with solver-specific controls.

    Returns:
        :class:`ApplyResult`.  If no drift is detected against the existing
        snapshot, returns ``status='up_to_date'`` without running the solver.

    Raises:
        ConvergenceError: When both direct solve and homotopy fallback fail.
    """
    s3_bucket = _resolve_s3_bucket(s3_bucket)
    s3_prefix = _resolve_s3_prefix(s3_prefix)
    if s3_bucket:
        validate_s3()

    config, source_name = _resolve_config(flowsheet)
    _check_providers(config, source_name)

    base_name = (
        flowsheet_basename(source_name)
        if isinstance(flowsheet, str)
        else config.get("metadata", {}).get("name", "flowsheet")
    )
    outputs_dir = _resolve_outputs_dir(outputs_dir)
    os.makedirs(outputs_dir, exist_ok=True)

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    if mode != "steady":
        raise FlowsheetValidationError(
            "apply_flowsheet is only supported for steady-state EO flowsheets."
        )

    archive, state = load_state_manager(outputs_dir, base_name)

    # Structural diff: detect topology changes.
    topology_changed = False
    if state is not None:
        diff = archive.detect_structural_diff(config, state)
        topology_changed = bool(diff.get("added") or diff.get("removed"))
        if topology_changed:
            logger.warning(
                "Topology changed (units added/removed). "
                "Homotopy requires identical topology — falling back to cold start."
            )

    # Parameter drift (only meaningful when topology is unchanged).
    drifted: list[str] = []
    current_metadata = build_run_metadata(config, tolerance, max_iter, backend or "scipy")
    if state is not None and not topology_changed:
        mismatches = archive.validate_metadata(current_metadata, state)
        if mismatches:
            logger.warning(f"Metadata mismatch: {mismatches}")
        drifted = archive.detect_drift(config, state)
        if not drifted:
            logger.info("No drift detected. System is already at the desired state.")
            return ApplyResult(
                run_id="",
                flowsheet_hash=flowsheet_hash(config),
                mode="steady",
                converged=True,
                status="up_to_date",
                archive_path=archive.path,
                manifest_path="",
                elapsed_s=0.0,
                snapshot_id=state.snapshot_id,
                backend=display_backend(config, backend or "scipy"),
            )
        logger.warning("Drift detected:")
        old_config = state.config if hasattr(state, "config") else state.get("config", {})
        for param in drifted:
            old_val = old_config
            new_val = config
            for part in param.split("."):
                old_val = old_val.get(part, "") if isinstance(old_val, dict) else ""
                new_val = new_val.get(part, "") if isinstance(new_val, dict) else ""
            logger.warning(f"  {param}: {old_val!r} → {new_val!r}")

    fs = EOFlowsheet(config, backend=backend)
    fs.saved_state = state if not topology_changed else None
    fs.solver_tol = tolerance
    fs.solver_max_iter = max_iter

    logger.info("=== Running Apply (Steady-State EO) ===")
    run_id = make_run_id()
    fhash = flowsheet_hash(config)
    t0 = time.perf_counter()
    fs.set_run_context(run_id, fhash)
    results = fs.run()
    elapsed = time.perf_counter() - t0

    logger.info(
        f"Direct solve completed in {elapsed:.2f}s (converged={fs.converged})."
    )

    if fs.converged:
        snapshot_id = _save_snapshot_safe(
            archive,
            config,
            fs.x_converged,
            fs.var_names,
            metadata=current_metadata,
            parent_snapshot_id=state.snapshot_id if state is not None and not topology_changed else None,
        )
        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
        persist_run(archive, fs, run_id, results, run_info, config, base_name, snapshot_id)
        logger.info("=== Apply Summary ===")
        logger.info("  Status       : CONVERGED")
        logger.info(
            f"  Final ||F||  : {fs.solver_stats.get('final_norm', '?'):.3e}"
        )
        logger.info(f"  Iterations   : {fs.solver_stats.get('iterations', '?')}")
        logger.info(f"  Backend      : {display_backend(config, fs.backend)}")
        logger.info(f"  Snapshot ID  : {snapshot_id}")
        logger.info(f"  Run ID       : {run_id}")
        logger.info(f"  Elapsed (s)  : {elapsed:.2f}")
        return _build_run_result(
            fs=fs,
            config=config,
            archive=archive,
            run_id=run_id,
            flowsheet_hash=fhash,
            mode="steady",
            converged=True,
            status="converged",
            elapsed_s=elapsed,
            snapshot_id=snapshot_id,
            s3_bucket=s3_bucket,
            s3_prefix=s3_prefix,
        )

    # Direct solve failed — try homotopy (only when topology is same and state exists).
    if state is not None and not topology_changed and drifted and not skip_homotopy:
        logger.warning("Direct solve failed. Attempting homotopy continuation...")
        solver = EOSolver(backend=fs.backend, tol=tolerance, max_iter=max_iter)
        run_id = make_run_id()
        tmp_fs = EOFlowsheet(config, backend=backend)
        tmp_fs.set_run_context(run_id, fhash)
        manager = tmp_fs._build()
        try:
            homotopy_result = solve_with_homotopy(tmp_fs, manager, solver, state, drifted)
            if homotopy_result.converged:
                tmp_fs.assemble_from_solution(
                    manager, homotopy_result.x_solution, True, homotopy_result.stats
                )
        finally:
            teardown_providers(tmp_fs._provider_map)

        if homotopy_result.converged:
            logger.info(
                f"Homotopy converged: ||F||="
                f"{homotopy_result.stats.get('final_norm', '?'):.3e}, "
                f"iterations={homotopy_result.stats.get('iterations', '?')}"
            )
            snapshot_id = _save_snapshot_safe(
                archive,
                config,
                homotopy_result.x_solution,
                fs.var_names,
                metadata=current_metadata,
                parent_snapshot_id=state.snapshot_id if state is not None and not topology_changed else None,
            )
            run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
            persist_run(archive, tmp_fs, run_id, tmp_fs.results, run_info, config, base_name, snapshot_id)
            logger.info("Homotopy apply succeeded. New snapshot saved.")
            return _build_run_result(
                fs=tmp_fs,
                config=config,
                archive=archive,
                run_id=run_id,
                flowsheet_hash=fhash,
                mode="steady",
                converged=True,
                status="converged",
                elapsed_s=time.perf_counter() - t0,
                snapshot_id=snapshot_id,
                s3_bucket=s3_bucket,
                s3_prefix=s3_prefix,
            )

        logger.error("Homotopy also failed to converge.")
        prev_id = state.snapshot_id if state is not None else "unknown"
        archive.rollback(1)
        logger.warning(f"Reverted .pfstate to snapshot before {prev_id}.")
        breakdown = log_residual_breakdown(fs)
        divergence = build_divergence_report(
            drifted_params=drifted,
            solver_stats=homotopy_result.stats,
            x_last=homotopy_result.x_solution,
            var_names=fs.var_names,
            breakdown=breakdown,
        )
    else:
        if not drifted and state is None:
            logger.error("Cold-start solve failed to converge (no prior snapshot).")
        elif skip_homotopy:
            logger.error("Cold-start solve failed to converge (--skip-homotopy).")
        else:
            logger.error("Cold-start solve failed to converge.")
        breakdown = log_residual_breakdown(fs)
        divergence = build_divergence_report(
            drifted_params=[],
            solver_stats=fs.solver_stats,
            x_last=getattr(fs, "x_converged", []),
            var_names=getattr(fs, "var_names", []),
            breakdown=breakdown,
        )

    write_divergence_report(outputs_dir, base_name, divergence)
    raise ConvergenceError(
        f"Flowsheet '{base_name}' failed to converge.",
        divergence_report=divergence,
    )
