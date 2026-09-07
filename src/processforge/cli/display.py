"""Display helpers for CLI output — DOF reports, unit mismatches, structural diffs."""

from __future__ import annotations

from typing import Any

import importlib.util

from loguru import logger


def print_dof_report(report: Any) -> None:
    """Pretty-print a Degrees of Freedom analysis report."""
    logger.info("=== Degrees of Freedom Analysis ===")
    if report.component_names:
        logger.info(f"Components: {', '.join(report.component_names)}  (N_c = {report.n_components})")
    else:
        logger.warning("No components found in feed streams — Degrees of Freedom analysis may be incomplete.")

    for r in report.per_unit:
        icon = "OK" if r.status == "determined" else ("WARN" if r.status == "under-specified" else "ERR")
        logger.info(
            f"  {r.unit_name:<8} {f'[{r.unit_type}]':<8} "
            f"variables={r.n_variables}  equations={r.n_equations}  "
            f"Degrees of Freedom={r.dof}  [{icon}] {r.status}"
        )
        for issue in r.issues:
            logger.warning(f"    -> {r.unit_name}: {issue}")

    logger.info(f"Feed stream specs:  {report.feed_stream_specs}")
    logger.info(f"Total variables:    {report.total_variables}")
    logger.info(
        f"Total equations:    {report.total_equations} (unit) "
        f"+ {report.feed_stream_specs} (feed) = "
        f"{report.total_equations + report.feed_stream_specs}"
    )

    if report.system_dof == 0:
        logger.info("System Degrees of Freedom: 0  Exactly determined — ready to solve")
    elif report.system_dof > 0:
        logger.warning(f"System Degrees of Freedom: {report.system_dof}  Under-specified — add {report.system_dof} more spec(s)")
    else:
        logger.warning(f"System Degrees of Freedom: {report.system_dof}  Over-specified — remove {-report.system_dof} spec(s)")


def print_unit_mismatches(mismatches: list[Any]) -> None:
    """Log each unit-annotation mismatch from a Pint consistency check."""
    if not mismatches:
        logger.debug("No unit mismatches found.")
        return
    for m in mismatches:
        if not m.compatible:
            logger.error(f"Unit mismatch — stream '{m.stream_name}'.{m.property_name}: {m.message}")
        else:
            logger.warning(f"Unit annotation — stream '{m.stream_name}'.{m.property_name}: {m.message}")


def print_provider_health(config: dict, strict: bool = False) -> list[str]:
    """Log per-provider reachability for flowsheets using Docker containers.

    For each provider declared in *config*:

    * containerized providers are probed via ``GET /health`` and reported with
      their URL, ``status``, and ``provider_type`` from the health payload;
    * pip-installable providers are reported as importable or not.

    Returns a list of failure messages — one per unreachable containerized
    provider (empty means every provider is healthy/available). When *strict*
    is True, a missing pip-installable provider is also appended to failures.
    """
    from ..providers.registry import is_containerized, _PROVIDER_CATALOG

    providers = config.get("providers", {})
    logger.info("=== Provider / Container Health ===")
    failures: list[str] = []
    if not providers:
        logger.info("  No providers declared.")
        return failures

    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            from .common import _resolve_provider_url, _ping_provider_health

            url = _resolve_provider_url(cfg, ptype)
            ok, info = _ping_provider_health(url, timeout=5)
            if ok:
                payload = info if isinstance(info, dict) else {}
                status = payload.get("status", "?")
                provider_type = payload.get("provider_type", "?")
                logger.info(
                    f"  [OK] {name} [{ptype}] {url} — status={status} provider_type={provider_type}"
                )
            else:
                msg = f"  [ERR] {name} [{ptype}] {url} — unreachable: {info}"
                logger.error(msg)
                failures.append(f"Provider '{name}' unreachable at {url}: {info}")
        else:
            catalog = _PROVIDER_CATALOG.get(ptype)
            module = catalog.module if catalog else ""
            try:
                importlib.util.find_spec(module)
                logger.info(f"  [OK] {name} [{ptype}] (pip — importable)")
            except (ModuleNotFoundError, ValueError):
                dep = catalog.optional_dep if catalog else None
                hint = f"pip install 'processforge[{dep}]'" if dep else "built-in"
                logger.warning(f"  [WARN] {name} [{ptype}] — not installed ({hint})")
                if strict:
                    failures.append(
                        f"Provider '{name}' not installed ({hint})"
                    )

    return failures


def print_structural_diff(diff: dict) -> None:
    """Print a +/~/- structural diff of units."""
    logger.info("=== Structural Diff vs. Saved State ===")
    for name, unit_type in diff.get("added", {}).items():
        logger.info(f"  + {name:<20} [{unit_type}]  (added)")
    for name, info in diff.get("modified", {}).items():
        unit_type = info.get("type", "?")
        changes = info.get("changes", [])
        changes_str = ", ".join(changes)
        logger.info(f"  ~ {name:<20} [{unit_type}]  {changes_str}")
    for name, unit_type in diff.get("removed", {}).items():
        logger.info(f"  - {name:<20} [{unit_type}]  (removed)")
    if not any(diff.get(k) for k in ("added", "modified", "removed")):
        logger.info("  (no structural changes)")


def _resolve_path(cfg: dict, path: str) -> object:
    """Resolve a dotted config path (e.g. ``simulation.tf``) to its value."""
    cur: object = cfg
    for part in path.split("."):
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            return "<missing>"
    return cur


def print_param_drift(drifted: list[str], old_config: dict, new_config: dict) -> None:
    """Log parameter drift paths with their old → new values."""
    logger.info("=== Parameter Drift vs. Saved State ===")
    if not drifted:
        logger.info("  (no parameter changes)")
        return
    for path in drifted:
        old_val = _resolve_path(old_config, path)
        new_val = _resolve_path(new_config, path)
        logger.info(f"  ~ {path}: {old_val!r} → {new_val!r}")
