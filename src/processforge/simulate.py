import argparse
import json
import os
import sys
import time
import hashlib
from loguru import logger

from . import __version__ as _pf_version
from .utils.validate_flowsheet import validate_flowsheet, validate_flowsheet_dict
from .flowsheet import Flowsheet
from .eo import EOFlowsheet
from .provenance import build_dynamic_x0, build_run_info
from .result import (
    plot_results,
    plot_timeseries,
    save_results_zarr,
)
from .utils.flowsheet_diagram import draw_flowsheet
from .state import StateManager

from .utils.mermaid_diagram import generate_mermaid

def _output_root() -> str:
    """Root directory for run outputs (zarr, provider artifacts).

    Defaults to ``outputs`` for local runs; the Docker image sets
    ``PROCESSFORGE_OUTPUT_DIR=/data`` so outputs land on the mounted volume.
    """
    return os.environ.get("PROCESSFORGE_OUTPUT_DIR", "outputs")


def _require_existing_file(path: str, label: str = "Flowsheet file") -> None:
    """Fail fast with a consistent message when an input path is missing."""
    if not os.path.exists(path):
        logger.error(f"{label} '{path}' not found.")
        raise SystemExit(1)


def _validate_runtime_flowsheet(path: str) -> dict:
    """Validate flowsheet config and attach source path for runtime providers."""
    try:
        config = validate_flowsheet(path)
    except Exception as e:
        logger.error(f"Failed to validate flowsheet file '{path}': {e}")
        raise SystemExit(1)

    # Added after schema validation so additionalProperties:false doesn't reject it.
    config["_config_path"] = path
    return config


def _build_run_metadata(config: dict, solver_tol: float, solver_max_iter: int, backend: str) -> dict:
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


def _extract_providers(flowsheet_path: str) -> dict:
    """Read flowsheet JSON and return the raw providers dict."""
    with open(flowsheet_path, "r", encoding="utf-8") as f:
        config = json.load(f)
    return config.get("providers", {})


def _cmd_init(args):
    """Initialise the .processforge/ project directory."""
    from .providers.registry import (
        is_containerized,
        get_provider_docker_image,
        get_provider_default_port,
        _PROVIDER_CATALOG,
    )
    from .lock import write_lock
    from .compose import generate_compose

    root = args.path or "."
    pf_dir = os.path.join(root, ".processforge")
    outputs_dir = os.path.join(root, "outputs")

    os.makedirs(pf_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)

    # Write config.json (always)
    config_path = os.path.join(pf_dir, "config.json")
    if not os.path.exists(config_path):
        default_config = {
            "version": 1,
            "default_backend": "scipy",
            "outputs_dir": "outputs",
        }
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(default_config, f, indent=2)
        logger.info(f"Created {config_path}")
    else:
        logger.info(f"{config_path} already exists — skipped.")

    # No flowsheet → scaffold only
    if not args.flowsheet:
        logger.info(".processforge/ initialised successfully.")
        logger.info("To set up providers: pf init <flowsheet.json>")
        return

    # Read providers from flowsheet
    flowsheet_path = args.flowsheet
    if not os.path.exists(flowsheet_path):
        logger.error(f"Flowsheet '{flowsheet_path}' not found.")
        raise SystemExit(1)

    providers = _extract_providers(flowsheet_path)
    logger.info(f"Reading providers from {flowsheet_path}...")

    # Categorize providers
    pip_providers = {}
    docker_providers = {}
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = cfg.get("url")
            if not url:
                port = get_provider_default_port(ptype) or 9000
                url = f"http://localhost:{port}"
            docker_providers[name] = {
                "type": ptype,
                "url": url,
                "docker_image": cfg.get("docker_image") or get_provider_docker_image(ptype),
                "port": get_provider_default_port(ptype),
            }
            logger.info(f"  {name}: type={ptype}, url={url} (Docker)")
        else:
            pip_providers[name] = {"type": ptype}
            logger.info(f"  {name}: type={ptype} (pip)")

    # Validate pip providers are importable
    import importlib

    for name, info in pip_providers.items():
        ptype = info["type"]
        catalog = _PROVIDER_CATALOG.get(ptype, {})
        module = catalog.get("module", "")
        try:
            importlib.util.find_spec(module)
            logger.info(f"  ✓ {name} — importable")
        except (ModuleNotFoundError, ValueError):
            dep = catalog.get("optional_dep")
            hint = f"pip install 'processforge[{dep}]'" if dep else "built-in"
            logger.warning(f"  ✗ {name} — not installed. Install with: {hint}")

    # Generate compose for Docker providers
    if docker_providers:
        compose_path = os.path.join(pf_dir, "docker-compose.yml")
        if os.path.exists(compose_path):
            logger.warning(
                f"⚠ Environment already initialized — reinitializing from {flowsheet_path}"
            )

        generate_compose(pf_dir, docker_providers, outputs_dir)
        logger.info(f"Generated {compose_path}")

        # Attempt docker compose pull
        import subprocess

        try:
            result = subprocess.run(
                ["docker", "compose", "-f", compose_path, "pull"],
                capture_output=True,
                text=True,
                timeout=120,
            )
            if result.returncode == 0:
                logger.info("Pulled Docker images.")
            else:
                logger.warning(f"docker compose pull failed: {result.stderr}")
        except FileNotFoundError:
            logger.warning(
                "Docker not found. Install Docker to use containerized providers."
            )
        except subprocess.TimeoutExpired:
            logger.warning("docker compose pull timed out.")
    else:
        logger.info("No containerized providers — skipping Docker setup.")

    # Write lock file
    lock_providers = {}
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            lock_providers[name] = {
                "docker_image": cfg.get("docker_image") or get_provider_docker_image(ptype),
                "url": cfg.get("url")
                or f"http://localhost:{get_provider_default_port(ptype) or 9000}",
            }
        else:
            lock_providers[name] = {
                "docker_image": None,
                "url": None,
            }

    from . import __version__ as pf_version

    write_lock(pf_dir, flowsheet_path, lock_providers, pf_version)
    logger.info(f"Wrote {os.path.join(pf_dir, 'lock.json')}")
    logger.info(".processforge/ initialised successfully.")


def _cmd_validate(args):
    """Check that providers and environment are ready for a flowsheet."""
    import importlib
    from .providers.registry import is_containerized, _PROVIDER_CATALOG

    fname = args.flowsheet
    _require_existing_file(fname, label="Flowsheet file")

    config = _validate_runtime_flowsheet(fname)
    providers = config.get("providers", {})

    all_ok = True

    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        catalog = _PROVIDER_CATALOG.get(ptype, {})

        if is_containerized(ptype):
            # Check Docker service reachability
            url = cfg.get("url")
            if not url:
                from .providers.registry import get_provider_default_port

                port = get_provider_default_port(ptype) or 9000
                url = f"http://localhost:{port}"

            import urllib.request
            import urllib.error

            try:
                urllib.request.urlopen(f"{url}/health", timeout=5)
                logger.info(f"✓ Provider '{name}' — reachable at {url}")
            except (urllib.error.URLError, OSError, TimeoutError):
                logger.error(f"✗ Provider '{name}' — unreachable at {url}")
                logger.info(
                    "  Start with: docker compose -f .processforge/docker-compose.yml up -d"
                )
                all_ok = False
        else:
            # Check pip importability
            module = catalog.get("module", "")
            try:
                importlib.util.find_spec(module)
                logger.info(f"✓ Provider '{name}' — importable (pip)")
            except (ModuleNotFoundError, ValueError):
                dep = catalog.get("optional_dep")
                hint = f"pip install 'processforge[{dep}]'" if dep else "built-in"
                logger.error(f"✗ Provider '{name}' — not installed. Install with: {hint}")
                all_ok = False

    if not all_ok:
        raise SystemExit(1)

    logger.info("All providers ready.")


def _check_providers(config: dict, flowsheet_path: str) -> None:
    """Verify all declared providers are reachable or importable. Exits on failure."""
    import importlib
    import urllib.request
    import urllib.error
    from .providers.registry import is_containerized, get_provider_default_port, _PROVIDER_CATALOG

    providers = config.get("providers", {})
    for name, cfg in providers.items():
        ptype = cfg.get("type", "")
        if is_containerized(ptype):
            url = cfg.get("url")
            if not url:
                port = get_provider_default_port(ptype) or 9000
                url = f"http://localhost:{port}"
            try:
                urllib.request.urlopen(f"{url}/health", timeout=5)
            except (urllib.error.URLError, OSError, TimeoutError):
                logger.error(
                    f"Provider '{name}' unreachable at {url}. "
                    f"Run: pf init {flowsheet_path}"
                )
                raise SystemExit(1)
        else:
            catalog = _PROVIDER_CATALOG.get(ptype, {})
            module = catalog.get("module", "")
            try:
                importlib.util.find_spec(module)
            except (ModuleNotFoundError, ValueError):
                logger.error(
                    f"Provider '{name}' not installed. "
                    f"Run: pf init {flowsheet_path}"
                )
                raise SystemExit(1)


def _cmd_run(args):
    """Run a process simulation from a flowsheet JSON file."""
    fname = args.flowsheet
    _require_existing_file(fname)
    config = _validate_runtime_flowsheet(fname)

    # Check provider availability
    _check_providers(config, fname)

    base_name = os.path.splitext(os.path.basename(fname))[0]
    outputs_dir = _output_root()

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    is_dynamic = mode == "dynamic"

    if is_dynamic:
        # Load .pfstate as t=0 if available
        state_path = os.path.join(outputs_dir, f"{base_name}.pfstate")
        sm = StateManager(state_path)
        state = sm.load_state()
        if state is not None:
            stream_inits = sm.state_to_stream_dicts(state)
            for s_name, s_vals in stream_inits.items():
                if s_name in config.get("streams", {}):
                    config["streams"][s_name].update(s_vals)
            logger.info("Using .pfstate converged state as dynamic t=0.")

        fs = Flowsheet(config)
        logger.info("=== Dynamic Results ===")
        results = fs.run()
        x0, var_names = build_dynamic_x0(config)
        run_info = build_run_info(config, x0=x0, var_names=var_names)
    else:
        # Pass None so EOFlowsheet resolves backend from config (with scipy default).
        fs = EOFlowsheet(config, backend=None)
        logger.info("=== Steady-State EO Results ===")
        results = fs.run()
        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)

    os.makedirs(outputs_dir, exist_ok=True)
    save_results_zarr(
        results,
        os.path.join(outputs_dir, f"{base_name}_results.zarr"),
        run_info=run_info,
    )
    if args.export_images:
        plot_results(results, fname=f"{base_name}_results.png")
        plot_timeseries(results, fname=f"{base_name}_timeseries.png")


def _cmd_export_fmu(args):
    """Export a flowsheet as an FMI 2.0 co-simulation FMU."""
    fname = args.flowsheet
    _require_existing_file(fname)

    from .fmu import build_fmu  # local import — pythonfmu is optional

    output_dir = args.output_dir or "outputs"
    backend = args.backend or "scipy"

    try:
        fmu_path = build_fmu(fname, output_dir=output_dir, backend=backend)
        logger.info(f"FMU written to: {fmu_path}")
    except Exception as e:
        logger.error(f"FMU export failed: {e}")
        raise SystemExit(1)


def _cmd_export_modelica(args):
    """Transpile a flowsheet to Modelica and optionally compile via OMPython."""
    fname = args.flowsheet
    _require_existing_file(fname)

    from .modelica import transpile, compile_modelica
    from .modelica.transpiler import _derive_model_name
    from .utils.validate_flowsheet import validate_flowsheet as _vf

    output_dir = args.output_dir or "outputs"

    try:
        mo_path = transpile(fname, output_dir=output_dir)
    except Exception as e:
        logger.error(f"Modelica transpilation failed: {e}")
        raise SystemExit(1)

    if not args.no_compile:
        try:
            config = _vf(fname)
            model_name = _derive_model_name(config, fname)
            fmu_path = compile_modelica(mo_path, model_name, output_dir=output_dir)
            logger.info(f"Model Exchange FMU written to: {fmu_path}")
        except Exception as e:
            logger.error(f"OMPython compilation failed: {e}")
            raise SystemExit(1)


def _cmd_diagram(args):
    """Generate a flowsheet diagram from a JSON file."""
    fname = args.flowsheet
    _require_existing_file(fname)

    with open(fname, "r", encoding="utf-8") as f:
        flowsheet_schema = json.load(f)

    output_dir = args.output_dir or "diagrams"
    fmt = args.format or "png"

    base_name = os.path.splitext(os.path.basename(fname))[0]

    output_path = draw_flowsheet(
        flowsheet_schema,
        output_directory=output_dir,
        output_filename=base_name,
        file_format=fmt,
    )
    logger.info(f"Diagram saved to {output_path}")


def _print_dof_report(report) -> None:
    logger.info("=== Degrees of Freedom Analysis ===")
    if report.component_names:
        logger.info(f"Components: {', '.join(report.component_names)}  (N_c = {report.n_components})")
    else:
        logger.warning("No components found in feed streams — Degrees of Freedom analysis may be incomplete.")

    for r in report.per_unit:
        icon = "✅" if r.status == "determined" else ("⚠️" if r.status == "under-specified" else "❌")
        logger.info(
            f"  {r.unit_name:<8} {f'[{r.unit_type}]':<8} "
            f"variables={r.n_variables}  equations={r.n_equations}  "
            f"Degrees of Freedom={r.dof}  {icon} {r.status}"
        )
        for issue in r.issues:
            logger.warning(f"    → {r.unit_name}: {issue}")

    logger.info(f"Feed stream specs:  {report.feed_stream_specs}")
    logger.info(f"Total variables:    {report.total_variables}")
    logger.info(f"Total equations:    {report.total_equations} (unit) + {report.feed_stream_specs} (feed) = {report.total_equations + report.feed_stream_specs}")

    if report.system_dof == 0:
        logger.info("System Degrees of Freedom: 0  ✅ Exactly determined — ready to solve")
    elif report.system_dof > 0:
        logger.warning(f"System Degrees of Freedom: {report.system_dof}  ⚠️  Under-specified — add {report.system_dof} more spec(s)")
    else:
        logger.warning(f"System Degrees of Freedom: {report.system_dof}  ❌ Over-specified — remove {-report.system_dof} spec(s)")


def _print_unit_mismatches(mismatches) -> None:
    for m in mismatches:
        if not m.compatible:
            logger.error(f"❌ Unit mismatch — stream '{m.stream_name}'.{m.property_name}: {m.message}")
        else:
            logger.warning(f"⚠️  Unit annotation — stream '{m.stream_name}'.{m.property_name}: {m.message}")


def _print_structural_diff(diff: dict) -> None:
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


def _cmd_plan(args):
    """Parse a PCL or JSON flowsheet, run validators, DOF analysis, structural diff, and emit a Mermaid diagram."""
    fname = args.flowsheet
    _require_existing_file(fname, label="File")

    # Step 1: Load config
    if fname.endswith(".pcl"):
        from .pcl import load_pcl, PCLCompileError
        try:
            json_config = load_pcl(fname)
        except PCLCompileError as e:
            logger.error(f"PCL compile error: {e}")
            raise SystemExit(1)
    else:
        try:
            with open(fname, "r", encoding="utf-8") as f:
                json_config = json.load(f)
        except json.JSONDecodeError as e:
            logger.error(f"JSON parse error in '{fname}': {e}")
            raise SystemExit(1)

    # Step 2: Pint unit consistency check (before _units are stripped)
    from .utils.unit_consistency import check_unit_consistency, strip_units_annotations
    mismatches = check_unit_consistency(json_config)
    _print_unit_mismatches(mismatches)

    # Step 3: Strip _units annotations before schema validation
    strip_units_annotations(json_config)

    # Step 4: Schema + connectivity validation
    try:
        config = validate_flowsheet_dict(json_config, source_name=fname)
    except SystemExit:
        raise

    # Step 5: DOF analysis
    from .analysis.dof import analyze_dof
    report = analyze_dof(config)
    _print_dof_report(report)

    # Step 5b: DOF fix suggestions
    if report.system_dof != 0:
        logger.warning("=== DOF Fix Suggestions ===")
        for r in report.per_unit:
            if r.dof > 0:
                for issue in r.issues:
                    logger.warning(f"  Unit '{r.unit_name}' [{r.unit_type}]: {issue}")
            elif r.dof < 0:
                logger.warning(
                    f"  Unit '{r.unit_name}' [{r.unit_type}]: "
                    f"possibly over-specified by {abs(r.dof)} equation(s)"
                )

    # Step 6: Structural diff vs. saved state
    base_name = os.path.splitext(os.path.basename(fname))[0]
    outputs_dir = _output_root()
    state_path = os.path.join(outputs_dir, f"{base_name}.pfstate")
    sm = StateManager(state_path)
    state = sm.load_state()
    diff = None
    if state is not None:
        diff = sm.detect_structural_diff(config, state)
        _print_structural_diff(diff)
    else:
        logger.info("=== Structural Diff vs. Saved State ===")
        logger.info("  No prior state found — this will be a cold start.")

    # Step 6b: Warm-start and homotopy eligibility
    logger.info("=== Warm-Start Status ===")
    if state is not None:        
        snap_id = state.snapshot_id
        snap_ts = state.timestamp
        logger.info(f"  Warm-start available : Yes  (snapshot {snap_id}, {snap_ts})")
        topology_ok = diff is None or not bool(diff.get("added") or diff.get("removed"))
        logger.info(
            f"  Homotopy eligible    : {'Yes' if topology_ok else 'No (topology changed)'}"
        )
    else:
        logger.info("  Warm-start available : No snapshot found")
        logger.info("  Homotopy eligible    : No (cold start)")

    # Step 7: Mermaid diagram
    if not args.no_diagram:
        output_dir = args.output_dir or "diagrams"
        out_path = os.path.join(output_dir, f"{base_name}_plan.mmd")
        generate_mermaid(config, output_path=out_path)
        logger.info(f"Mermaid diagram → {out_path}")

    # Exit non-zero on hard errors
    hard_errors = [m for m in mismatches if not m.compatible]
    if hard_errors or report.system_dof < 0:
        raise SystemExit(1)


def _cmd_apply(args):
    """Apply flowsheet: drift detection, warm-start, homotopy fallback, convergence guardrails."""
    fname = args.flowsheet
    _require_existing_file(fname)
    config = _validate_runtime_flowsheet(fname)

    # Check provider availability
    _check_providers(config, fname)

    base_name = os.path.splitext(os.path.basename(fname))[0]

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    if mode != "steady":
        logger.error("pf apply is only supported for steady-state EO flowsheets.")
        raise SystemExit(1)

    outputs_dir = _output_root()
    os.makedirs(outputs_dir, exist_ok=True)
    state_path = os.path.join(outputs_dir, f"{base_name}.pfstate")
    sm = StateManager(state_path)

    state = sm.load_state()

    # Structural diff: detect topology changes
    topology_changed = False
    if state is not None:
        diff = sm.detect_structural_diff(config, state)
        _print_structural_diff(diff)
        topology_changed = bool(diff.get("added") or diff.get("removed"))
        if topology_changed:
            logger.warning(
                "Topology changed (units added/removed). "
                "Homotopy requires identical topology — falling back to cold start."
            )

    # Parameter drift (only meaningful when topology is unchanged)
    drifted: list[str] = []
    current_metadata = _build_run_metadata(config, args.tolerance, args.max_iter, args.backend)
    if state is not None and not topology_changed:
        mismatches = sm.validate_metadata(current_metadata, state)
        if mismatches:
            logger.warning(f"⚠️ Metadata mismatch: {mismatches}")
        drifted = sm.detect_drift(config, state)
        if not drifted:
            logger.info("✅ No drift detected. System is already at the desired state.")
            return
        logger.info(f"⚠️ Drift detected: {drifted}")

    # Build flowsheet; attach saved state for warm-start unless topology changed
    fs = EOFlowsheet(config, backend=args.backend)
    fs.saved_state = state if not topology_changed else None
    fs.solver_tol = args.tolerance
    fs.solver_max_iter = args.max_iter

    logger.info("=== Running Apply (Steady-State EO) ===")
    t0 = time.perf_counter()
    results = fs.run()
    elapsed = time.perf_counter() - t0

    if fs.converged:
        snapshot_id = sm.save_state(
            config, fs.x_converged, fs.var_names,
            metadata=current_metadata,
            parent_snapshot_id=state.snapshot_id if state is not None and not topology_changed else None,
        )
        run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
        zarr_path = os.path.join(outputs_dir, f"{base_name}_results.zarr")
        save_results_zarr(results, zarr_path, run_info=run_info)
        logger.info("=== Apply Summary ===")
        logger.info("  Status       : CONVERGED")
        logger.info(f"  Final ||F||  : {fs.solver_stats.get('final_norm', '?'):.3e}")
        logger.info(f"  Iterations   : {fs.solver_stats.get('iterations', '?')}")
        logger.info(f"  Backend      : {fs.backend}")
        logger.info(f"  Snapshot ID  : {snapshot_id}")
        logger.info(f"  Results zarr : {zarr_path}")
        logger.info(f"  Elapsed (s)  : {elapsed:.2f}")
        return

    # Direct solve failed — try homotopy (only when topology is same and state exists)
    if state is not None and not topology_changed and drifted and not args.skip_homotopy:
        logger.warning("Direct solve failed. Attempting homotopy continuation...")
        from .eo.solver import EOSolver, solve_with_homotopy

        solver = EOSolver(backend=fs.backend, tol=args.tolerance, max_iter=args.max_iter)
        from .eo.flowsheet import EOFlowsheet as _EO
        tmp_fs = _EO(config, backend=args.backend)
        from processforge.providers.manager import teardown_providers
        manager = tmp_fs._build()
        try:
            x_hom, hom_converged, hom_stats = solve_with_homotopy(
                tmp_fs, manager, solver, state, drifted
            )
        finally:
            teardown_providers(tmp_fs._provider_map)

        if hom_converged:
            sm.save_state(
                config, x_hom, fs.var_names,
                metadata=current_metadata,
                parent_snapshot_id=state.snapshot_id if state is not None and not topology_changed else None,
            )
            run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
            save_results_zarr(
                results,
                os.path.join(outputs_dir, f"{base_name}_results.zarr"),
                run_info=run_info,
            )
            logger.info("✅ Homotopy apply succeeded. New snapshot saved.")
            return

        # Both failed — auto-revert and write divergence report
        logger.error("Homotopy also failed to converge.")
        sm.rollback(1)
        logger.warning("Reverted .pfstate to last good snapshot.")
        breakdown = getattr(fs, "residual_breakdown", [])
        if breakdown:
            logger.error("Top residual violators:")
            for entry in breakdown:
                logger.error(
                    f"  [{entry['index']:4d}] {entry['var_name']:<35s}  |F| = {entry['residual']:.4e}"
                )
        import datetime
        divergence = {
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "drifted_params": drifted,
            "final_norm": hom_stats.get("final_norm"),
            "solver_stats": hom_stats,
            "x_last": x_hom.tolist(),
            "var_names": fs.var_names,
            "top_violators": breakdown,
        }
    else:
        # Cold start also failed; no rollback (nothing to revert to)
        logger.error("Cold-start solve failed to converge.")
        breakdown = getattr(fs, "residual_breakdown", [])
        if breakdown:
            logger.error("Top residual violators:")
            for entry in breakdown:
                logger.error(
                    f"  [{entry['index']:4d}] {entry['var_name']:<35s}  |F| = {entry['residual']:.4e}"
                )
        import datetime
        divergence = {
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "drifted_params": [],
            "final_norm": fs.solver_stats.get("final_norm"),
            "solver_stats": fs.solver_stats,
            "x_last": fs.x_converged.tolist() if hasattr(fs, "x_converged") else [],
            "var_names": fs.var_names if hasattr(fs, "var_names") else [],
            "top_violators": breakdown,
        }

    div_path = os.path.join(outputs_dir, f"{base_name}_divergence.json")
    with open(div_path, "w", encoding="utf-8") as f:
        json.dump(divergence, f, indent=2)
    logger.error(f"Divergence report written to {div_path}")
    raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Processforge - Chemical Process Simulation",
        prog="processforge",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Enable detailed debug tracebacks for errors"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # processforge init
    init_parser = subparsers.add_parser(
        "init",
        help="Initialise .processforge/ project directory and provider environment",
    )
    init_parser.add_argument(
        "flowsheet", nargs="?", default=None,
        help="Flowsheet JSON to initialise environment for (omit for scaffold only)",
    )
    init_parser.add_argument(
        "--path", default=".",
        help="Root directory to initialise in (default: current directory)",
    )

    # processforge validate
    validate_parser = subparsers.add_parser(
        "validate",
        help="Check providers and environment are ready for a flowsheet",
    )
    validate_parser.add_argument(
        "flowsheet", help="Path to the flowsheet JSON file",
    )

    # processforge apply
    apply_parser = subparsers.add_parser(
        "apply",
        help="Apply flowsheet using state-based warm start and homotopy fallback",
    )
    apply_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    apply_parser.add_argument(
        "--backend", choices=["scipy", "pyomo", "casadi"], default=None,
        help="Override the flowsheet's simulation.backend",
    )
    apply_parser.add_argument(
        "--tolerance", type=float, default=1e-6,
        help="Newton solver convergence tolerance (default: 1e-6)",
    )
    apply_parser.add_argument(
        "--max-iter", type=int, default=50,
        help="Max Newton iterations (default: 50)",
    )
    apply_parser.add_argument(
        "--skip-homotopy", action="store_true",
        help="Disable homotopy fallback; cold-start only",
    )

    # processforge run
    run_parser = subparsers.add_parser("run", help="Run a process simulation")
    run_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    run_parser.add_argument(
        "--export-images",
        action="store_true",
        help="Generate PNG plots for simulation outputs",
    )

    # processforge plan
    plan_parser = subparsers.add_parser(
        "plan",
        help="Validate a flowsheet, run DOF analysis, structural diff, and generate a Mermaid diagram",
    )
    plan_parser.add_argument("flowsheet", help="Path to .pcl or .json flowsheet file")
    plan_parser.add_argument(
        "--output-dir", "-o", default="diagrams",
        help="Output directory for the Mermaid diagram (default: diagrams)",
    )
    plan_parser.add_argument(
        "--no-diagram", action="store_true",
        help="Skip Mermaid diagram generation",
    )

    # processforge diagram
    diagram_parser = subparsers.add_parser(
        "diagram", help="Generate a flowsheet diagram"
    )
    diagram_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    diagram_parser.add_argument(
        "--output-dir",
        "-o",
        default="diagrams",
        help="Output directory (default: diagrams)",
    )
    diagram_parser.add_argument(
        "--format",
        "-f",
        default="png",
        choices=["png", "svg", "pdf"],
        help="Output format (default: png)",
    )

    # processforge export-modelica
    mo_parser = subparsers.add_parser(
        "export-modelica",
        help="Transpile flowsheet to Modelica .mo and compile via OMPython",
    )
    mo_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    mo_parser.add_argument(
        "--output-dir",
        "-o",
        default="outputs",
        help="Directory for the .mo and .fmu outputs (default: outputs/)",
    )
    mo_parser.add_argument(
        "--no-compile",
        action="store_true",
        help="Write the .mo file only; skip OMPython/omc compilation",
    )

    # processforge export-fmu
    fmu_parser = subparsers.add_parser(
        "export-fmu",
        help="Export flowsheet as FMI 2.0 co-simulation FMU",
    )
    fmu_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
    fmu_parser.add_argument(
        "--output-dir",
        "-o",
        default="outputs",
        help="Directory for the output FMU (default: outputs/)",
    )
    fmu_parser.add_argument(
        "--backend",
        choices=["scipy", "pyomo", "casadi"],
        default="scipy",
        help="EO solver backend for steady-state mode (default: scipy)",
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        raise SystemExit(1)

    commands = {
        "init": _cmd_init,
        "validate": _cmd_validate,
        "run": _cmd_run,
        "apply": _cmd_apply,
        "plan": _cmd_plan,
        "diagram": _cmd_diagram,
        "export-fmu": _cmd_export_fmu,
        "export-modelica": _cmd_export_modelica,
    }
    
    try:
        commands[args.command](args)
    except Exception as e:
        if getattr(args, "debug", False):
            logger.exception("An unhandled exception occurred during command execution:")
            raise SystemExit(1)
        else:
            logger.error(f"{type(e).__name__}: {str(e)}")
            logger.info("Run with --debug for the full traceback.")
            raise SystemExit(1)


if __name__ == "__main__":
    main()
