import argparse
import json
import os

from loguru import logger

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

def _cmd_apply(args):
    """Apply a flowsheet change using state drift detection and warm-starting."""
    import os
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    try:
        from .utils.validate_flowsheet import validate_flowsheet
        config = validate_flowsheet(fname)
    except SystemExit:
        raise
    except Exception as e:
        logger.error(f"Failed to validate flowsheet file '{fname}': {e}")
        raise SystemExit(1)

    config["_config_path"] = fname
    base_name = os.path.splitext(os.path.basename(fname))[0]
    
    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    if mode != "steady":
        logger.error("pf apply is only supported for steady-state EO flows.")
        raise SystemExit(1)

    state_path = os.path.join("outputs", f"{base_name}.pfstate")
    manager = StateManager(state_path)
    
    state = manager.load_state()
    if state is not None:
        drifted = manager.detect_drift(config, state)
        if not drifted:
            logger.info("✅ No drift detected. System is already at the desired state.")
            return
        logger.info(f"⚠️ Drift detected in parameters: {drifted}")
    else:
        logger.info("No prior state found. Running cold start...")
        drifted = []

    from .eo import EOFlowsheet
    from .provenance import build_run_info
    from .result import save_results_zarr
    
    fs = EOFlowsheet(config, backend=None)
    fs.state_manager = manager
    fs.saved_state = state
    fs.drifted_params = drifted
    
    logger.info("=== Running Apply (Steady-State EO) ===")
    results = fs.run()
    
    if hasattr(fs, "x_converged"):
        manager.save_state(config, fs.x_converged, fs.var_names)
        
    run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
    save_results_zarr(
        results,
        os.path.join("outputs", f"{base_name}_results.zarr"),
        run_info=run_info,
    )

def _cmd_run(args):
    """Run a process simulation from a flowsheet JSON file."""
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    try:
        config = validate_flowsheet(fname)
    except SystemExit:
        raise
    except Exception as e:
        logger.error(f"Failed to validate flowsheet file '{fname}': {e}")
        raise SystemExit(1)

    # Stash the config file path for providers that need it (e.g. ModelicaProvider).
    # Added after schema validation so additionalProperties:false doesn't reject it.
    config["_config_path"] = fname

    base_name = os.path.splitext(os.path.basename(fname))[0]

    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    is_dynamic = mode == "dynamic"

    if is_dynamic:
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

    zarr_path = save_results_zarr(
        results,
        os.path.join("outputs", f"{base_name}_results.zarr"),
        run_info=run_info,
    )
    if args.export_images:
        plot_results(results, fname=f"{base_name}_results.png")
        plot_timeseries(results, fname=f"{base_name}_timeseries.png")



def _cmd_export_fmu(args):
    """Export a flowsheet as an FMI 2.0 co-simulation FMU."""
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

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
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

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
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    with open(fname, "r") as f:
        flowsheet_schema = json.load(f)

    output_dir = args.output_dir or "."
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
    from .analysis.dof import SystemDOFReport
    logger.info("=== DOF Analysis ===")
    if report.component_names:
        logger.info(f"Components: {', '.join(report.component_names)}  (N_c = {report.n_components})")
    else:
        logger.warning("No components found in feed streams — DOF analysis may be incomplete.")

    for r in report.per_unit:
        icon = "✅" if r.status == "determined" else ("⚠️" if r.status == "under-specified" else "❌")
        logger.info(
            f"  {r.unit_name:<16} [{r.unit_type:<18}] "
            f"unknowns={r.n_unknowns}  equations={r.n_equations}  "
            f"DOF={r.dof}  {icon} {r.status}"
        )
        for issue in r.issues:
            logger.warning(f"    → {r.unit_name}: {issue}")

    logger.info(f"Feed stream specs:  {report.feed_stream_specs}")
    logger.info(f"Total unknowns:     {report.total_unknowns}")
    logger.info(f"Total equations:    {report.total_equations} (unit) + {report.feed_stream_specs} (feed) = {report.total_equations + report.feed_stream_specs}")

    if report.system_dof == 0:
        logger.info("System DOF: 0  ✅ Exactly determined — ready to solve")
    elif report.system_dof > 0:
        logger.warning(f"System DOF: {report.system_dof}  ⚠️  Under-specified — add {report.system_dof} more spec(s)")
    else:
        logger.error(f"System DOF: {report.system_dof}  ❌ Over-specified — remove {abs(report.system_dof)} spec(s)")


def _print_unit_mismatches(mismatches) -> None:
    for m in mismatches:
        if not m.compatible:
            logger.error(f"❌ Unit mismatch — stream '{m.stream_name}'.{m.property_name}: {m.message}")
        else:
            logger.warning(f"⚠️  Unit annotation — stream '{m.stream_name}'.{m.property_name}: {m.message}")


def _cmd_plan(args):
    """Parse a PCL or JSON flowsheet, run validators, DOF analysis, and emit a Mermaid diagram."""
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"File '{fname}' not found.")
        raise SystemExit(1)

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

    # Step 6: Mermaid diagram
    if not args.no_diagram:
        from .utils.mermaid_diagram import generate_mermaid
        output_dir = args.output_dir or "."
        base_name = os.path.splitext(os.path.basename(fname))[0]
        out_path = os.path.join(output_dir, f"{base_name}_plan.mmd")
        generate_mermaid(config, output_path=out_path)
        logger.info(f"Mermaid diagram → {out_path}")

    # Exit non-zero on hard errors
    hard_errors = [m for m in mismatches if not m.compatible]
    if hard_errors or report.system_dof < 0:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(
        description="ProcessForge - Chemical Process Simulation",
        prog="processforge",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # processforge apply
    apply_parser = subparsers.add_parser("apply", help="Apply flowsheet using state-based warm start and homotopy")
    apply_parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")
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
        help="Validate a flowsheet, run DOF analysis, and generate a Mermaid diagram",
    )
    plan_parser.add_argument("flowsheet", help="Path to .pcl or .json flowsheet file")
    plan_parser.add_argument(
        "--output-dir", "-o", default=".",
        help="Output directory for the Mermaid diagram (default: current directory)",
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
        default=".",
        help="Output directory (default: current directory)",
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
        "run": _cmd_run,
        "apply": _cmd_apply,
        "plan": _cmd_plan,
        "diagram": _cmd_diagram,
        "export-fmu": _cmd_export_fmu,
        "export-modelica": _cmd_export_modelica,
    }
    commands[args.command](args)


if __name__ == "__main__":
    main()

def _cmd_apply(args):
    """Apply the closest existing state, detect drift, block or warm-start, fallback to homotopy."""
    from .state import StateManager
    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Flowsheet file '{fname}' not found.")
        raise SystemExit(1)

    try:
        config = validate_flowsheet(fname)
    except SystemExit:
        raise
    
    config["_config_path"] = fname
    base_name = os.path.splitext(os.path.basename(fname))[0]
    
    sim_cfg = config.get("simulation", {})
    mode = sim_cfg.get("mode", "steady")
    if mode != "steady":
        logger.error("pf apply is only supported for steady-state EO flows.")
        raise SystemExit(1)
        
    state_path = os.path.join("outputs", f"{base_name}.pfstate")
    manager = StateManager(state_path)
    
    # Drift Detection
    state = manager.load_state()
    if state is not None:
        drifted = manager.detect_drift(config, state)
        if not drifted:
            logger.info("✅ No drift detected. System is already at the desired state.")
            return
        else:
            logger.info(f"⚠️ Drift detected in parameters: {drifted}")
    else:
        logger.info("No prior state found. Running cold start...")
        drifted = []

    # Run standard flowsheet with possible warm-start
    fs = EOFlowsheet(config, backend=None)
    fs.state_manager = manager
    fs.saved_state = state
    
    logger.info("=== Running Apply... ===")
    results = fs.run()
    
    # Ensure save state
    if hasattr(fs, "x_converged"):
        manager.save_state(config, fs.x_converged, fs.var_names)
        
    # Also save outputs.zarr
    run_info = build_run_info(config, x0=fs.x0, var_names=fs.var_names)
    save_results_zarr(
        results,
        os.path.join("outputs", f"{base_name}_results.zarr"),
        run_info=run_info,
    )
