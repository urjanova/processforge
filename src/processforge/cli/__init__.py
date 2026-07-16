"""CLI subpackage — one module per ``pf`` subcommand."""


def register_commands(subparsers) -> None:
    """Register all subcommands on *subparsers* (an ``argparse._SubParsersAction``)."""
    from .init import cmd_init, add_init_args
    from .validate import cmd_validate, add_validate_args
    from .run import cmd_run, add_run_args
    from .apply import cmd_apply, add_apply_args
    from .plan import cmd_plan, add_plan_args
    from .diagram import cmd_diagram, add_diagram_args
    from .export_fmu import cmd_export_fmu, add_export_fmu_args
    from .export_modelica import cmd_export_modelica, add_export_modelica_args

    _COMMANDS = [
        ("init", add_init_args, cmd_init,
         "Initialise .processforge/ project directory and provider environment"),
        ("validate", add_validate_args, cmd_validate,
         "Check providers and environment are ready for a flowsheet"),
        ("run", add_run_args, cmd_run,
         "Run a process simulation"),
        ("apply", add_apply_args, cmd_apply,
         "Apply flowsheet using state-based warm start and homotopy fallback"),
        ("plan", add_plan_args, cmd_plan,
         "Validate a flowsheet, run DOF analysis, structural diff, and generate a Mermaid diagram"),
        ("diagram", add_diagram_args, cmd_diagram,
         "Generate a flowsheet diagram"),
        ("export-modelica", add_export_modelica_args, cmd_export_modelica,
         "Transpile flowsheet to Modelica .mo and compile via OMPython"),
        ("export-fmu", add_export_fmu_args, cmd_export_fmu,
         "Export flowsheet as FMI 2.0 co-simulation FMU"),
    ]

    for name, add_args, func, help_text in _COMMANDS:
        p = subparsers.add_parser(name, help=help_text)
        add_args(p)
        p.set_defaults(func=func)
