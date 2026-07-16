"""CLI subpackage — one module per ``pf`` subcommand."""

from __future__ import annotations

import typer


def register_commands(app: typer.Typer) -> None:
    """Register all subcommands on *app*."""
    from .init import init
    from .validate import validate
    from .run import run
    from .apply import apply
    from .plan import plan
    from .diagram import diagram
    from .export_fmu import export_fmu
    from .export_modelica import export_modelica

    app.command(
        "init",
        help="Initialise .processforge/ project directory and provider environment",
    )(init)
    app.command(
        "validate",
        help="Check providers and environment are ready for a flowsheet",
    )(validate)
    app.command("run", help="Run a process simulation")(run)
    app.command(
        "apply",
        help="Apply flowsheet using state-based warm start and homotopy fallback",
    )(apply)
    app.command(
        "plan",
        help="Validate a flowsheet, run DOF analysis, structural diff, and generate a Mermaid diagram",
    )(plan)
    app.command("diagram", help="Generate a flowsheet diagram")(diagram)
    app.command(
        "export-modelica",
        help="Transpile flowsheet to Modelica .mo and compile via OMPython",
    )(export_modelica)
    app.command(
        "export-fmu",
        help="Export flowsheet as FMI 2.0 co-simulation FMU",
    )(export_fmu)
