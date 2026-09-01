"""``pf export-fmu`` — export flowsheet as FMI 2.0 co-simulation FMU."""

from __future__ import annotations

import os
from typing import Literal

import typer
from loguru import logger

from .common import require_existing_file


def export_fmu(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
    output_dir: str = typer.Option(
        "outputs",
        "--output-dir",
        "-o",
        help="Directory for the output FMU (default: outputs/)",
    ),
    backend: Literal["scipy", "pyomo", "casadi"] = typer.Option(
        "scipy",
        "--backend",
        help="EO solver backend for steady-state mode (default: scipy)",
    ),
    pathsim: bool = typer.Option(
        False,
        "--pathsim",
        help="Also emit the fmu_interface.json manifest and a PathSim wiring example",
    ),
) -> None:
    """Export a flowsheet as an FMI 2.0 co-simulation FMU."""
    require_existing_file(flowsheet)

    from ..fmu import build_fmu  # local import — pythonfmu is optional

    output_dir = output_dir or "outputs"
    backend = backend or "scipy"

    try:
        fmu_path = build_fmu(
            flowsheet,
            output_dir=output_dir,
            backend=backend,
            write_manifest=pathsim,
        )
        logger.info(f"FMU written to: {fmu_path}")
        if pathsim:
            base = os.path.splitext(os.path.basename(fmu_path))[0]
            manifest_path = os.path.join(output_dir, f"{base}_fmu_interface.json")
            example_path = os.path.join(output_dir, f"{base}_pathsim.py")
            _write_pathsim_example(example_path, fmu_path, manifest_path)
            logger.info(f"Manifest written to: {manifest_path}")
            logger.info(f"PathSim example written to: {example_path}")
    except Exception as e:
        logger.error(f"FMU export failed: {type(e).__name__}: {e}")
        logger.debug("FMU export traceback:", exc_info=True)
        raise SystemExit(1)


def _write_pathsim_example(example_path: str, fmu_path: str, manifest_path: str) -> None:
    """Write a minimal PathSim wiring example next to the exported FMU."""
    content = f'''"""Auto-generated PathSim wiring example."""
import pathsim as ps
from processforge.pathsim import ProcessForgeFMU

plant = ProcessForgeFMU(
    fmu_path={fmu_path!r},
    manifest_path={manifest_path!r},
)

# Inspect available ports:
# print(plant.inputs)
# print(plant.outputs)

# Example: drive an input with a constant setpoint and read an output.
setpoint = ps.Constant(val=1.0)
sim = ps.Simulation(
    blocks=[setpoint, plant],
    connections=[
        # ps.Connection(setpoint, plant["<unit>.<input>"]),
    ],
)
sim.run(time_span=(0, 10))
'''
    with open(example_path, "w", encoding="utf-8") as f:
        f.write(content)
