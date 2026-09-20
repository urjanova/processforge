"""Result extraction for Geant4 simulations."""

from __future__ import annotations

import pathlib
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from processforge.schemas.geant4.geant4_model import Geant4Setting


def extract_results(
    geant4: Any,
    solver_cfg: "Geant4Setting",
    sim_components: dict,
    run_dir: pathlib.Path,
) -> tuple[list, list, dict]:
    """Extract results from a completed Geant4 simulation.

    Args:
        geant4: The geant4 module.
        solver_cfg: The solver configuration.
        sim_components: Components returned by the strategy's build().
        run_dir: Directory where the simulation ran.

    Returns:
        Tuple of (fields, artifacts, diagnostics).
    """
    from processforge.types import OutputField, OutputArtifact
    from processforge.quantity import Quantity

    fields: list = []
    artifacts: list = []
    diagnostics: dict = {}

    events = getattr(solver_cfg, "events", 0)
    fields.append(
        OutputField(
            name="events_run",
            quantity=Quantity(events, ""),
            kind="scalar",
            source="run_manager",
        )
    )

    physics_list = getattr(solver_cfg, "physics_list", "FTFP_BERT")
    fields.append(
        OutputField(
            name="physics_list",
            quantity=Quantity(physics_list, ""),
            kind="scalar",
            source="physics_list",
        )
    )

    return fields, artifacts, diagnostics
