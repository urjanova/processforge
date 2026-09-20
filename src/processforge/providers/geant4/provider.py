"""Geant4 provider for Monte Carlo particle transport simulations."""

from __future__ import annotations

import threading
from typing import TYPE_CHECKING, Optional

from loguru import logger

from processforge.providers.base_simulation_provider import BaseSimulationProvider
from processforge.providers.errors import (
    ProviderNotAvailableError,
    classify_run_error,
    make_failed_output,
)
from processforge.providers.registry import register_provider
from processforge.types import EngineOutput, OutputProvenance

from .build_helpers import Geant4BuildHelpers
from .result_extraction import extract_results
from .strategies import get_registered_sim_types

if TYPE_CHECKING:
    from processforge.types import (
        FlowsheetConfig,
        Geant4ProviderConfig,
        MaterialDef,
        UnitConfig,
    )


_GEANT4_RUN_LOCK = threading.Lock()


class Geant4Provider(BaseSimulationProvider):
    """Monte Carlo particle transport provider backed by Geant4.

    Responsibilities
    ----------------
    * **Material registry** — ``initialize()`` loads all material defs from the
      flowsheet and indexes them by name.
    * **Working-directory management** — each run executes in
      ``<provider_output_dir>`` via ``G4RunManager`` — no process-wide ``chdir``.
      Falls back to a temp dir when the configured output dir is not writable.
    * **Physics list** — solver config provides the physics list name (default FTFP_BERT).
    * **Simulation dispatch** — ``run_simulation()`` looks up the strategy from
      the shared Geant4 strategy registry by ``sim_type``, calls
      ``strategy.build()``, runs via ``G4RunManager``, then parses the results.
      Failures are returned as ``EngineOutput(status="failed")`` with ``run_dir``
      and error metadata.
    * **Result extraction** — extracts per-shell energy deposition and transmission.
    """

    def __init__(self):
        super().__init__()
        self._provider_output_dir: str = "outputs/geant4"

    def initialize(
        self,
        provider_config: "Geant4ProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        """Verify Geant4 is installed, store config, build material registry."""
        try:
            import geant4  # noqa: F401
        except ImportError as exc:
            raise ProviderNotAvailableError(
                "Geant4 is not installed. Install with: pip install geant4-pybind"
            ) from exc

        self._provider_output_dir = self._expand_output_dir(provider_config.output_dir)

        self._materials = dict(flowsheet_config.materials.items())

        self._initialized = True
        n_mats = len(self._materials)
        logger.info(
            f"Geant4Provider initialized with {n_mats} material(s). "
            f"Registered sim_types: {sorted(get_registered_sim_types())}"
        )

    def teardown(self) -> None:
        """Release provider state."""
        super().teardown()

    @classmethod
    def validate_material(cls, mat_name: str, mat_def, unit_cfg) -> list:
        """Validate Geant4-specific material properties."""
        errors = []

        if mat_def.density is None and mat_def.density_units != "sum":
            errors.append(
                f"Material '{mat_name}' is missing 'density' (required for Geant4; "
                "only 'sum' density_units may omit it)."
            )
        if mat_def.density_units is None:
            errors.append(
                f"Material '{mat_name}' is missing 'density_units' (required for Geant4)."
            )

        return errors

    def run_simulation(self, unit_config: "UnitConfig", inlet: dict) -> "EngineOutput":
        """Build and run a Geant4 simulation from a typed ``UnitConfig``."""
        import geant4

        sim_type = unit_config.sim_type
        strategy_cls = get_registered_sim_types().get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"Geant4Provider: unknown sim_type '{sim_type}'. "
                f"Register with register_geant4_sim_type(). "
                f"Built-in types: {sorted(get_registered_sim_types())}"
            )

        from processforge.schemas.geant4.geant4_model import Geant4Setting

        solver_cfg = Geant4Setting.model_validate(unit_config.solver_config or {})

        geometry_cfg = None
        cfg_model = getattr(strategy_cls, "config_model", None)
        if cfg_model is not None:
            geometry_cfg = cfg_model.model_validate(
                unit_config.geometry_config or {},
                context={"materials": set(self._materials.keys())},
            )

        helpers = Geant4BuildHelpers()
        materials_map: dict = {}
        for key, mdef in self._materials.items():
            g4_mat = helpers.build_material(geant4, key, mdef)
            materials_map[key] = g4_mat

        logger.info(
            f"Geant4Provider: running sim_type='{sim_type}' "
            f"with strategy {strategy_cls.__name__}, "
            f"{solver_cfg.events} events"
        )

        run_dir = self._resolve_run_dir()

        with _GEANT4_RUN_LOCK:
            try:
                sim_components = strategy_cls().build(
                    geant4, solver_cfg, geometry_cfg, materials_map, helpers
                )
                helpers.run_simulation(geant4, sim_components, solver_cfg)
                logger.info(f"Geant4Provider: '{sim_type}' completed in '{run_dir}'")
            except Exception as exc:
                logger.exception(
                    f"Geant4Provider: '{sim_type}' failed in '{run_dir}': {exc}"
                )
                err = classify_run_error("geant4", exc)
                return make_failed_output("geant4", sim_type, run_dir, err)

        fields, artifacts, diagnostics = extract_results(
            geant4, solver_cfg, sim_components, run_dir
        )

        return EngineOutput(
            status="completed",
            engine="geant4",
            sim_type=sim_type,
            fields=fields,
            artifacts=artifacts,
            diagnostics=diagnostics,
            provenance=OutputProvenance(),
            run_dir=str(run_dir),
        )


register_provider("geant4", Geant4Provider)
