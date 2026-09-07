"""FESTIM provider — in-process hydrogen transport FEM solver."""
from __future__ import annotations

import os
from typing import TYPE_CHECKING

from loguru import logger

from processforge.providers.base_simulation_provider import BaseSimulationProvider
from processforge.providers.errors import ProviderNotAvailableError, classify_run_error, make_failed_output
from processforge.providers.registry import register_provider
from processforge.schemas.festim.festim_model import FestimModel
from processforge.types import EngineOutput, OutputProvenance

from .build_helpers import FestimBuildHelpers
from .result_extraction import extract_results
from .strategies import get_registered_sim_types

if TYPE_CHECKING:
    from processforge.types import FestimProviderConfig, FlowsheetConfig, MaterialDef, UnitConfig


class FestimProvider(BaseSimulationProvider):
    """In-process hydrogen transport provider backed by FESTIM."""

    def __init__(self):
        super().__init__()
        self._provider_output_dir: str = "outputs/festim"

    def initialize(
        self,
        provider_config: "FestimProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        """Verify FESTIM is installed, store config, build material registry."""
        try:
            import festim  # noqa: F401
        except ImportError as exc:
            raise ProviderNotAvailableError(
                "FESTIM is not installed. Install with: pip install festim"
            ) from exc

        self._provider_output_dir = self._expand_output_dir(provider_config.output_dir)
        self._materials = dict(flowsheet_config.materials.items())

        self._initialized = True
        logger.info(
            f"FestimProvider initialized with {len(self._materials)} material(s). "
            f"Registered sim_types: {sorted(get_registered_sim_types())}"
        )

    @classmethod
    def validate_material(cls, mat_name: str, mat_def: "MaterialDef", unit_cfg) -> list:
        """Validate FESTIM-specific material properties."""
        errors = []
        extra = mat_def.extra or {}

        d0 = extra.get("D_0")
        if d0 is None:
            errors.append(
                f"Material '{mat_name}' is missing 'extra.D_0' "
                "(diffusion pre-exponential, required for FESTIM)."
            )
        elif not isinstance(d0, (int, float)):
            errors.append(
                f"Material '{mat_name}' 'extra.D_0' must be numeric, "
                f"got {type(d0).__name__}."
            )
        elif d0 <= 0:
            errors.append(
                f"Material '{mat_name}' 'extra.D_0' must be positive "
                f"(diffusion pre-exponential, m2 s-1), got {d0}."
            )

        ed = extra.get("E_D")
        if ed is None:
            errors.append(
                f"Material '{mat_name}' is missing 'extra.E_D' "
                "(diffusion activation energy, required for FESTIM)."
            )
        elif not isinstance(ed, (int, float)):
            errors.append(
                f"Material '{mat_name}' 'extra.E_D' must be numeric, "
                f"got {type(ed).__name__}."
            )
        elif ed < 0:
            errors.append(
                f"Material '{mat_name}' 'extra.E_D' must be non-negative "
                f"(activation energy, eV), got {ed}."
            )

        has_k = "K_S_0" in extra
        has_e = "E_K_S" in extra
        if has_k != has_e:
            errors.append(
                f"Material '{mat_name}' must define 'extra.K_S_0' and "
                "'extra.E_K_S' together."
            )

        return errors

    def run_simulation(
        self, unit_config: "UnitConfig", inlet: dict
    ) -> "EngineOutput":
        """Build and run a FESTIM simulation from a typed ``UnitConfig``."""
        if not self._initialized:
            raise RuntimeError(
                "FestimProvider has not been initialized — call initialize() first."
            )

        sim_type = unit_config.sim_type
        strategy_cls = get_registered_sim_types().get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"FestimProvider: unknown sim_type '{sim_type}'. "
                f"Register with register_festim_sim_type(). "
                f"Built-in types: {sorted(get_registered_sim_types())}"
            )

        import festim

        festim_model = FestimModel.model_validate(unit_config.solver_config or {})

        helpers = FestimBuildHelpers()
        materials_map: dict = {}
        for key, mdef in self._materials.items():
            materials_map[key] = helpers.build_material(festim, key, mdef)

        logger.info(
            f"FestimProvider: running sim_type='{sim_type}' "
            f"with strategy {strategy_cls.__name__}"
        )

        problem_kwargs = strategy_cls().build(
            festim, festim_model, materials_map, helpers
        )

        run_dir = self._resolve_run_dir()

        prev_cwd = os.getcwd()
        os.chdir(run_dir)
        try:
            problem = festim.HydrogenTransportProblem(**problem_kwargs)
            problem.initialise()
            problem.run()
            logger.info(f"FestimProvider: '{sim_type}' completed in '{run_dir}'")
        except Exception as exc:  # noqa: BLE001
            logger.exception(
                f"FestimProvider: '{sim_type}' failed in '{run_dir}': {exc}"
            )
            err = classify_run_error("festim", exc)
            return make_failed_output("festim", sim_type, run_dir, err)
        finally:
            os.chdir(prev_cwd)

        fields, artifacts, diagnostics = extract_results(
            problem_kwargs["exports"], run_dir
        )
        return EngineOutput(
            status="completed",
            engine="festim",
            sim_type=sim_type,
            fields=fields,
            artifacts=artifacts,
            diagnostics=diagnostics,
            provenance=OutputProvenance(),
            run_dir=str(run_dir),
        )


register_provider("festim", FestimProvider)
