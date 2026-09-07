"""OpenMC provider for Monte Carlo neutronics simulations."""
from __future__ import annotations

import os
import pathlib
import threading
from typing import TYPE_CHECKING, Optional

from loguru import logger

from processforge.providers.base_simulation_provider import BaseSimulationProvider
from processforge.providers.errors import ProviderNotAvailableError, classify_run_error, make_failed_output
from processforge.providers.registry import register_provider
from processforge.schemas.openmc.openmc_model import OpenMCSetting
from processforge.types import EngineOutput, OutputProvenance

from .build_helpers import OpenMCBuildHelpers
from .result_extraction import extract_results
from .strategies import _resolved_tally_cfgs, get_registered_sim_types

if TYPE_CHECKING:
    from processforge.types import FlowsheetConfig, MaterialDef, OpenMCProviderConfig, UnitConfig


def _resolve_omc_path(path: Optional[str]) -> Optional[str]:
    """Expand environment variables (e.g. ``${OPENMC_DATA_ROOT}``) in a path string."""
    if path is None:
        return None
    return os.path.expandvars(path)


# Backward-compat alias (kept briefly during migration; ``OpenMCSetting`` is the
# canonical runtime settings model).
OpenMCSolverConfig = OpenMCSetting


# Serializes OpenMC runs. Every run mutates process-global state (the
# OPENMC_CROSS_SECTIONS env var read by the spawned `openmc` subprocess), so
# concurrent runs sharing a process — e.g. the provider HTTP server threadpool —
# must not interleave.
_OPENMC_RUN_LOCK = threading.Lock()


# Material validation constants
_VALID_DENSITY_UNITS = frozenset(
    {"g/cm3", "g/cc", "kg/m3", "atom/b-cm", "atom/cm3", "sum", "macro"}
)


class OpenMCProvider(BaseSimulationProvider):
    """Monte Carlo neutronics provider backed by OpenMC.

    Responsibilities
    ----------------
    * **Material registry** — ``initialize()`` loads all material defs from the
      flowsheet and indexes them by ``id`` and name.
    * **Working-directory management** — each run executes in
      ``<provider_output_dir>`` via ``openmc.model.Model.run(cwd=...)`` — no
      process-wide ``chdir``. Falls back to a temp dir when the configured
      output dir is not writable.
    * **Cross-section override** — If ``provider_config.cross_sections`` or
      ``solver_cfg.cross_sections`` is set, ``OPENMC_CROSS_SECTIONS`` is
      temporarily overridden. Runs are serialized behind a module-wide lock so
      the process-global env var cannot race across concurrent requests.
    * **Simulation dispatch** — ``run_simulation()`` looks up the strategy from
      the shared OpenMC strategy registry by ``sim_type``, calls
      ``strategy.build()``, runs via ``Model.run(cwd=...)``, then parses the
      returned statepoint. Failures are returned as ``EngineOutput(status="failed")``
      with ``run_dir`` and error metadata.
    * **Result extraction** — extracts ``k_eff`` (eigenvalue mode) and per-tally
      aggregate scalars. Extraction issues surface in ``metadata["tally_warnings"]``.
    """

    def __init__(self):
        super().__init__()
        self._provider_output_dir: str = "outputs/openmc"
        self._cross_sections: Optional[str] = None

    def initialize(
        self,
        provider_config: "OpenMCProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        """Verify OpenMC is installed, store config, build material registry."""
        try:
            import openmc  # noqa: F401
        except ImportError as exc:
            raise ProviderNotAvailableError(
                "OpenMC is not installed. Install with: pip install openmc"
            ) from exc

        # Resolve the output dir. Expand ${...} (mirrors cross_sections handling),
        # then anchor a relative path under the run output root so container/CI
        # runs land on the mounted volume rather than inside the working dir.
        self._provider_output_dir = self._expand_output_dir(provider_config.output_dir)
        self._cross_sections = provider_config.cross_sections

        # Material IDs must be unique — duplicate IDs were previously dropped
        # silently, causing missing-material KeyErrors deep in the run.
        id_to_name: dict = {}
        duplicates = []
        for mat_name, mat_def in flowsheet_config.materials.items():
            mid = mat_def.id
            if mid in id_to_name:
                duplicates.append(
                    f"id={mid} shared by '{id_to_name[mid]}' and '{mat_name}'"
                )
            else:
                id_to_name[mid] = mat_name
        if duplicates:
            raise ValueError(
                "Flowsheet materials have duplicate OpenMC ids: "
                + "; ".join(duplicates)
            )

        # Fail fast on a bad cross-section path instead of deep inside openmc.
        resolved_xs = _resolve_omc_path(self._cross_sections)
        if resolved_xs:
            self._check_cross_sections_path(resolved_xs)

        self._materials = dict(flowsheet_config.materials.items())

        self._initialized = True
        n_mats = len({v.id for v in self._materials.values()})
        logger.info(
            f"OpenMCProvider initialized with {n_mats} material(s). "
            f"Registered sim_types: {sorted(get_registered_sim_types())}"
        )

    def teardown(self) -> None:
        """Release provider state."""
        super().teardown()

    @classmethod
    def validate_material(cls, mat_name: str, mat_def, unit_cfg) -> list:
        """Validate OpenMC-specific material properties."""
        errors = []

        if mat_def.density is None and mat_def.density_units != "sum":
            errors.append(
                f"Material '{mat_name}' is missing 'density' (required for OpenMC; "
                "only 'sum' density_units may omit it)."
            )
        if mat_def.density_units is None:
            errors.append(
                f"Material '{mat_name}' is missing 'density_units' (required for OpenMC)."
            )
        elif mat_def.density_units not in _VALID_DENSITY_UNITS:
            errors.append(
                f"Material '{mat_name}' has invalid density_units "
                f"'{mat_def.density_units}'. Valid: {sorted(_VALID_DENSITY_UNITS)}"
            )

        has_nuclides = bool(mat_def.nuclides)
        has_elements = bool(mat_def.extra.get("elements"))
        if not has_nuclides and not has_elements:
            errors.append(
                f"Material '{mat_name}' has no nuclides or elements. "
                "Add at least one 'nuclides' or 'elements' entry."
            )

        errors.extend(
            cls._validate_component_entries(mat_name, "nuclides", mat_def.nuclides)
        )
        errors.extend(
            cls._validate_component_entries(
                mat_name, "elements", mat_def.extra.get("elements", [])
            )
        )

        return errors

    @staticmethod
    def _validate_component_entries(mat_name: str, kind: str, entries: list) -> list:
        """Validate each nuclide/element component entry in a material."""
        label_key = "name" if kind == "nuclides" else "element"
        errors = []
        for i, entry in enumerate(entries or []):
            if not isinstance(entry, dict):
                errors.append(
                    f"Material '{mat_name}' {kind}[{i}] must be an object, "
                    f"got {type(entry).__name__}."
                )
                continue
            if not entry.get(label_key):
                errors.append(
                    f"Material '{mat_name}' {kind}[{i}] is missing '{label_key}'."
                )
            if not isinstance(entry.get("percent"), (int, float)):
                errors.append(
                    f"Material '{mat_name}' {kind}[{i}] must have a numeric 'percent'."
                )
            ptype = entry.get("percent_type", "ao")
            if ptype not in ("ao", "wo"):
                errors.append(
                    f"Material '{mat_name}' {kind}[{i}] percent_type must be "
                    f"'ao' or 'wo', got '{ptype}'."
                )
        return errors

    def run_simulation(
        self, unit_config: "UnitConfig", inlet: dict
    ) -> "EngineOutput":
        """Build and run an OpenMC simulation from a typed ``UnitConfig``."""
        import openmc

        sim_type = unit_config.sim_type
        strategy_cls = get_registered_sim_types().get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"OpenMCProvider: unknown sim_type '{sim_type}'. "
                f"Register with register_openmc_sim_type(). "
                f"Built-in types: {sorted(get_registered_sim_types())}"
            )

        solver_cfg = OpenMCSetting.model_validate(unit_config.solver_config or {})

        # Validate the per-strategy geometry_config (if the strategy declares one).
        geometry_cfg = None
        cfg_model = getattr(strategy_cls, "config_model", None)
        if cfg_model is not None:
            geometry_cfg = cfg_model.model_validate(
                unit_config.geometry_config or {},
                context={"materials": set(self._materials.keys())},
            )

        # Build openmc.Material objects for all registry materials.
        helpers = OpenMCBuildHelpers()
        materials_map: dict = {}
        seen_ids: set = set()
        for key, mdef in self._materials.items():
            mid = mdef.id
            if mid in seen_ids:
                continue
            seen_ids.add(mid)
            omc_mat = helpers.build_material(openmc, key, mdef)
            materials_map[key] = omc_mat

        logger.info(
            f"OpenMCProvider: running sim_type='{sim_type}' "
            f"with strategy {strategy_cls.__name__}, "
            f"{solver_cfg.batches} batches × {solver_cfg.particles} particles"
        )

        omc_materials, geometry, settings, tallies = strategy_cls().build(
            openmc, solver_cfg, geometry_cfg, materials_map, helpers
        )
        resolved_tally_cfgs = _resolved_tally_cfgs(solver_cfg, geometry_cfg)

        run_dir = self._resolve_run_dir()

        xs_path = _resolve_omc_path(solver_cfg.cross_sections or self._cross_sections)
        if xs_path:
            self._check_cross_sections_path(xs_path)

        # Runs spawn a fresh `openmc` subprocess that reads OPENMC_CROSS_SECTIONS
        # from the environment. That env var is process-global, so the whole
        # export+run window is serialized.
        with _OPENMC_RUN_LOCK:
            original_xs = os.environ.get("OPENMC_CROSS_SECTIONS")
            if xs_path:
                os.environ["OPENMC_CROSS_SECTIONS"] = str(xs_path)
            try:
                model = openmc.model.Model(
                    geometry=geometry,
                    materials=omc_materials,
                    settings=settings,
                    tallies=tallies,
                )
                statepoint_path = model.run(cwd=str(run_dir), export_model_xml=False)
                logger.info(f"OpenMCProvider: '{sim_type}' completed in '{run_dir}'")
            except Exception as exc:  # noqa: BLE001
                logger.exception(
                    f"OpenMCProvider: '{sim_type}' failed in '{run_dir}': {exc}"
                )
                err = classify_run_error("openmc", exc)
                return make_failed_output("openmc", sim_type, run_dir, err)
            finally:
                if xs_path:
                    if original_xs is None:
                        os.environ.pop("OPENMC_CROSS_SECTIONS", None)
                    else:
                        os.environ["OPENMC_CROSS_SECTIONS"] = original_xs

        fields, artifacts, diagnostics = extract_results(
            openmc, solver_cfg, resolved_tally_cfgs, statepoint_path, run_dir
        )

        if statepoint_path is None:
            logger.warning(
                f"OpenMCProvider: '{sim_type}' produced no statepoint in '{run_dir}'."
            )
            diagnostics.setdefault("error", "run produced no statepoint file")
            return EngineOutput(
                status="failed",
                engine="openmc",
                sim_type=sim_type,
                fields=fields,
                artifacts=artifacts,
                diagnostics=diagnostics,
                run_dir=str(run_dir),
            )

        return EngineOutput(
            status="completed",
            engine="openmc",
            sim_type=sim_type,
            fields=fields,
            artifacts=artifacts,
            diagnostics=diagnostics,
            provenance=OutputProvenance(),
            run_dir=str(run_dir),
        )

    @staticmethod
    def _check_cross_sections_path(xs_path: str) -> None:
        """Fail fast if a cross-section path does not exist."""
        if pathlib.Path(xs_path).is_file():
            return
        raise RuntimeError(
            f"OpenMC cross-sections file not found at '{xs_path}'. "
            "Set a valid OPENMC_CROSS_SECTIONS / provider 'cross_sections' / "
            "solver_config 'cross_sections' path."
        )


register_provider("openmc", OpenMCProvider)
