"""OpenMC provider for Monte Carlo neutronics simulations.

Architecture (three-layer strategy pattern)
-------------------------------------------
1. **Pydantic models** (``SourcePoint``, ``MeshTallyConfig``, ``SolverConfig``) in
   :mod:`processforge.schemas.openmc.openmc_model` parse the opaque
   ``solver_config`` JSON dict into typed, validated objects.

2. **Strategy registry** (``OpenMCSimStrategy`` + ``register_openmc_sim_type``)
   Each ``sim_type`` string maps to a strategy class.  New simulation types are
   added by subclassing ``OpenMCSimStrategy`` and calling
   ``register_openmc_sim_type(name, MyStrategy)`` — no changes needed here.

3. **OpenMCProvider** (``AbstractProvider`` subclass)
   Owns material lookup, material validation, working-directory management,
   and simulation dispatch.  All OpenMC object construction lives in
   ``OpenMCBuildHelpers`` so strategy classes are concise and independently
   testable.

Adding a new sim_type::

    from processforge.providers.openmc_provider import (
        OpenMCSimStrategy, SolverConfig, OpenMCBuildHelpers,
        register_openmc_sim_type,
    )

    class MyCustomSim(OpenMCSimStrategy):
        def build(self, openmc, solver_cfg, materials_map, helpers):
            ...
            return openmc.Materials(...), openmc.Geometry(...), settings, tallies

    register_openmc_sim_type("my_custom_sim", MyCustomSim)
"""

from __future__ import annotations

import math
import os
import pathlib
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional

from loguru import logger

from .base import AbstractProvider
from .registry import register_provider
from processforge.schemas.openmc.openmc_model import MeshTallyConfig, SolverConfig, SourcePoint


def _resolve_omc_path(path: Optional[str]) -> Optional[str]:
    """Expand environment variables (e.g. ``${OPENMC_DATA_ROOT}``) in a path string."""
    if path is None:
        return None
    return os.path.expandvars(path)

if TYPE_CHECKING:
    from processforge.types import (
        FlowsheetConfig,
        MaterialDef,
        OpenMCProviderConfig,
        SimulationResult,
        UnitConfig,
    )


# Backward-compat alias
OpenMCSolverConfig = SolverConfig


# ---------------------------------------------------------------------------
# Build helpers
# ---------------------------------------------------------------------------

class OpenMCBuildHelpers:
    """Shared OpenMC model-building utilities injected into sim strategies.

    All methods are ``@staticmethod`` so strategies can call them without a
    provider instance, making them unit-testable in isolation.
    """

    @staticmethod
    def build_material(openmc, mat_name: str, mat_def) -> object:
        """Construct an ``openmc.Material`` from a :class:`~processforge.types.MaterialDef`.

        Reads ``density``, ``density_units``, ``temperature`` from typed fields.
        Calls ``add_nuclide()`` for each entry in ``mat_def.nuclides`` (framework field).
        Calls ``add_element()`` for each entry in ``mat_def.extra["elements"]`` (OpenMC-specific).
        """
        mat = openmc.Material(material_id=mat_def.id, name=mat_name)

        if mat_def.density is not None:
            units = mat_def.density_units or "g/cm3"
            mat.set_density(units, mat_def.density)

        if mat_def.temperature is not None:
            mat.temperature = mat_def.temperature

        # Nuclides from the typed framework field
        for nuc in mat_def.nuclides:
            mat.add_nuclide(nuc["name"], nuc["percent"], nuc.get("percent_type", "ao"))

        # Elements from provider-specific extra field
        for elem in mat_def.extra.get("elements", []):
            mat.add_element(
                elem["element"],
                elem["percent"],
                elem.get("percent_type", "ao"),
            )

        if mat_def.depletable:
            mat.depletable = True

        return mat

    @staticmethod
    def build_mesh_tally(openmc, tally_cfg: MeshTallyConfig) -> object:
        """Construct an ``openmc.Tally`` with a ``MeshFilter`` over a ``RegularMesh``."""
        mesh = openmc.RegularMesh()
        mesh.dimension = tally_cfg.dimension
        mesh.lower_left = tally_cfg.lower_left
        mesh.upper_right = tally_cfg.upper_right

        mesh_filter = openmc.MeshFilter(mesh)
        tally = openmc.Tally(tally_id=tally_cfg.tally_id, name=tally_cfg.name)
        tally.filters = [mesh_filter]
        tally.scores = tally_cfg.scores
        if tally_cfg.nuclides:
            tally.nuclides = tally_cfg.nuclides
        if tally_cfg.estimator:
            tally.estimator = tally_cfg.estimator
        return tally

    @staticmethod
    def build_point_source(openmc, source_point: SourcePoint) -> object:
        """Construct an ``openmc.IndependentSource`` at a single point."""
        source = openmc.IndependentSource()
        source.space = openmc.stats.Point(source_point.xyz)
        source.angle = openmc.stats.Isotropic()
        if source_point.energy_eV is not None:
            source.energy = openmc.stats.Discrete([source_point.energy_eV], [1.0])
        else:
            source.energy = openmc.stats.Watt()
        return source

    @staticmethod
    def build_settings(openmc, solver_cfg: SolverConfig, source: object) -> object:
        """Construct an ``openmc.Settings`` object from solver config."""
        settings = openmc.Settings()
        settings.batches = solver_cfg.batches
        settings.inactive = solver_cfg.inactive
        settings.particles = solver_cfg.particles
        settings.run_mode = solver_cfg.run_mode
        settings.source = [source]
        if solver_cfg.temperature_default is not None:
            settings.temperature = {"default": solver_cfg.temperature_default}
        return settings


# ---------------------------------------------------------------------------
# Simulation type strategy registry
# ---------------------------------------------------------------------------

class OpenMCSimStrategy(ABC):
    """Base class for an OpenMC simulation type.

    Each subclass encapsulates the model setup logic for one ``sim_type`` value.
    The ``build`` method receives all necessary objects and returns a 4-tuple
    ``(materials, geometry, settings, tallies)`` ready for XML export.

    Register subclasses with :func:`register_openmc_sim_type`.

    Example::

        class MyFixedSourceCSG(OpenMCSimStrategy):
            def build(self, openmc, solver_cfg, materials_map, helpers):
                ...
                return omc_materials, geometry, settings, tallies

        register_openmc_sim_type("fixed_source_csg", MyFixedSourceCSG)
    """

    @abstractmethod
    def build(
        self,
        openmc,
        solver_cfg: SolverConfig,
        materials_map: dict,
        helpers: OpenMCBuildHelpers,
    ) -> tuple:
        """Set up the OpenMC model and return ``(materials, geometry, settings, tallies)``.

        Args:
            openmc:        The ``openmc`` module.
            solver_cfg:    Typed solver config from the flowsheet JSON.
            materials_map: ``{mat_name: openmc.Material}`` built by the provider.
            helpers:       Shared :class:`OpenMCBuildHelpers` instance.

        Returns:
            4-tuple: ``(omc_materials, geometry, settings, tallies)``
        """


_SIM_TYPE_REGISTRY: dict = {}


def register_openmc_sim_type(name: str, strategy_cls: type) -> None:
    """Register an OpenMC simulation type by name.

    Args:
        name:         The ``sim_type`` string used in the flowsheet JSON.
        strategy_cls: A subclass of :class:`OpenMCSimStrategy`.
    """
    _SIM_TYPE_REGISTRY[name] = strategy_cls


def get_registered_sim_types() -> dict[str, type]:
    """Return a view of the currently registered OpenMC sim_type → strategy mapping.

    Use this function (rather than accessing ``_SIM_TYPE_REGISTRY`` directly) so
    that callers are insulated from internal implementation changes.
    """
    return dict(_SIM_TYPE_REGISTRY)


# ---------------------------------------------------------------------------
# Built-in strategies
# ---------------------------------------------------------------------------

class _FixedSourcePointStrategy(OpenMCSimStrategy):
    """Fixed-source point-source approximation using a simple CSG sphere geometry.

    No DAGMC file required. Models the plant as a single point source inside
    a homogeneous sphere of a given radius and material — useful for broad-scale
    dose/flux estimates where full CAD geometry is not needed.
    """

    def build(
        self,
        openmc,
        solver_cfg: SolverConfig,
        materials_map: dict,
        helpers: OpenMCBuildHelpers,
    ) -> tuple:
        if solver_cfg.source_point is None:
            raise ValueError(
                "sim_type='fixed_source_point' requires 'source_point' in solver_config"
            )

        fill_name = solver_cfg.point_source_material
        if fill_name is not None and fill_name not in materials_map:
            raise ValueError(
                f"point_source_material='{fill_name}' not found in materials map. "
                f"Available: {sorted(materials_map)}"
            )

        omc_materials = openmc.Materials(list(materials_map.values()))

        sphere = openmc.Sphere(r=solver_cfg.point_source_sphere_radius, boundary_type="vacuum")
        fill_mat = materials_map[fill_name] if fill_name else None
        geometry = openmc.Geometry([
            openmc.Cell(fill=fill_mat, region=-sphere),
            openmc.Cell(region=+sphere),
        ])

        source = helpers.build_point_source(openmc, solver_cfg.source_point)
        fixed_cfg = solver_cfg.model_copy(update={"run_mode": "fixed source"})
        settings = helpers.build_settings(openmc, fixed_cfg, source)

        tally_objs = [helpers.build_mesh_tally(openmc, t) for t in solver_cfg.mesh_tallies]
        tallies = openmc.Tallies(tally_objs) if tally_objs else openmc.Tallies()

        return omc_materials, geometry, settings, tallies


class _EigenvalueCSGStrategy(OpenMCSimStrategy):
    """Eigenvalue (criticality) simulation using a simple CSG sphere geometry.

    Computes k_eff for a homogeneous sphere of fissile material.  The initial
    fission source is seeded from ``source_point``; OpenMC converges it over
    ``inactive`` batches before accumulating statistics.

    No DAGMC file required — useful for quick criticality estimates of molten
    salt or other homogeneous fissile compositions.
    """

    def build(
        self,
        openmc,
        solver_cfg: SolverConfig,
        materials_map: dict,
        helpers: OpenMCBuildHelpers,
    ) -> tuple:
        if solver_cfg.source_point is None:
            raise ValueError(
                "sim_type='eigenvalue_csg' requires 'source_point' in solver_config"
            )

        fill_name = solver_cfg.point_source_material
        if fill_name is not None and fill_name not in materials_map:
            raise ValueError(
                f"point_source_material='{fill_name}' not found in materials map. "
                f"Available: {sorted(materials_map)}"
            )

        omc_materials = openmc.Materials(list(materials_map.values()))

        sphere = openmc.Sphere(r=solver_cfg.point_source_sphere_radius, boundary_type="vacuum")
        fill_mat = materials_map[fill_name] if fill_name else None
        geometry = openmc.Geometry([
            openmc.Cell(fill=fill_mat, region=-sphere),
            openmc.Cell(region=+sphere),
        ])

        source = helpers.build_point_source(openmc, solver_cfg.source_point)
        eigenvalue_cfg = solver_cfg.model_copy(update={"run_mode": "eigenvalue"})
        settings = helpers.build_settings(openmc, eigenvalue_cfg, source)

        tally_objs = [helpers.build_mesh_tally(openmc, t) for t in solver_cfg.mesh_tallies]
        tallies = openmc.Tallies(tally_objs) if tally_objs else openmc.Tallies()

        return omc_materials, geometry, settings, tallies


# Register built-ins at module load
register_openmc_sim_type("fixed_source_point", _FixedSourcePointStrategy)
register_openmc_sim_type("eigenvalue_csg", _EigenvalueCSGStrategy)


# ---------------------------------------------------------------------------
# Material validation constants
# ---------------------------------------------------------------------------

_VALID_DENSITY_UNITS = frozenset({
    "g/cm3", "g/cc", "kg/m3", "atom/b-cm", "atom/cm3", "sum", "macro"
})


# ---------------------------------------------------------------------------
# OpenMCProvider
# ---------------------------------------------------------------------------

class OpenMCProvider(AbstractProvider):
    """Monte Carlo neutronics provider backed by OpenMC.

    Responsibilities
    ----------------
    * **Material registry** — ``initialize()`` loads all material defs from the
      flowsheet and indexes them by ``id`` and name.
    * **Working-directory management** — OpenMC writes XML and HDF5 output to
      the current working directory.  The provider changes CWD to
            ``<provider_output_dir>`` before each run and
      restores it afterwards in a ``finally`` block.
    * **Cross-section override** — If ``provider_config.cross_sections`` or
      ``solver_cfg.cross_sections`` is set, ``OPENMC_CROSS_SECTIONS`` is
      temporarily overridden.
    * **Simulation dispatch** — ``run_simulation()`` looks up the strategy from
      ``_SIM_TYPE_REGISTRY`` by ``sim_type``, calls ``strategy.build()``, exports
      XML, calls ``openmc.run()``, then parses the statepoint file.
    * **Result extraction** — extracts ``k_eff`` (eigenvalue mode) and per-tally
      aggregate scalars (mean totals, representative std devs).
    """

    def __init__(self):
        self._materials: dict = {}
        self._provider_output_dir: str = "outputs/openmc"
        self._cross_sections: Optional[str] = None
        self._initialized: bool = False

    def initialize(
        self,
        provider_config: "OpenMCProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        """Verify OpenMC is installed, store config, build material registry.

        Material registry is keyed by material name.
        """
        try:
            import openmc  # noqa: F401
        except ImportError as exc:
            raise RuntimeError(
                "OpenMC is not installed. Install with: pip install openmc"
            ) from exc

        # Resolve the output dir. Expand ${...} (mirrors cross_sections handling),
        # then anchor a relative path under the run output root so container/CI
        # runs land on the mounted volume rather than inside the working dir.
        out_dir = os.path.expandvars(provider_config.output_dir)
        if not os.path.isabs(out_dir):
            root = os.environ.get("PROCESSFORGE_OUTPUT_DIR", "outputs")
            out_dir = os.path.join(root, out_dir)
        self._provider_output_dir = out_dir
        self._cross_sections = provider_config.cross_sections
        self._materials = {}

        for mat_name, mat_def in flowsheet_config.materials.items():
            self._materials[mat_name] = mat_def

        self._initialized = True
        n_mats = len({v.id for v in self._materials.values()})
        logger.info(
            f"OpenMCProvider initialized with {n_mats} material(s). "
            f"Registered sim_types: {sorted(_SIM_TYPE_REGISTRY)}"
        )

    def get_thermo_properties(self, stream: dict) -> dict:
        raise NotImplementedError("OpenMCProvider does not support stream thermodynamics.")

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        """Return ``None`` — OpenMC uses ``run_simulation`` via ``SolverUnit``."""
        return None

    def teardown(self) -> None:
        """Release provider state."""
        self._initialized = False

    @classmethod
    def validate_material(cls, mat_name: str, mat_def, unit_cfg) -> list:
        """Validate OpenMC-specific material properties.

        Rules
        -----
        * ``density`` and ``density_units`` are required.
        * ``density_units`` must be one of the OpenMC-recognised strings.
        * At least one of ``nuclides`` or ``elements`` must be non-empty.

        Returns:
            List of error strings (empty = valid).
        """
        errors = []

        if mat_def.density is None:
            errors.append(
                f"Material '{mat_name}' is missing 'density' (required for OpenMC)."
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

        return errors

    def run_simulation(self, unit_config: "UnitConfig", inlet: dict) -> "SimulationResult":
        """Build and run an OpenMC simulation from a typed ``UnitConfig``.

        Execution flow
        --------------
        1. Parse ``solver_config`` → ``OpenMCSolverConfig``
        2. Look up strategy from ``_SIM_TYPE_REGISTRY``
        3. Build ``openmc.Material`` objects for all materials in the registry
        4. Call ``strategy.build()`` → ``(materials, geometry, settings, tallies)``
        5. Create output directory and ``chdir`` into it
        6. Apply cross-section env override if needed
        7. Export XML and call ``openmc.run()``
        8. Parse statepoint → extract k_eff and tally scalars
        9. Restore CWD and return ``SimulationResult``
        """
        import openmc
        from processforge.types import SimulationResult

        sim_type = unit_config.sim_type
        strategy_cls = _SIM_TYPE_REGISTRY.get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"OpenMCProvider: unknown sim_type '{sim_type}'. "
                f"Register with register_openmc_sim_type(). "
                f"Built-in types: {sorted(_SIM_TYPE_REGISTRY)}"
            )

        solver_cfg = SolverConfig.model_validate(unit_config.solver_config or {})

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

        omc_materials, geometry, settings, tallies = (
            strategy_cls().build(openmc, solver_cfg, materials_map, helpers)
        )

        # --- Working directory management ---
        run_dir = pathlib.Path(self._provider_output_dir)
        run_dir.mkdir(parents=True, exist_ok=True)
        original_cwd = pathlib.Path.cwd()

        # --- Cross-section override ---
        xs_path = _resolve_omc_path(solver_cfg.cross_sections or self._cross_sections)
        original_xs = os.environ.get("OPENMC_CROSS_SECTIONS")
        if xs_path:
            os.environ["OPENMC_CROSS_SECTIONS"] = str(xs_path)

        try:
            os.chdir(run_dir)

            # Export model XML files
            omc_materials.export_to_xml()
            geometry.export_to_xml()
            settings.export_to_xml()
            if tallies:
                tallies.export_to_xml()

            openmc.run()
            logger.info(f"OpenMCProvider: '{sim_type}' completed in '{run_dir}'")

            scalars, metadata = self._extract_results(openmc, solver_cfg, run_dir)

        finally:
            os.chdir(original_cwd)
            if xs_path:
                if original_xs is None:
                    os.environ.pop("OPENMC_CROSS_SECTIONS", None)
                else:
                    os.environ["OPENMC_CROSS_SECTIONS"] = original_xs

        return SimulationResult(
            status="completed",
            sim_type=sim_type,
            scalars=scalars,
            metadata=metadata,
        )

    def _extract_results(
        self,
        openmc,
        solver_cfg: SolverConfig,
        run_dir: pathlib.Path,
    ) -> tuple:
        """Parse the statepoint file and return ``(scalars_dict, metadata_dict)``.

        Scalar extraction
        -----------------
        * ``k_eff``, ``k_eff_std_dev`` — eigenvalue mode only (``sp.keff.n/.s``)
        * Per mesh tally and score:
          ``tally_{id}_{score}_mean_total`` — sum of all bin means
          ``tally_{id}_{score}_std_dev``    — RMS of per-bin std devs

        ``metadata["statepoint_path"]`` holds the absolute HDF5 path for
        downstream post-processing.
        """
        import glob as _glob

        scalars: dict = {}
        metadata: dict = {"run_dir": str(run_dir.resolve())}

        sp_pattern = str(run_dir / f"statepoint.{solver_cfg.batches}.h5")
        sp_files = _glob.glob(sp_pattern)
        if not sp_files:
            logger.warning(
                f"OpenMCProvider: statepoint file not found at '{sp_pattern}'. "
                "Scalars will be empty."
            )
            return scalars, metadata

        sp_path = sp_files[0]
        metadata["statepoint_path"] = str(pathlib.Path(sp_path).resolve())

        sp = openmc.StatePoint(sp_path)

        # k_eff — only present in eigenvalue mode
        if sp.keff is not None:
            scalars["k_eff"] = float(sp.keff.n)
            scalars["k_eff_std_dev"] = float(sp.keff.s)
            logger.info(
                f"OpenMCProvider: k_eff = {scalars['k_eff']:.6f} "
                f"+/- {scalars['k_eff_std_dev']:.6f}"
            )

        # Tally scalars
        for tally_cfg in solver_cfg.mesh_tallies:
            try:
                tally = sp.get_tally(id=tally_cfg.tally_id)
            except Exception:  # noqa: BLE001
                logger.warning(
                    f"OpenMCProvider: tally id={tally_cfg.tally_id} not found "
                    "in statepoint. Skipping."
                )
                continue

            for score in tally_cfg.scores:
                try:
                    df = tally.get_pandas_dataframe(scores=[score])
                    means = df["mean"].values
                    std_devs = df["std. dev."].values
                    key_prefix = f"tally_{tally_cfg.tally_id}_{score}"
                    scalars[f"{key_prefix}_mean_total"] = float(means.sum())
                    scalars[f"{key_prefix}_std_dev"] = (
                        float(math.sqrt((std_devs ** 2).mean())) if len(std_devs) > 0 else 0.0
                    )
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        f"OpenMCProvider: could not extract score '{score}' "
                        f"from tally {tally_cfg.tally_id}: {exc}"
                    )

        del sp
        return scalars, metadata


# Register the provider
register_provider("openmc", OpenMCProvider)
