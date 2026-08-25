"""OpenMC provider for Monte Carlo neutronics simulations.

Architecture (three-layer strategy pattern)
-------------------------------------------
1. **Pydantic models** (``SourcePoint``, ``MeshTallyConfig``, ``OpenMCSetting``) in
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
        OpenMCSimStrategy, OpenMCSetting, OpenMCBuildHelpers,
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
import tempfile
import threading
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional

from loguru import logger

from .base import AbstractProvider
from .registry import register_provider
from processforge.schemas.openmc.openmc_model import (
    MeshTallyConfig,
    OpenMCSetting,
    PointSourceGeometryConfig,
    ReactorCoreGeometryConfig,
    RunMode,
    SourcePoint,
)
from processforge.types import (
    EngineOutput,
    OutputArtifact,
    OutputField,
    OutputProvenance,
)
from processforge.providers.errors import classify_run_error, make_failed_output
from processforge.units import Quantity


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
        UnitConfig,
    )


# Backward-compat alias (kept briefly during migration; ``OpenMCSetting`` is the
# canonical runtime settings model).
OpenMCSolverConfig = OpenMCSetting


# Serializes OpenMC runs. Every run mutates process-global state (the
# OPENMC_CROSS_SECTIONS env var read by the spawned `openmc` subprocess), so
# concurrent runs sharing a process — e.g. the provider HTTP server threadpool —
# must not interleave.
_OPENMC_RUN_LOCK = threading.Lock()


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
    def build_reactor_core(openmc, geo_cfg, materials_map: dict) -> object:
        """Build a nested-cylinder reactor-core geometry from a ``ReactorCoreGeometryConfig``.

        Concentric cylindrical shells (core → reflector → vessel → gap → structure)
        sharing the core height, each filled with its declared material.  The
        outermost cylindrical surface and the top/bottom planes are vacuum; inner
        surfaces are transmission so particles cross between regions.
        """
        z0 = -geo_cfg.core_height / 2.0
        z1 = geo_cfg.core_height / 2.0
        zbot = openmc.ZPlane(z0=z0)
        ztop = openmc.ZPlane(z0=z1)

        mats = [
            geo_cfg.core_material,
            geo_cfg.reflector_material,
            geo_cfg.vessel_material,
            geo_cfg.gap_material,
            geo_cfg.structure_material,
        ]
        thicks = [
            0.0,
            geo_cfg.reflector_thickness,
            geo_cfg.vessel_thickness,
            geo_cfg.gap_thickness,
            geo_cfg.structure_thickness,
        ]

        # Cumulative outer radii for the shells that actually have a material.
        radii: list = []
        mat_list: list = []
        for mat, th in zip(mats, thicks):
            if mat is None:
                continue
            r_outer = geo_cfg.core_radius if not radii else radii[-1] + th
            radii.append(r_outer)
            mat_list.append(mat)

        if not radii:
            raise ValueError("reactor_core geometry defines no material shells")

        cells = []
        prev_cyl = None
        for mat, r_outer in zip(mat_list, radii):
            zcyl = openmc.ZCylinder(r=r_outer, boundary_type="transmission")
            if prev_cyl is None:
                region = -zcyl & -ztop & +zbot
            else:
                region = +prev_cyl & -zcyl & -ztop & +zbot
            cells.append(openmc.Cell(fill=materials_map[mat], region=region))
            prev_cyl = zcyl

        # Vacuum boundary on the outermost surfaces.
        zbot.boundary_type = "vacuum"
        ztop.boundary_type = "vacuum"
        prev_cyl.boundary_type = "vacuum"
        outer_region = +prev_cyl | +ztop | -zbot
        cells.append(openmc.Cell(region=outer_region))

        return openmc.Geometry(cells)

    @staticmethod
    def build_settings(openmc, solver_cfg: OpenMCSetting, source: object) -> object:
        """Construct an ``openmc.Settings`` object from solver config."""
        settings = openmc.Settings()
        settings.batches = solver_cfg.batches
        settings.inactive = solver_cfg.inactive
        settings.particles = solver_cfg.particles
        # `run_mode` may be a RunMode enum (validated config) or a raw string
        # (strategies apply it via `model_copy(update=...)`, which pydantic v2
        # does not re-validate). Normalise to a plain string for OpenMC.
        run_mode = solver_cfg.run_mode
        settings.run_mode = (
            run_mode.value if isinstance(run_mode, RunMode) else str(run_mode)
        )
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

    Subclasses declare ``config_model`` — the Pydantic model validating the
    per-strategy ``geometry_config`` block of the flowsheet.  The provider
    validates that block (with the flowsheet materials in context) before
    calling ``build``.

    Register subclasses with :func:`register_openmc_sim_type`.

    Example::

        class MyFixedSourceCSG(OpenMCSimStrategy):
            config_model = MyGeometryConfig
            def build(self, openmc, solver_cfg, geometry_cfg, materials_map, helpers):
                ...
                return omc_materials, geometry, settings, tallies

        register_openmc_sim_type("fixed_source_csg", MyFixedSourceCSG)
    """

    #: Pydantic model validating the unit's ``geometry_config`` block.  ``None``
    #: means the strategy takes no geometry config.
    config_model: type | None = None

    @abstractmethod
    def build(
        self,
        openmc,
        solver_cfg: OpenMCSetting,
        geometry_cfg,
        materials_map: dict,
        helpers: OpenMCBuildHelpers,
    ) -> tuple:
        """Set up the OpenMC model and return ``(materials, geometry, settings, tallies)``.

        Args:
            openmc:        The ``openmc`` module.
            solver_cfg:    Typed shared settings from ``solver_config``.
            geometry_cfg:  Validated per-strategy geometry config (or ``None``).
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


def _resolved_tally_cfgs(solver_cfg: OpenMCSetting, geometry_cfg) -> list:
    """Return the mesh-tally configs to actually build.

    Uses ``solver_cfg.mesh_tallies`` when provided, otherwise synthesises a
    single default tally spanning the geometry bounding box so flowsheets don't
    have to specify one.
    """
    if solver_cfg.mesh_tallies:
        return list(solver_cfg.mesh_tallies)
    if isinstance(geometry_cfg, ReactorCoreGeometryConfig):
        outer_r = (
            geometry_cfg.core_radius
            + geometry_cfg.reflector_thickness
            + geometry_cfg.vessel_thickness
            + geometry_cfg.gap_thickness
            + geometry_cfg.structure_thickness
        )
        half_z = geometry_cfg.core_height / 2.0
        ll = [-outer_r, -outer_r, -half_z]
        ur = [outer_r, outer_r, half_z]
    else:
        r = getattr(geometry_cfg, "sphere_radius", 500.0) if geometry_cfg else 500.0
        ll = [-r, -r, -r]
        ur = [r, r, r]
    return [
        MeshTallyConfig(
            tally_id=1,
            name="flux_default",
            lower_left=ll,
            upper_right=ur,
            dimension=[50, 50, 50],
            scores=["flux", "fission"],
        )
    ]


# ---------------------------------------------------------------------------
# Built-in strategies
# ---------------------------------------------------------------------------


class _PointSourceSphereStrategy(OpenMCSimStrategy):
    """Shared build logic for CSG-sphere point-source simulations.

    Models the plant as a single point source inside a homogeneous sphere of a
    given radius and material. Subclasses differ only in ``sim_type``
    and ``run_mode``. No DAGMC file required.
    """

    config_model = PointSourceGeometryConfig
    sim_type: str = "point_source_csg"
    run_mode: str = "eigenvalue"

    def build(
        self,
        openmc,
        solver_cfg: OpenMCSetting,
        geometry_cfg: PointSourceGeometryConfig,
        materials_map: dict,
        helpers: OpenMCBuildHelpers,
    ) -> tuple:
        if geometry_cfg is None or geometry_cfg.source_point is None:
            raise ValueError(
                f"sim_type='{self.sim_type}' requires 'source_point' in geometry_config"
            )

        fill_name = geometry_cfg.sphere_material
        if fill_name is not None and fill_name not in materials_map:
            raise ValueError(
                f"sphere_material='{fill_name}' not found in materials map. "
                f"Available: {sorted(materials_map)}"
            )

        omc_materials = openmc.Materials(list(materials_map.values()))

        sphere = openmc.Sphere(r=geometry_cfg.sphere_radius, boundary_type="vacuum")
        fill_mat = materials_map[fill_name] if fill_name else None
        geometry = openmc.Geometry(
            [
                openmc.Cell(fill=fill_mat, region=-sphere),
                openmc.Cell(region=+sphere),
            ]
        )

        source = helpers.build_point_source(openmc, geometry_cfg.source_point)
        run_cfg = solver_cfg.model_copy(update={"run_mode": self.run_mode})
        settings = helpers.build_settings(openmc, run_cfg, source)

        tally_objs = [
            helpers.build_mesh_tally(openmc, t)
            for t in _resolved_tally_cfgs(solver_cfg, geometry_cfg)
        ]
        tallies = openmc.Tallies(tally_objs)

        return omc_materials, geometry, settings, tallies


class _FixedSourcePointStrategy(_PointSourceSphereStrategy):
    """Fixed-source point-source approximation using a simple CSG sphere geometry.

    No DAGMC file required. Models the plant as a single point source inside
    a homogeneous sphere of a given radius and material — useful for broad-scale
    dose/flux estimates where full CAD geometry is not needed.
    """

    sim_type = "fixed_source_point"
    run_mode = "fixed source"


class _EigenvalueCSGStrategy(_PointSourceSphereStrategy):
    """Eigenvalue (criticality) simulation using a simple CSG sphere geometry.

    Computes k_eff for a homogeneous sphere of fissile material.  The initial
    fission source is seeded from ``source_point``; OpenMC converges it over
    ``inactive`` batches before accumulating statistics.

    No DAGMC file required — useful for quick criticality estimates of molten
    salt or other homogeneous fissile compositions.
    """

    sim_type = "eigenvalue_csg"
    run_mode = "eigenvalue"


class _ReactorCoreStrategy(OpenMCSimStrategy):
    """Approximate reactor-core simulation from a ``reactor_core`` geometry block.

    Builds nested cylindrical shells (core + reflector + vessel + gap +
    structure) from the declared flowsheet materials, so all of them
    participate — giving a neutronics result much closer to a real (DAGMC)
    reactor than the single homogeneous sphere, without any CAD asset.
    """

    config_model = ReactorCoreGeometryConfig
    sim_type: str = "reactor_core"
    run_mode: str = "eigenvalue"

    def build(
        self,
        openmc,
        solver_cfg: OpenMCSetting,
        geometry_cfg: ReactorCoreGeometryConfig,
        materials_map: dict,
        helpers: OpenMCBuildHelpers,
    ) -> tuple:
        if geometry_cfg is None:
            raise ValueError(
                f"sim_type='{self.sim_type}' requires a 'geometry_config' block"
            )
        if geometry_cfg.source_point is None:
            raise ValueError(
                f"sim_type='{self.sim_type}' requires 'source_point' in geometry_config"
            )

        omc_materials = openmc.Materials(list(materials_map.values()))
        geometry = helpers.build_reactor_core(openmc, geometry_cfg, materials_map)

        source = helpers.build_point_source(openmc, geometry_cfg.source_point)
        run_cfg = solver_cfg.model_copy(update={"run_mode": self.run_mode})
        settings = helpers.build_settings(openmc, run_cfg, source)

        tally_objs = [
            helpers.build_mesh_tally(openmc, t)
            for t in _resolved_tally_cfgs(solver_cfg, geometry_cfg)
        ]
        tallies = openmc.Tallies(tally_objs)

        return omc_materials, geometry, settings, tallies


class _FixedSourceReactorCoreStrategy(_ReactorCoreStrategy):
    """Fixed-source variant of the approximate reactor-core geometry."""

    sim_type = "fixed_source_reactor_core"
    run_mode = "fixed source"


# Register built-ins at module load
register_openmc_sim_type("fixed_source_point", _FixedSourcePointStrategy)
register_openmc_sim_type("eigenvalue_csg", _EigenvalueCSGStrategy)
register_openmc_sim_type("eigenvalue_reactor", _ReactorCoreStrategy)
register_openmc_sim_type("fixed_source_reactor_core", _FixedSourceReactorCoreStrategy)


# ---------------------------------------------------------------------------
# Material validation constants
# ---------------------------------------------------------------------------

_VALID_DENSITY_UNITS = frozenset(
    {"g/cm3", "g/cc", "kg/m3", "atom/b-cm", "atom/cm3", "sum", "macro"}
)


# ---------------------------------------------------------------------------
# OpenMCProvider
# ---------------------------------------------------------------------------


class OpenMCProvider(AbstractProvider):
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
      ``_SIM_TYPE_REGISTRY`` by ``sim_type``, calls ``strategy.build()``, runs
      via ``Model.run(cwd=...)``, then parses the returned statepoint.
      Failures are returned as ``SimulationResult(status="failed")`` with
      ``run_dir`` and error metadata rather than silently bubbling up.
    * **Result extraction** — extracts ``k_eff`` (eigenvalue mode) and per-tally
      aggregate scalars (mean totals, representative std devs). Extraction
      issues surface in ``metadata["tally_warnings"]``.
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
            f"Registered sim_types: {sorted(_SIM_TYPE_REGISTRY)}"
        )

    def get_thermo_properties(self, stream: dict) -> dict:
        raise NotImplementedError(
            "OpenMCProvider does not support stream thermodynamics."
        )

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
        * ``density`` is required unless ``density_units`` is ``"sum"``.
        * ``density_units`` must be one of the OpenMC-recognised strings.
        * At least one of ``nuclides`` or ``elements`` must be non-empty.
        * Each ``nuclides`` entry needs ``name`` + numeric ``percent``;
          each ``elements`` entry needs ``element`` + numeric ``percent``.
        * ``percent_type`` must be ``"ao"`` or ``"wo"``.

        Returns:
            List of error strings (empty = valid).
        """
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
        """Build and run an OpenMC simulation from a typed ``UnitConfig``.

        Execution flow
        --------------
        1. Parse ``solver_config`` → ``OpenMCSetting``
        2. Look up strategy from ``_SIM_TYPE_REGISTRY``
        3. Build ``openmc.Material`` objects for all materials in the registry
        4. Call ``strategy.build()`` → ``(materials, geometry, settings, tallies)``
        5. Resolve the run directory (with temp-dir fallback)
        6. Apply the cross-section override inside a module-wide lock
        7. Run via ``openmc.model.Model.run(cwd=run_dir, ...)`` (no ``chdir``)
        8. Parse the returned statepoint → extract k_eff and tally scalars
        9. Return :class:`EngineOutput`; failures surface as ``status="failed"``
           with ``run_dir`` and error diagnostics in ``diagnostics``
        """
        import openmc
        from processforge.types import (
            EngineOutput,
            OutputArtifact,
            OutputField,
            OutputProvenance,
            Quantity,
        )

        sim_type = unit_config.sim_type
        strategy_cls = _SIM_TYPE_REGISTRY.get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"OpenMCProvider: unknown sim_type '{sim_type}'. "
                f"Register with register_openmc_sim_type(). "
                f"Built-in types: {sorted(_SIM_TYPE_REGISTRY)}"
            )

        solver_cfg = OpenMCSetting.model_validate(unit_config.solver_config or {})

        # Validate the per-strategy geometry_config (if the strategy declares one).
        # Flowsheet materials are supplied via context so the config_model can
        # check that every referenced material actually exists.
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
        # export+run window is serialized (see module-level comment).
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

        fields, artifacts, diagnostics = self._extract_results(
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
            )

        return EngineOutput(
            status="completed",
            engine="openmc",
            sim_type=sim_type,
            fields=fields,
            artifacts=artifacts,
            diagnostics=diagnostics,
            provenance=OutputProvenance(),
        )

    def _resolve_run_dir(self) -> pathlib.Path:
        """Create the run output directory, falling back to a temp dir on failure.

        The configured output dir (often the mounted /data volume) may not be
        writable by the container user. Fall back to a temp dir so the
        simulation can still run; only the persisted XML artifacts are lost.
        """
        run_dir = pathlib.Path(self._provider_output_dir)
        try:
            run_dir.mkdir(parents=True, exist_ok=True)
        except (PermissionError, OSError) as exc:
            run_dir = pathlib.Path(tempfile.mkdtemp(prefix="processforge_openmc_"))
            logger.warning(
                f"OpenMCProvider: cannot use output dir "
                f"'{self._provider_output_dir}' ({exc}); "
                f"falling back to '{run_dir}'."
            )
        return run_dir

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

    # Score → SI-ish unit string, so outputs are comparable with other engines
    # (CoolProp J/mol, FESTIM W/m², …) via the shared ``Quantity``/pint layer.
    _SCORE_UNITS = {
        "flux": "n/cm^2/s",
        "fission": "1/cm^3/s",
        "heating": "W",
        "kappa_fission": "W",
        "absorption": "1/cm^3/s",
        "scatter": "1/cm^3/s",
        "elastic": "1/cm^3/s",
        "capture": "1/cm^3/s",
        "total": "1/cm^3/s",
        "current": "1/cm^2/s",
        "delayed_nu_fission": "1/cm^3/s",
    }
    # Energy per fission used to convert a fission rate into a comparable power.
    _MEV_PER_FISSION = 200.0
    _J_PER_MEV = 1.602176634e-13  # 1 MeV = 1e6 eV * 1.602e-19 J/eV

    def _extract_results(
        self,
        openmc,
        solver_cfg: OpenMCSetting,
        tally_cfgs: list,
        statepoint_path,
        run_dir: pathlib.Path,
    ) -> tuple:
        """Parse the statepoint and return ``(fields, artifacts, diagnostics)``.

        Standardized, unit-bearing fields
        ----------------------------------
        * ``k_eff`` (dimensionless) with ``std_dev`` — eigenvalue mode only
        * Per mesh tally + score:
          ``tally_{id}_{score}_mean``       — volume-averaged bin mean
          ``tally_{id}_{score}_integrated`` — Σ(bin mean × cell volume)
          both tagged with the score's unit (e.g. ``flux`` → ``n/cm^2/s``).
        * ``power`` (W) — derived from the total fission rate
          (``fission`` score) using 200 MeV/fission, so OpenMC outputs are
          directly comparable with CoolProp/FESTIM thermal quantities.

        The raw per-voxel tally dataframe is also written to a CSV
        :class:`OutputArtifact`, and the statepoint HDF5 is registered as an
        artifact. Extraction issues are collected in ``diagnostics["tally_warnings"]``.
        """
        from processforge.types import OutputArtifact, OutputField, Quantity

        fields: list = []
        artifacts: list = []
        diagnostics: dict = {"run_dir": str(run_dir.resolve())}
        warnings: list = []

        if statepoint_path is None:
            warnings.append("run produced no statepoint file")
            diagnostics["tally_warnings"] = warnings
            return fields, artifacts, diagnostics

        sp_path = str(pathlib.Path(statepoint_path).resolve())
        artifacts.append(OutputArtifact(
            name="statepoint",
            kind="statepoint",
            local_path=sp_path,
            source="local",
        ))
        diagnostics["statepoint_path"] = sp_path

        sp = openmc.StatePoint(statepoint_path)
        try:
            # k_eff — only present in eigenvalue mode
            if sp.keff is not None:
                fields.append(OutputField(
                    name="k_eff",
                    quantity=Quantity(
                        value=float(sp.keff.n), unit="", std_dev=float(sp.keff.s)
                    ),
                    kind="scalar",
                    source="keff",
                ))
                logger.info(
                    f"OpenMCProvider: k_eff = {float(sp.keff.n):.6f} "
                    f"+/- {float(sp.keff.s):.6f}"
                )
            elif solver_cfg.run_mode == RunMode.eigenvalue:
                warnings.append("eigenvalue run produced no k_eff in statepoint")

            for tally_cfg in tally_cfgs:
                try:
                    tally = sp.get_tally(id=tally_cfg.tally_id)
                except Exception as exc:  # noqa: BLE001
                    warnings.append(
                        f"tally id={tally_cfg.tally_id} not found in statepoint: {exc}"
                    )
                    continue

                # Cell volume (area for 2D meshes) for volume-weighted aggregation.
                try:
                    ll = [float(x) for x in tally.mesh.lower_left]
                    ur = [float(x) for x in tally.mesh.upper_right]
                    dim = [float(x) for x in tally.mesh.dimension]
                    cell_vol = 1.0
                    for a, b, d in zip(ll, ur, dim):
                        cell_vol *= (b - a) / d if d else 1.0
                except Exception:  # noqa: BLE001
                    cell_vol = 1.0

                for score in tally_cfg.scores:
                    try:
                        df = tally.get_pandas_dataframe(scores=[score])
                        means = [float(x) for x in df["mean"].values]
                        std_devs = [float(x) for x in df["std. dev."].values]
                        unit = self._SCORE_UNITS.get(score, "")
                        key_prefix = f"tally_{tally_cfg.tally_id}_{score}"

                        mean_val = sum(means) / len(means) if means else 0.0
                        integrated = sum(m * cell_vol for m in means)
                        integrated_std = math.sqrt(
                            sum((s * cell_vol) ** 2 for s in std_devs)
                        ) if std_devs else 0.0

                        fields.append(OutputField(
                            name=f"{key_prefix}_mean",
                            quantity=Quantity(value=mean_val, unit=unit),
                            kind="scalar",
                            source=f"tally_{tally_cfg.tally_id}/{score}",
                        ))
                        fields.append(OutputField(
                            name=f"{key_prefix}_integrated",
                            quantity=Quantity(
                                value=integrated, unit=unit, std_dev=integrated_std
                            ),
                            kind="scalar",
                            source=f"tally_{tally_cfg.tally_id}/{score}",
                        ))

                        # Total fission rate → comparable power (W).
                        if score == "fission" and integrated > 0:
                            power_W = (
                                integrated
                                * self._MEV_PER_FISSION
                                * self._J_PER_MEV
                            )
                            power_std = (
                                power_W * (integrated_std / integrated)
                                if integrated_std else 0.0
                            )
                            fields.append(OutputField(
                                name="power",
                                quantity=Quantity(
                                    value=power_W, unit="W", std_dev=power_std
                                ),
                                kind="scalar",
                                source="fission_rate",
                            ))
                            diagnostics.setdefault("notes", []).append(
                                "power derived from total fission rate assuming "
                                f"{self._MEV_PER_FISSION} MeV/fission."
                            )
                    except Exception as exc:  # noqa: BLE001
                        warnings.append(
                            f"could not extract score '{score}' from tally "
                            f"{tally_cfg.tally_id}: {exc}"
                        )

                # Raw per-voxel field as a CSV artifact for downstream post-processing.
                try:
                    csv_name = f"tally_{tally_cfg.tally_id}_{tally_cfg.name or tally_cfg.tally_id}.csv"
                    csv_path = pathlib.Path(run_dir) / csv_name
                    df_all = tally.get_pandas_dataframe()
                    df_all.to_csv(csv_path)
                    artifacts.append(OutputArtifact(
                        name=csv_name,
                        kind="csv",
                        local_path=str(csv_path.resolve()),
                        source="local",
                    ))
                except Exception as exc:  # noqa: BLE001
                    warnings.append(
                        f"could not export tally {tally_cfg.tally_id} CSV: {exc}"
                    )
        finally:
            del sp

        if warnings:
            diagnostics["tally_warnings"] = warnings
        return fields, artifacts, diagnostics


# Register the provider
register_provider("openmc", OpenMCProvider)
