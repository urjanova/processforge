"""OpenMC provider for Monte Carlo neutronics simulations.

Architecture (mirrors festim_provider.py three layers)
------------------------------------------------------
1. **Dataclasses** (``OpenMCNuclide``, ``OpenMCElement``, ``OpenMCMeshTally``,
   ``OpenMCSourceBox``, ``OpenMCSolverConfig``) parse the opaque ``solver_config``
   JSON dict into typed, self-documenting objects.

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
        OpenMCSimStrategy, OpenMCSolverConfig, OpenMCBuildHelpers,
        register_openmc_sim_type,
    )

    class MyVolumeSim(OpenMCSimStrategy):
        def build(self, openmc, solver_cfg, materials_map, helpers):
            ...
            return openmc.Materials(...), openmc.Geometry(...), settings, tallies

    register_openmc_sim_type("volume_dagmc", MyVolumeSim)
"""

from __future__ import annotations

import dataclasses
import math
import os
import pathlib
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, List, Optional

from loguru import logger

from .base import AbstractProvider
from .registry import register_provider

if TYPE_CHECKING:
    from processforge.types import (
        FlowsheetConfig,
        MaterialDef,
        OpenMCProviderConfig,
        SimulationResult,
        UnitConfig,
    )


# ---------------------------------------------------------------------------
# OpenMC-specific solver config dataclasses
# ---------------------------------------------------------------------------

@dataclass
class OpenMCNuclide:
    """A single nuclide entry used with ``material.add_nuclide()``.

    name:         OpenMC nuclide identifier, e.g. ``"Li7"``, ``"U235"``, ``"He4"``
    percent:      atom or weight fraction (non-negative)
    percent_type: ``"ao"`` (atom) | ``"wo"`` (weight)
    """

    name: str
    percent: float
    percent_type: str = "ao"

    @classmethod
    def from_dict(cls, d: dict) -> "OpenMCNuclide":
        return cls(name=d["name"], percent=d["percent"], percent_type=d.get("percent_type", "ao"))


@dataclass
class OpenMCElement:
    """A single element entry used with ``material.add_element()``.

    element:      Periodic-table symbol, e.g. ``"C"``, ``"Ni"``, ``"Fe"``
    percent:      atom or weight fraction (non-negative)
    percent_type: ``"ao"`` | ``"wo"``
    """

    element: str
    percent: float
    percent_type: str = "ao"

    @classmethod
    def from_dict(cls, d: dict) -> "OpenMCElement":
        return cls(
            element=d["element"],
            percent=d["percent"],
            percent_type=d.get("percent_type", "ao"),
        )


@dataclass
class OpenMCMeshTally:
    """A ``RegularMesh`` tally definition.

    tally_id:    Integer identifier — must be unique across all tallies.
    lower_left:  ``[x, y, z]`` of mesh lower-left corner in cm.
    upper_right: ``[x, y, z]`` of mesh upper-right corner in cm.
    dimension:   ``[nx, ny, nz]`` voxel counts.
    scores:      List of score strings, e.g. ``["flux", "fission"]``.
    nuclides:    Optional list of nuclide identifiers to score over.
    estimator:   ``"tracklength"`` (default) | ``"analog"`` | ``"collision"``
    """

    tally_id: int
    lower_left: List[float]
    upper_right: List[float]
    dimension: List[int]
    scores: List[str] = field(default_factory=lambda: ["flux"])
    name: Optional[str] = None
    nuclides: Optional[List[str]] = None
    estimator: str = "tracklength"

    @classmethod
    def from_dict(cls, d: dict) -> "OpenMCMeshTally":
        return cls(
            tally_id=d["tally_id"],
            lower_left=d["lower_left"],
            upper_right=d["upper_right"],
            dimension=d["dimension"],
            scores=d.get("scores", ["flux"]),
            name=d.get("name"),
            nuclides=d.get("nuclides"),
            estimator=d.get("estimator", "tracklength"),
        )


@dataclass
class OpenMCSourceBox:
    """Box source definition corresponding to ``openmc.stats.Box``.

    lower_left:       ``[x, y, z]`` lower corner in cm.
    upper_right:      ``[x, y, z]`` upper corner in cm.
    only_fissionable: If ``True``, only sample source sites in fissionable regions.
    """

    lower_left: List[float]
    upper_right: List[float]
    only_fissionable: bool = True

    @classmethod
    def from_dict(cls, d: dict) -> "OpenMCSourceBox":
        return cls(
            lower_left=d["lower_left"],
            upper_right=d["upper_right"],
            only_fissionable=d.get("only_fissionable", True),
        )


@dataclass
class OpenMCSolverConfig:
    """Typed representation of the ``solver_config`` block for OpenMC ``SolverUnit`` units.

    Parsed from the opaque ``solver_config`` JSON dict by
    ``OpenMCProvider.run_simulation()``.

    Example flowsheet JSON::

        "solver_config": {
            "dagmc_path": "geometry/msre.h5m",
            "batches": 20,
            "inactive": 5,
            "particles": 20000,
            "run_mode": "eigenvalue",
            "temperature_default": 900.0,
            "source_box": {
                "lower_left":  [-125, -125, 0],
                "upper_right": [ 125,  125, 500],
                "only_fissionable": true
            },
            "mesh_tallies": [
                {
                    "tally_id": 1,
                    "name": "flux_tally",
                    "lower_left":  [-125, -125, 90],
                    "upper_right": [ 125,  125, 110],
                    "dimension": [400, 400, 1],
                    "scores": ["flux", "fission"]
                }
            ]
        }
    """

    batches: int = 20
    inactive: int = 5
    particles: int = 1000
    run_mode: str = "eigenvalue"
    dagmc_path: Optional[str] = None
    source_box: Optional[OpenMCSourceBox] = None
    mesh_tallies: List[OpenMCMeshTally] = field(default_factory=list)
    temperature_default: Optional[float] = None
    cross_sections: Optional[str] = None

    @classmethod
    def from_dict(cls, d: dict) -> "OpenMCSolverConfig":
        """Parse a raw ``solver_config`` dict into typed OpenMC dataclasses."""
        source_box = None
        if "source_box" in d:
            source_box = OpenMCSourceBox.from_dict(d["source_box"])

        mesh_tallies = [OpenMCMeshTally.from_dict(t) for t in d.get("mesh_tallies", [])]

        return cls(
            batches=d.get("batches", 20),
            inactive=d.get("inactive", 5),
            particles=d.get("particles", 1000),
            run_mode=d.get("run_mode", "eigenvalue"),
            dagmc_path=d.get("dagmc_path"),
            source_box=source_box,
            mesh_tallies=mesh_tallies,
            temperature_default=d.get("temperature_default"),
            cross_sections=d.get("cross_sections"),
        )


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
    def build_mesh_tally(openmc, tally_cfg: OpenMCMeshTally) -> object:
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
    def build_source(openmc, source_box: OpenMCSourceBox) -> object:
        """Construct an ``openmc.IndependentSource`` with a ``Box`` spatial distribution."""
        space = openmc.stats.Box(
            source_box.lower_left,
            source_box.upper_right,
            only_fissionable=source_box.only_fissionable,
        )
        source = openmc.IndependentSource()
        source.space = space
        source.angle = openmc.stats.Isotropic()
        return source

    @staticmethod
    def build_settings(openmc, solver_cfg: OpenMCSolverConfig, source: object) -> object:
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
        openmc: object,
        solver_cfg: OpenMCSolverConfig,
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


# ---------------------------------------------------------------------------
# Built-in strategies
# ---------------------------------------------------------------------------

class _EigenvalueDAGMCStrategy(OpenMCSimStrategy):
    """Eigenvalue (k-eff) simulation using DAGMC CAD geometry from a ``.h5m`` file.

    Corresponds to the MSRE reference script:
    - Geometry: ``openmc.DAGMCUniverse(h5m_filepath)``
    - Source:   ``openmc.stats.Box`` + ``Isotropic`` angle
    - Settings: eigenvalue mode
    - Tallies:  ``RegularMesh`` tallies

    ``solver_config`` must contain ``dagmc_path`` and ``source_box``.
    """

    def build(
        self,
        openmc,
        solver_cfg: OpenMCSolverConfig,
        materials_map: dict,
        helpers: OpenMCBuildHelpers,
    ) -> tuple:
        if not solver_cfg.dagmc_path:
            raise ValueError(
                "sim_type='eigenvalue_dagmc' requires 'dagmc_path' in solver_config"
            )
        if solver_cfg.source_box is None:
            raise ValueError(
                "sim_type='eigenvalue_dagmc' requires 'source_box' in solver_config"
            )

        omc_materials = openmc.Materials(list(materials_map.values()))
        dagmc_univ = openmc.DAGMCUniverse(solver_cfg.dagmc_path)
        geometry = openmc.Geometry(dagmc_univ)
        source = helpers.build_source(openmc, solver_cfg.source_box)
        settings = helpers.build_settings(openmc, solver_cfg, source)

        tally_objs = [helpers.build_mesh_tally(openmc, t) for t in solver_cfg.mesh_tallies]
        tallies = openmc.Tallies(tally_objs) if tally_objs else openmc.Tallies()

        return omc_materials, geometry, settings, tallies


class _FixedSourceDAGMCStrategy(OpenMCSimStrategy):
    """Fixed-source simulation using DAGMC CAD geometry.

    Identical setup to :class:`_EigenvalueDAGMCStrategy` but forces
    ``run_mode = "fixed source"`` regardless of the ``solver_config`` value.
    """

    def build(
        self,
        openmc,
        solver_cfg: OpenMCSolverConfig,
        materials_map: dict,
        helpers: OpenMCBuildHelpers,
    ) -> tuple:
        if not solver_cfg.dagmc_path:
            raise ValueError(
                "sim_type='fixed_source_dagmc' requires 'dagmc_path' in solver_config"
            )
        if solver_cfg.source_box is None:
            raise ValueError(
                "sim_type='fixed_source_dagmc' requires 'source_box' in solver_config"
            )

        omc_materials = openmc.Materials(list(materials_map.values()))
        dagmc_univ = openmc.DAGMCUniverse(solver_cfg.dagmc_path)
        geometry = openmc.Geometry(dagmc_univ)
        source = helpers.build_source(openmc, solver_cfg.source_box)

        # Force run_mode to "fixed source" regardless of JSON value
        fixed_cfg = dataclasses.replace(solver_cfg, run_mode="fixed source")
        settings = helpers.build_settings(openmc, fixed_cfg, source)

        tally_objs = [helpers.build_mesh_tally(openmc, t) for t in solver_cfg.mesh_tallies]
        tallies = openmc.Tallies(tally_objs) if tally_objs else openmc.Tallies()

        return omc_materials, geometry, settings, tallies


# Register built-ins at module load
register_openmc_sim_type("eigenvalue_dagmc", _EigenvalueDAGMCStrategy)
register_openmc_sim_type("fixed_source_dagmc", _FixedSourceDAGMCStrategy)


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

        self._provider_output_dir = provider_config.output_dir
        self._cross_sections = provider_config.cross_sections
        self._materials = {}

        for mat_name, mat_def in flowsheet_config.materials.items():
            self._materials[mat_name] = mat_def

        self._initialized = True
        n_mats = len({id(v) for v in self._materials.values()})
        logger.info(
            f"OpenMCProvider initialized with {n_mats} material(s). "
            f"Registered sim_types: {sorted(_SIM_TYPE_REGISTRY)}"
        )

    def get_thermo_properties(self, stream: dict) -> dict:
        """Delegate stream thermodynamics to CoolProp.

        OpenMC handles neutronics only; enthalpy, Cp, and K-values use CoolProp.
        """
        from processforge.thermo import get_Cp_molar, get_enthalpy_molar, get_K_values

        z, T, P = stream["z"], stream["T"], stream["P"]
        return {
            "H": get_enthalpy_molar(z, T, P),
            "Cp": get_Cp_molar(z, T, P),
            "K_values": get_K_values(list(z.keys()), T, P),
        }

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

        solver_cfg = OpenMCSolverConfig.from_dict(unit_config.solver_config or {})

        # Build openmc.Material objects for ALL registry materials — DAGMC
        # assigns materials by name so every declared material must be present.
        helpers = OpenMCBuildHelpers()
        materials_map: dict = {}
        seen_ids: set = set()
        for key, mdef in self._materials.items():
            mid = id(mdef)
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
        xs_path = solver_cfg.cross_sections or self._cross_sections
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
        solver_cfg: OpenMCSolverConfig,
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

        sp._f.close()  # explicitly release HDF5 file handle
        return scalars, metadata


# Register the provider
register_provider("openmc", OpenMCProvider)
