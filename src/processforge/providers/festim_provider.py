"""FESTIM provider — in-process hydrogen transport FEM solver.

Architecture (three-layer strategy pattern)
------------------------------------------
1. **FestimBuildHelpers** — shared utilities for constructing FESTIM objects
   (``Material``, ``Mesh1D``, ``Species``, ``Stepsize``, …) from the typed
   ``FestimModel`` schema.  Structured profiles (``GaussianProfile``,
   ``RampProfile``, ``AfterTProfile``) are converted back into UFL expressions
   / Python callables here.

2. **Strategy registry** (``FestimSimStrategy`` + ``register_festim_sim_type``)
   Each ``sim_type`` string maps to a strategy class.  New simulation types are
   added by subclassing ``FestimSimStrategy`` and calling
   ``register_festim_sim_type(name, MyStrategy)`` — no changes needed here.

3. **FestimProvider** (``AbstractProvider`` subclass)
   Imports ``festim`` lazily, builds the model parts via the strategy's
   ``build`` method, then runs ``HydrogenTransportProblem.initialise()`` +
   ``.run()`` in the resolved run directory.  FESTIM is a Docker-only provider
   (like OpenMC): this class is the engine executed inside the container by
   ``processforge.api.serve``, and the CLI only ever reaches it through
   ``ContainerProviderClient``.

Adding a new sim_type::

    from processforge.providers.festim_provider import (
        FestimSimStrategy, FestimBuildHelpers, register_festim_sim_type,
    )

    class MyCustomSim(FestimSimStrategy):
        def build(self, festim, festim_model, materials_map, helpers):
            ...
            return problem_kwargs

    register_festim_sim_type("my_custom_sim", MyCustomSim)
"""

from __future__ import annotations

import os
import pathlib
import tempfile
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from loguru import logger

from processforge.schemas.festim.festim_model import (
    AfterTProfile,
    FestimModel,
    FixedConcentrationBCConfig,
    GaussianProfile,
    HenrysBCConfig,
    ParticleFluxBCConfig,
    Profile1DExportConfig,
    RampProfile,
    SievertsBCConfig,
    SurfaceFluxExportConfig,
    TotalVolumeExportConfig,
)

from .base import AbstractProvider
from .registry import register_provider
from processforge.types import (
    EngineOutput,
    OutputArtifact,
    OutputField,
    OutputProvenance,
)
from processforge.units import Quantity

if TYPE_CHECKING:
    from processforge.types import (
        FestimProviderConfig,
        FlowsheetConfig,
        MaterialDef,
        UnitConfig,
    )


# ---------------------------------------------------------------------------
# Value-profile → callable builders
# ---------------------------------------------------------------------------


def _build_ramp_callable(profile: RampProfile):
    """Turn a :class:`RampProfile` into a ``t -> float`` callable."""

    def ramp(t):
        if t <= profile.start_t:
            return profile.T0
        return profile.T0 + profile.ramp_rate * (t - profile.start_t)

    return ramp


def _build_max_stepsize(value: float | AfterTProfile):
    """Turn a ``max_stepsize`` into a float or a ``t -> float | None`` callable."""
    if isinstance(value, AfterTProfile):

        def max_stepsize(t):
            return value.value if t > value.after_t else None

        return max_stepsize
    return float(value)


def _build_gaussian_value(profile: GaussianProfile):
    """Build a ``(x, t) -> UFL expr`` callable for a Gaussian particle source.

    Matches the FESTIM TDS tutorial: the gaussian is weighted by the incident
    flux amplitude and optionally gated off after ``active_until`` seconds via
    a ``ufl.conditional``.  ``ufl`` is imported lazily so the provider module
    stays importable without FESTIM/UFL installed.
    """

    def value(x, t):
        import ufl

        expr = (
            profile.amplitude
            / (profile.width * (2 * ufl.pi) ** 0.5)
            * ufl.exp(-0.5 * ((x[0] - profile.center) / profile.width) ** 2)
        )
        if profile.active_until is not None:
            expr = expr * ufl.conditional(ufl.le(t, profile.active_until), 1.0, 0.0)
        return expr

    return value


def _build_expression_value(
    value: float | GaussianProfile | RampProfile,
):
    """Resolve a structured value into a FESTIM-friendly float or callable."""
    if isinstance(value, GaussianProfile):
        return _build_gaussian_value(value)
    if isinstance(value, RampProfile):
        return _build_ramp_callable(value)
    return float(value)


# ---------------------------------------------------------------------------
# Build helpers
# ---------------------------------------------------------------------------


class FestimBuildHelpers:
    """Shared FESTIM object-construction utilities injected into sim strategies.

    All methods are ``@staticmethod`` so strategies can call them without a
    provider instance, making them unit-testable in isolation.  The ``festim``
    module is always passed in as the first argument (imported lazily by the
    provider), which is what lets tests substitute a ``fake_festim``.
    """

    # Optional FESTIM material properties read from ``extra``. FESTIM accepts
    # ``None`` for all of these *except* ``solubility_law`` (its setter raises on
    # ``None``), so that one needs a real default. Every field missing from
    # ``extra`` is logged as a WARNING so users know a default was applied.
    _MATERIAL_FIELD_DEFAULTS = {
        "solubility_law": "none",
        "K_S_0": None,
        "E_K_S": None,
        "thermal_conductivity": None,
        "density": None,
        "heat_capacity": None,
    }

    @staticmethod
    def build_material(festim, mat_name: str, mat_def: "MaterialDef"):
        """Construct a ``festim.Material`` from a flowsheet ``MaterialDef``."""
        extra = mat_def.extra or {}
        for field, default in FestimBuildHelpers._MATERIAL_FIELD_DEFAULTS.items():
            if field not in extra:
                logger.warning(
                    f"FestimProvider: material '{mat_name}' is missing "
                    f"'extra.{field}'; defaulting to {default!r}."
                )
        return festim.Material(
            D_0=extra["D_0"],
            E_D=extra["E_D"],
            K_S_0=extra.get(
                "K_S_0", FestimBuildHelpers._MATERIAL_FIELD_DEFAULTS["K_S_0"]
            ),
            E_K_S=extra.get(
                "E_K_S", FestimBuildHelpers._MATERIAL_FIELD_DEFAULTS["E_K_S"]
            ),
            thermal_conductivity=extra.get(
                "thermal_conductivity",
                FestimBuildHelpers._MATERIAL_FIELD_DEFAULTS["thermal_conductivity"],
            ),
            density=extra.get(
                "density", FestimBuildHelpers._MATERIAL_FIELD_DEFAULTS["density"]
            ),
            heat_capacity=extra.get(
                "heat_capacity",
                FestimBuildHelpers._MATERIAL_FIELD_DEFAULTS["heat_capacity"],
            ),
            name=mat_name,
            solubility_law=extra.get(
                "solubility_law",
                FestimBuildHelpers._MATERIAL_FIELD_DEFAULTS["solubility_law"],
            ),
        )

    @staticmethod
    def build_mesh(festim, mesh_cfg):
        """Construct a ``festim.Mesh1D`` from a :class:`MeshConfig`."""
        import numpy as np

        vertices = mesh_cfg.vertices
        if vertices is None:
            blocks = [
                np.linspace(seg.start, seg.stop, seg.num) for seg in mesh_cfg.segments
            ]
            vertices = np.concatenate(blocks)
        return festim.Mesh1D(vertices)

    @staticmethod
    def build_species(festim, model: FestimModel) -> tuple:
        """Construct ``festim.Species`` / ``festim.ImplicitSpecies`` objects.

        Returns:
            ``(species_list, all_map)`` where ``species_list`` contains only
            the real :class:`Species` objects (what
            ``HydrogenTransportProblem`` accepts) and ``all_map`` maps every
            name (including implicit species) to its object.
        """
        species_map: dict = {}
        for s in model.species:
            species_map[s.name] = festim.Species(name=s.name, mobile=s.mobile)

        implicit_map: dict = {}
        for im in model.implicit_species:
            others = [species_map[o] for o in im.others]
            implicit_map[im.name] = festim.ImplicitSpecies(
                n=im.n, others=others, name=im.name
            )

        return list(species_map.values()), {**species_map, **implicit_map}

    @staticmethod
    def build_subdomains(festim, model: FestimModel, materials_map: dict) -> tuple:
        """Construct volume + surface subdomains.

        Returns:
            ``(subdomain_list, volume_map, surface_map)``.
        """
        volume_map: dict = {}
        for v in model.subdomains.volume:
            mat = materials_map.get(v.material)
            if mat is None:
                raise ValueError(
                    f"volume subdomain id={v.id} references material "
                    f"'{v.material}', which is not defined in the flowsheet "
                    f"materials registry; available: {sorted(materials_map)}"
                )
            volume_map[v.id] = festim.VolumeSubdomain1D(
                id=v.id, borders=v.borders, material=mat
            )

        surface_map: dict = {}
        for s in model.subdomains.surface:
            surface_map[s.id] = festim.SurfaceSubdomain1D(id=s.id, x=s.x)

        subdomains = list(volume_map.values()) + list(surface_map.values())
        return subdomains, volume_map, surface_map

    @staticmethod
    def build_reactions(festim, model: FestimModel, all_map: dict, volume_map: dict):
        """Construct ``festim.Reaction`` objects."""
        reactions = []
        for r in model.reactions:
            reactant = [all_map[name] for name in r.reactant]
            kwargs = {
                "reactant": reactant,
                "k_0": r.k_0,
                "E_k": r.E_k,
                "volume": volume_map[r.volume],
            }
            if r.product:
                kwargs.update(
                    product=[all_map[name] for name in r.product],
                    p_0=r.p_0,
                    E_p=r.E_p,
                )
            reactions.append(festim.Reaction(**kwargs))
        return reactions

    @staticmethod
    def build_sources(festim, model: FestimModel, all_map: dict, volume_map: dict):
        """Construct ``festim.ParticleSource`` objects."""
        sources = []
        for src in model.sources:
            sources.append(
                festim.ParticleSource(
                    value=_build_expression_value(src.value),
                    volume=volume_map[src.volume],
                    species=all_map[src.species],
                )
            )
        return sources

    @staticmethod
    def build_initial_conditions(
        festim, model: FestimModel, all_map: dict, volume_map: dict
    ):
        """Construct ``festim.InitialConcentration`` objects."""
        ics = []
        for ic in model.initial_conditions:
            ics.append(
                festim.InitialConcentration(
                    value=ic.value,
                    volume=volume_map[ic.volume],
                    species=all_map[ic.species],
                )
            )
        return ics

    @staticmethod
    def build_boundary_conditions(
        festim, model: FestimModel, all_map: dict, surface_map: dict
    ):
        """Construct FESTIM boundary-condition objects from the schema."""
        bcs = []
        for bc in model.boundary_conditions:
            surface = surface_map[bc.subdomain]
            species = all_map[bc.species]
            if isinstance(bc, FixedConcentrationBCConfig):
                bcs.append(
                    festim.FixedConcentrationBC(
                        subdomain=surface,
                        value=_build_expression_value(bc.value),
                        species=species,
                    )
                )
            elif isinstance(bc, ParticleFluxBCConfig):
                bcs.append(
                    festim.ParticleFluxBC(
                        subdomain=surface,
                        value=_build_expression_value(bc.value),
                        species=species,
                    )
                )
            elif isinstance(bc, SievertsBCConfig):
                bcs.append(
                    festim.SievertsBC(
                        subdomain=surface,
                        S_0=bc.S_0,
                        E_S=bc.E_S,
                        pressure=bc.pressure,
                        species=species,
                    )
                )
            elif isinstance(bc, HenrysBCConfig):
                bcs.append(
                    festim.HenrysBC(
                        subdomain=surface,
                        H_0=bc.H_0,
                        E_H=bc.E_H,
                        pressure=bc.pressure,
                        species=species,
                    )
                )
        return bcs

    @staticmethod
    def build_temperature(festim, model: FestimModel):
        """Resolve the temperature into a float or a ``t -> float`` callable."""
        temp = model.temperature
        if isinstance(temp, float):
            return temp
        if temp.type == "constant":
            return temp.value
        return _build_ramp_callable(temp)

    @staticmethod
    def build_exports(
        festim, model: FestimModel, all_map: dict, volume_map: dict, surface_map: dict
    ):
        """Construct FESTIM derived-quantity / profile export objects."""
        exports = []
        for ex in model.exports:
            if isinstance(ex, SurfaceFluxExportConfig):
                exports.append(
                    festim.SurfaceFlux(
                        field=all_map[ex.field],
                        surface=surface_map[ex.surface],
                        filename=ex.filename,
                    )
                )
            elif isinstance(ex, TotalVolumeExportConfig):
                exports.append(
                    festim.TotalVolume(
                        field=all_map[ex.field],
                        volume=volume_map[ex.volume],
                        filename=ex.filename,
                    )
                )
            elif isinstance(ex, Profile1DExportConfig):
                exports.append(
                    festim.Profile1DExport(
                        field=all_map[ex.field],
                        subdomain=volume_map[ex.volume],
                        times=ex.times,
                    )
                )
        return exports

    @staticmethod
    def build_settings(festim, model: FestimModel):
        """Construct ``festim.Stepsize`` + ``festim.Settings``."""
        stepsize = None
        if model.stepsize is not None:
            ss = model.stepsize
            kwargs: dict = {"initial_value": ss.initial_value}
            if ss.growth_factor is not None:
                kwargs.update(
                    growth_factor=ss.growth_factor,
                    cutback_factor=ss.cutback_factor,
                    target_nb_iterations=ss.target_nb_iterations,
                )
            if ss.max_stepsize is not None:
                kwargs["max_stepsize"] = _build_max_stepsize(ss.max_stepsize)
            if ss.milestones is not None:
                kwargs.update(
                    milestones=ss.milestones,
                    milestone_tolerance=ss.milestone_tolerance,
                )
            stepsize = festim.Stepsize(**kwargs)

        return festim.Settings(
            atol=model.atol,
            rtol=model.rtol,
            max_iterations=model.max_iterations,
            transient=model.transient,
            final_time=model.final_time,
            element_degree=model.element_degree,
            stepsize=stepsize,
            convergence_criterion=model.convergence_criterion,
        )


# ---------------------------------------------------------------------------
# Simulation type strategy registry
# ---------------------------------------------------------------------------


class FestimSimStrategy(ABC):
    """Base class for a FESTIM simulation type.

    Each subclass encapsulates the model-setup logic for one ``sim_type`` value.
    The ``build`` method receives the lazily-imported ``festim`` module, the
    typed solver config, the material map, and the shared helpers, and returns
    the keyword arguments for ``festim.HydrogenTransportProblem(**kwargs)``.

    Register subclasses with :func:`register_festim_sim_type`.

    Example::

        class MyTDS(FestimSimStrategy):
            def build(self, festim, festim_model, materials_map, helpers):
                ...
                return problem_kwargs

        register_festim_sim_type("my_tds", MyTDS)
    """

    @abstractmethod
    def build(
        self,
        festim,
        festim_model: FestimModel,
        materials_map: dict,
        helpers: FestimBuildHelpers,
    ) -> dict:
        """Set up the FESTIM model and return ``HydrogenTransportProblem`` kwargs.

        Args:
            festim:       The ``festim`` module (imported lazily by the provider).
            festim_model: Typed solver config from the flowsheet JSON.
            materials_map: ``{mat_name: festim.Material}`` built by the provider.
            helpers:      Shared :class:`FestimBuildHelpers` instance.

        Returns:
            Keyword arguments for ``festim.HydrogenTransportProblem(...)``.
        """


_SIM_TYPE_REGISTRY: dict[str, type] = {}


def register_festim_sim_type(name: str, strategy_cls: type) -> None:
    """Register a FESTIM simulation type by name.

    Args:
        name:         The ``sim_type`` string used in the flowsheet JSON.
        strategy_cls: A subclass of :class:`FestimSimStrategy`.
    """
    _SIM_TYPE_REGISTRY[name] = strategy_cls


def get_registered_sim_types() -> dict[str, type]:
    """Return a view of the currently registered FESTIM sim_type → strategy mapping.

    Use this function (rather than accessing ``_SIM_TYPE_REGISTRY`` directly) so
    that callers are insulated from internal implementation changes.
    """
    return dict(_SIM_TYPE_REGISTRY)


# ---------------------------------------------------------------------------
# Built-in strategies
# ---------------------------------------------------------------------------


class _TDSTrappingStrategy(FestimSimStrategy):
    """1-D hydrogen transport with trapping/detrapping (TDS-type simulation).

    Builds a full ``HydrogenTransportProblem`` from the typed
    :class:`FestimModel`: refined 1-D mesh, tungsten material, mobile +
    immobile (trapped) species with implicit empty-trap species, Arrhenius
    trapping reactions, Gaussian implantation source, surface BCs, temperature
    ramp, adaptive stepsize with milestones, and derived-quantity exports.
    """

    sim_type: str = "hydrogen_transport_tds"

    def build(
        self,
        festim,
        festim_model: FestimModel,
        materials_map: dict,
        helpers: FestimBuildHelpers,
    ) -> dict:
        mesh = helpers.build_mesh(festim, festim_model.mesh)
        subdomains, volume_map, surface_map = helpers.build_subdomains(
            festim, festim_model, materials_map
        )
        species, all_map = helpers.build_species(festim, festim_model)
        reactions = helpers.build_reactions(festim, festim_model, all_map, volume_map)
        sources = helpers.build_sources(festim, festim_model, all_map, volume_map)
        initial_conditions = helpers.build_initial_conditions(
            festim, festim_model, all_map, volume_map
        )
        boundary_conditions = helpers.build_boundary_conditions(
            festim, festim_model, all_map, surface_map
        )
        temperature = helpers.build_temperature(festim, festim_model)
        exports = helpers.build_exports(
            festim, festim_model, all_map, volume_map, surface_map
        )
        settings = helpers.build_settings(festim, festim_model)

        return {
            "mesh": mesh,
            "subdomains": subdomains,
            "species": species,
            "reactions": reactions,
            "sources": sources,
            "initial_conditions": initial_conditions,
            "boundary_conditions": boundary_conditions,
            "temperature": temperature,
            "exports": exports,
            "settings": settings,
        }


# Register built-ins at module load
register_festim_sim_type("hydrogen_transport_tds", _TDSTrappingStrategy)


# ---------------------------------------------------------------------------
# FestimProvider
# ---------------------------------------------------------------------------


class FestimProvider(AbstractProvider):
    """In-process hydrogen transport provider backed by FESTIM.

    Responsibilities
    ----------------
    * **Installation check** — ``initialize()`` verifies ``festim`` is
      importable before the first run.
    * **Material registry** — ``initialize()`` loads all material defs from the
      flowsheet so strategies can build ``festim.Material`` objects.
    * **Working-directory management** — each run executes in
      ``<provider_output_dir>`` (``os.chdir`` around ``initialise()``/``run()``
      so FESTIM's CSV exports land there). Falls back to a temp dir when the
      configured output dir is not writable.
    * **Simulation dispatch** — ``run_simulation()`` parses ``solver_config``
      into a :class:`FestimModel`, looks up the strategy from
      ``_SIM_TYPE_REGISTRY``, calls ``strategy.build()``, then runs
      ``HydrogenTransportProblem.initialise()`` + ``.run()``. Failures are
      returned as ``SimulationResult(status="failed")`` with ``run_dir`` and
      error metadata rather than silently bubbling up.
    * **Result extraction** — per-export scalar values (surface fluxes, volume
      inventories) become ``scalars``; full time-series and profile data are
      captured in ``metadata["series"]``.
    """

    def __init__(self):
        self._materials: dict = {}
        self._provider_output_dir: str = "outputs/festim"
        self._initialized: bool = False

    def initialize(
        self,
        provider_config: "FestimProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        """Verify FESTIM is installed, store config, build material registry."""
        try:
            import festim  # noqa: F401
        except ImportError as exc:
            raise RuntimeError(
                "FESTIM is not installed. Install with: pip install festim"
            ) from exc

        out_dir = os.path.expandvars(provider_config.output_dir)
        if not os.path.isabs(out_dir):
            root = os.environ.get("PROCESSFORGE_OUTPUT_DIR", "outputs")
            out_dir = os.path.join(root, out_dir)
        self._provider_output_dir = out_dir

        self._materials = dict(flowsheet_config.materials.items())

        self._initialized = True
        logger.info(
            f"FestimProvider initialized with {len(self._materials)} material(s). "
            f"Registered sim_types: {sorted(_SIM_TYPE_REGISTRY)}"
        )

    def get_thermo_properties(self, stream: dict) -> dict:
        raise NotImplementedError(
            "FestimProvider does not support stream thermodynamics."
        )

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        """Return ``None`` — FESTIM uses ``run_simulation`` via ``SolverUnit``."""
        return None

    def teardown(self) -> None:
        """Release provider state."""
        self._initialized = False

    @classmethod
    def validate_material(cls, mat_name: str, mat_def, unit_cfg) -> list:
        """Validate FESTIM-specific material properties.

        Rules
        -----
        * ``extra.D_0`` is required and must be positive.
        * ``extra.E_D`` is required and must be non-negative.
        * ``extra.K_S_0`` and ``extra.E_K_S`` must be provided together.

        Returns:
            List of error strings (empty = valid).
        """
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
        """Build and run a FESTIM simulation from a typed ``UnitConfig``.

        Execution flow
        --------------
        1. Parse ``solver_config`` → ``FestimModel``
        2. Look up strategy from ``_SIM_TYPE_REGISTRY`` by ``sim_type``
        3. Build ``festim.Material`` objects for all registry materials
        4. Call ``strategy.build()`` → ``HydrogenTransportProblem`` kwargs
        5. Resolve the run directory (with temp-dir fallback)
        6. ``os.chdir(run_dir)``, then ``initialise()`` + ``run()``
        7. Extract per-export scalars and time-series (reduced + artifacts)
        8. Return :class:`EngineOutput`; failures surface as ``status="failed"``
           with ``run_dir`` and error diagnostics.
        """
        if not self._initialized:
            raise RuntimeError(
                "FestimProvider has not been initialized — call initialize() first."
            )

        sim_type = unit_config.sim_type
        strategy_cls = _SIM_TYPE_REGISTRY.get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"FestimProvider: unknown sim_type '{sim_type}'. "
                f"Register with register_festim_sim_type(). "
                f"Built-in types: {sorted(_SIM_TYPE_REGISTRY)}"
            )

        import festim
        from processforge.types import (
            EngineOutput,
            OutputArtifact,
            OutputField,
            OutputProvenance,
            Quantity,
        )

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
            return EngineOutput(
                status="failed",
                engine="festim",
                sim_type=sim_type,
                diagnostics={
                    "run_dir": str(run_dir.resolve()),
                    "error": str(exc),
                },
            )
        finally:
            os.chdir(prev_cwd)

        fields, artifacts, diagnostics = self._extract_results(
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
        )

    def _resolve_run_dir(self) -> pathlib.Path:
        """Create the run output directory, falling back to a temp dir on failure.

        The configured output dir (often the mounted /data volume) may not be
        writable by the current user. Fall back to a temp dir so the simulation
        can still run; only the persisted CSV artifacts are lost.
        """
        run_dir = pathlib.Path(self._provider_output_dir)
        try:
            run_dir.mkdir(parents=True, exist_ok=True)
        except (PermissionError, OSError) as exc:
            run_dir = pathlib.Path(tempfile.mkdtemp(prefix="processforge_festim_"))
            logger.warning(
                f"FestimProvider: cannot use output dir "
                f"'{self._provider_output_dir}' ({exc}); "
                f"falling back to '{run_dir}'."
            )
        return run_dir

    @staticmethod
    def _extract_results(exports, run_dir: pathlib.Path) -> tuple:
        """Return ``(fields, artifacts, diagnostics)`` from FESTIM exports.

        * Every export exposing a ``title`` + ``value`` (surface fluxes, total
          volumes) contributes a scalar :class:`OutputField`.
        * Every export exposing ``t`` + ``data`` (incl. ``Profile1DExport``)
          contributes a reduced ``{title}_mean_total`` / ``{title}_std_dev``
          :class:`OutputField` (mirroring the OpenMC tally convention) plus the
          full series written to a content-addressed CSV :class:`OutputArtifact`.
          The field carries a ``timeseries`` ``Domain`` so downstream code can
          reconstruct the shape without the file.

        ``diagnostics`` always carries ``run_dir``.
        """
        import csv

        from processforge.types import Axis, Domain, OutputArtifact, OutputField, Quantity

        fields: list = []
        artifacts: list = []
        diagnostics: dict = {"run_dir": str(run_dir.resolve())}

        for export in exports:
            title = getattr(export, "title", None)
            value = getattr(export, "value", None)
            if title is not None and value is not None:
                fields.append(OutputField(
                    name=str(title),
                    quantity=Quantity(value=[float(value)], unit=""),
                    kind="scalar",
                    source=str(title),
                ))

            t = getattr(export, "t", None)
            data = getattr(export, "data", None)
            if not (t and data):
                continue

            key = str(title if title is not None else getattr(export, "field", ""))
            if not key:
                continue

            # Flatten all time steps into a single value array for statistics.
            flat: list[float] = []
            rows: list[list[float]] = []
            for step in data:
                if hasattr(step, "tolist"):
                    arr = step.tolist()
                else:
                    arr = [float(step)]
                if not isinstance(arr, list):
                    arr = [arr]
                rows.append(arr)
                flat.extend(arr)

            import numpy as np

            flat_arr = np.asarray(flat, dtype=float)
            mean_total = float(flat_arr.mean()) if flat_arr.size else 0.0
            std_dev = float(flat_arr.std()) if flat_arr.size else 0.0

            domain = Domain(
                type="timeseries",
                axes=[Axis(name="t", unit="s", coordinates=list(map(float, t)))],
                shape=[len(rows), max((len(r) for r in rows), default=0)],
            )
            fields.append(OutputField(
                name=f"{key}_mean_total",
                quantity=Quantity(value=mean_total, unit=""),
                kind="timeseries",
                source=f"{key}/mean",
                domain=domain,
            ))
            fields.append(OutputField(
                name=f"{key}_std_dev",
                quantity=Quantity(value=std_dev, unit=""),
                kind="timeseries",
                source=f"{key}/std_dev",
                domain=domain,
            ))

            # Persist the full series as a CSV artifact (not inlined over HTTP).
            csv_name = f"{key}.csv"
            csv_path = pathlib.Path(run_dir) / csv_name
            try:
                with open(csv_path, "w", newline="", encoding="utf-8") as fh:
                    writer = csv.writer(fh)
                    writer.writerow(["t"] + [f"v{i}" for i in range(max((len(r) for r in rows), default=0))])
                    for ti, row in zip(t, rows):
                        writer.writerow([float(ti)] + list(row))
                artifacts.append(OutputArtifact(
                    name=csv_name,
                    kind="csv",
                    local_path=str(csv_path.resolve()),
                    source="local",
                ))
            except Exception as exc:  # noqa: BLE001
                logger.warning(f"FestimProvider: failed to write series CSV '{csv_name}': {exc}")

        return fields, artifacts, diagnostics


# Register the provider
register_provider("festim", FestimProvider)
