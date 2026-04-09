"""FESTIM provider for hydrogen transport and heat transfer simulations.

Architecture
------------
This module wires FESTIM into the processforge provider system using three layers:

1. **Dataclasses** (``FestimBoundaryCondition``, ``FestimSolverConfig``, etc.)
   Parse the opaque ``solver_config`` JSON dict into typed, self-documenting objects.

2. **Strategy registry** (``FestimSimStrategy`` + ``register_festim_sim_type``)
   Each ``sim_type`` string maps to a strategy class that builds a ``F.Simulation``.
   New simulation types are added by subclassing ``FestimSimStrategy`` and calling
   ``register_festim_sim_type(name, MyStrategy)`` — no changes to this file needed.

3. **FestimProvider** (``AbstractProvider`` subclass)
   Owns material lookup (``initialize``), material validation (``validate_material``),
   and simulation dispatch (``run_simulation``).  All mesh/BC/trap building lives in
   ``FestimBuildHelpers`` so strategy classes are concise and testable in isolation.

Adding a new sim_type example::

    from processforge.providers.festim_provider import (
        FestimSimStrategy, FestimSolverConfig, FestimBuildHelpers,
        register_festim_sim_type,
    )

    class MyCustomSim(FestimSimStrategy):
        def build(self, model, F, solver_cfg, mat_def, helpers):
            # configure model ...
            return {}

    register_festim_sim_type("my_custom_3d", MyCustomSim)
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional, Union

from loguru import logger

from .base import AbstractProvider
from .registry import register_provider


# ---------------------------------------------------------------------------
# Festim-specific config dataclasses
# ---------------------------------------------------------------------------

@dataclass
class FestimBoundaryCondition:
    """One boundary condition entry from ``solver_config["boundary_conditions"]``.

    type:     "DirichletBC" | "ConvectiveFlux" | "FluxBC" | "SievertsBC"
    surfaces: FESTIM surface marker integer (1-based, matches mesh marker values)
    field:    "T" for temperature BCs; 0 (int) for hydrogen solute (FESTIM convention)
    value:    constant float or a FEniCS/sympy-parseable expression string
    h_coeff:  convective heat transfer coefficient — ``ConvectiveFlux`` only
    T_ext:    external temperature — ``ConvectiveFlux`` only
    S_0, E_S, pressure: Sieverts law parameters — ``SievertsBC`` only
    """
    type: str
    surfaces: int
    value: Optional[Union[float, str]] = None
    field: Optional[Union[str, int]] = "T"
    h_coeff: Optional[Union[float, str]] = None
    T_ext: Optional[Union[float, str]] = None
    S_0: Optional[float] = None
    E_S: Optional[float] = None
    pressure: Optional[float] = None

    @classmethod
    def from_dict(cls, d: dict) -> "FestimBoundaryCondition":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class FestimSource:
    """A volumetric source term (heat generation, particle source, etc.)."""
    value: Union[float, str]       # constant or expression string
    volume: int                    # volume marker integer
    field: Union[str, int] = "T"   # "T" for heat, 0 for hydrogen solute


@dataclass
class FestimTrap:
    """A hydrogen trap site definition.

    See FESTIM ``F.Trap`` for full semantics.

    k_0:     pre-exponential trapping rate [m³/s]
    E_k:     trapping activation energy [eV]
    p_0:     pre-exponential detrapping rate [1/s]
    E_p:     detrapping activation energy [eV]
    density: trap site density [traps/m³]
    """
    k_0: float
    E_k: float
    p_0: float
    E_p: float
    density: float


@dataclass
class FestimExport:
    """One XDMF field export entry."""
    field: str       # "solute" | "retention" | "T"
    folder: str = "output"


@dataclass
class FestimMesh1D:
    """1D mesh config: a uniform line from 0 to ``length``."""
    length: float = 1.0
    num_vertices: int = 1000


@dataclass
class FestimMesh2D:
    """2D mesh config: a unit square with ``nx`` × ``ny`` cells.

    Surface markers follow the simple_2D.py convention:
      left  (x=0) → 1,  top+bottom → 2,  right (x=1) → 3
    """
    nx: int = 10
    ny: int = 10


@dataclass
class FestimSolverConfig:
    """Typed representation of the ``solver_config`` block for Festim ``SolverUnit`` units.

    Parsed from the opaque ``solver_config`` JSON dict by
    ``FestimProvider.run_simulation()``.  All fields map 1-to-1 to FESTIM API concepts.

    Example JSON::

        "solver_config": {
            "mesh": {"nx": 20, "ny": 20},
            "transient": false,
            "boundary_conditions": [
                {"type": "DirichletBC", "field": "T", "value": 300, "surfaces": 1}
            ],
            "traps": [{"k_0": 3.8e-17, "E_k": 0.39, "p_0": 8.4e12, "E_p": 0.9, "density": 1e25}],
            "exports": [{"field": "solute"}, {"field": "T"}],
            "export_folder": "output/my_sim"
        }
    """
    boundary_conditions: list = field(default_factory=list)   # list[FestimBoundaryCondition]
    sources: list = field(default_factory=list)               # list[FestimSource]
    traps: list = field(default_factory=list)                 # list[FestimTrap]
    exports: list = field(default_factory=list)               # list[FestimExport]
    mesh: Optional[Union[FestimMesh1D, FestimMesh2D]] = None
    transient: bool = False
    T_fixed: Optional[float] = None
    export_folder: str = "output"

    @classmethod
    def from_dict(cls, d: dict) -> "FestimSolverConfig":
        """Parse a raw ``solver_config`` dict into typed Festim dataclasses."""
        bcs = [FestimBoundaryCondition.from_dict(bc)
               for bc in d.get("boundary_conditions", [])]
        sources = [FestimSource(**s) for s in d.get("sources", [])]
        traps = [FestimTrap(**t) for t in d.get("traps", [])]
        exports = [FestimExport(**e) for e in d.get("exports", [])]

        mesh_raw = d.get("mesh", {})
        if "nx" in mesh_raw:
            mesh: Optional[Union[FestimMesh1D, FestimMesh2D]] = FestimMesh2D(**mesh_raw)
        elif mesh_raw:
            mesh = FestimMesh1D(**mesh_raw)
        else:
            mesh = None

        return cls(
            boundary_conditions=bcs,
            sources=sources,
            traps=traps,
            exports=exports,
            mesh=mesh,
            transient=d.get("transient", False),
            T_fixed=d.get("T_fixed"),
            export_folder=d.get("export_folder", "output"),
        )


# ---------------------------------------------------------------------------
# Build helpers — shared mesh/material/BC builders used by all strategies
# ---------------------------------------------------------------------------

class FestimBuildHelpers:
    """Shared FESTIM model-building utilities injected into sim strategies.

    Separating these from the provider keeps strategy classes concise and
    allows strategies to be unit-tested without instantiating a full provider.
    """

    @staticmethod
    def build_material(F, mat_def) -> object:
        """Construct ``F.Material`` from a :class:`~processforge.types.MaterialDef`.

        ``thermal_cond``:
          - ``float``: passed directly
          - ``str``: parsed via ``sympy.lambdify`` (falls back to restricted ``eval``
            if sympy is unavailable).  The expression must be valid Python in ``T``,
            e.g. ``"3 + 0.1 * T"``.
        """
        tc_raw = mat_def.extra.get("thermal_cond")
        if isinstance(tc_raw, str):
            try:
                import sympy
                T_sym = sympy.Symbol("T")
                expr = sympy.sympify(tc_raw)
                tc = sympy.lambdify(T_sym, expr, modules="numpy")
            except ImportError:
                # sympy not available — fall back to restricted eval
                tc = lambda T, _expr=tc_raw: eval(_expr, {"T": T, "__builtins__": {}})  # noqa: E731
        else:
            tc = tc_raw  # float constant or None

        kwargs = {
            "id": 1,
            "D_0": mat_def.extra.get("D_0"),
            "E_D": mat_def.extra.get("E_D"),
        }
        if tc is not None:
            kwargs["thermal_cond"] = tc
        # Filter out None values — FESTIM raises on unexpected None kwargs
        return F.Material(**{k: v for k, v in kwargs.items() if v is not None})

    @staticmethod
    def build_mesh_1d(F, mesh: FestimMesh1D) -> object:
        """Build a 1D mesh from equally-spaced vertices."""
        import numpy as np
        return F.MeshFromVertices(np.linspace(0, mesh.length, mesh.num_vertices))

    @staticmethod
    def build_mesh_2d(mesh: FestimMesh2D):
        """Build a FEniCS ``UnitSquareMesh`` with volume and surface markers.

        Surface convention (matches ``simple_2D.py``):
          left  (x=0)    → marker 1
          top + bottom   → marker 2
          right (x=1)    → marker 3

        Returns:
            tuple: ``(fenics_mesh, volume_markers, surface_markers)``
        """
        from fenics import CompiledSubDomain, MeshFunction, UnitSquareMesh
        m = UnitSquareMesh(mesh.nx, mesh.ny)

        vol = MeshFunction("size_t", m, m.topology().dim())
        vol.set_all(1)

        surf = MeshFunction("size_t", m, m.topology().dim() - 1)
        surf.set_all(0)
        CompiledSubDomain("on_boundary && near(x[0], 0, 1e-14)").mark(surf, 1)
        CompiledSubDomain("on_boundary && near(x[0], 1, 1e-14)").mark(surf, 3)
        CompiledSubDomain(
            "on_boundary && (near(x[1], 0, 1e-14) || near(x[1], 1, 1e-14))"
        ).mark(surf, 2)
        return m, vol, surf

    @staticmethod
    def build_bcs(F, bcs: list) -> list:
        """Convert a list of :class:`FestimBoundaryCondition` to FESTIM BC objects."""
        result = []
        for bc in bcs:
            if bc.type == "DirichletBC":
                result.append(
                    F.DirichletBC(value=bc.value, field=bc.field, surfaces=bc.surfaces)
                )
            elif bc.type == "ConvectiveFlux":
                result.append(
                    F.ConvectiveFlux(h_coeff=bc.h_coeff, T_ext=bc.T_ext, surfaces=bc.surfaces)
                )
            elif bc.type == "FluxBC":
                result.append(
                    F.FluxBC(value=bc.value, field=bc.field, surfaces=bc.surfaces)
                )
            elif bc.type == "SievertsBC":
                result.append(
                    F.SievertsBC(
                        surfaces=bc.surfaces, S_0=bc.S_0, E_S=bc.E_S, pressure=bc.pressure
                    )
                )
            else:
                raise ValueError(
                    f"Unknown Festim BC type '{bc.type}'. "
                    "Valid: DirichletBC, ConvectiveFlux, FluxBC, SievertsBC"
                )
        return result

    @staticmethod
    def build_sources(F, sources: list) -> list:
        """Convert a list of :class:`FestimSource` to ``F.Source`` objects."""
        return [F.Source(value=s.value, volume=s.volume, field=s.field) for s in sources]

    @staticmethod
    def build_traps(F, traps: list, mat: object) -> object:
        """Build a ``F.Traps`` collection from :class:`FestimTrap` dataclasses."""
        return F.Traps([
            F.Trap(
                k_0=t.k_0, E_k=t.E_k, p_0=t.p_0, E_p=t.E_p,
                density=t.density, materials=mat,
            )
            for t in traps
        ])

    @staticmethod
    def build_exports(F, exports: list, folder: str) -> object:
        """Build a ``F.Exports`` collection from :class:`FestimExport` dataclasses."""
        return F.Exports([F.XDMFExport(e.field, folder=folder) for e in exports])


# ---------------------------------------------------------------------------
# Simulation type strategy registry
# ---------------------------------------------------------------------------

class FestimSimStrategy(ABC):
    """Base class for a Festim simulation type.

    Each subclass encapsulates the ``F.Simulation`` configuration logic for
    one ``sim_type`` value.  Register subclasses with
    :func:`register_festim_sim_type` — the provider's ``run_simulation``
    dispatcher never needs to change to support new types.

    Example::

        class MyCustomSim(FestimSimStrategy):
            def build(self, model, F, solver_cfg, mat_def, helpers):
                model.mesh = helpers.build_mesh_1d(F, solver_cfg.mesh)
                ...
                return {}   # scalar results (populated after model.run())

        register_festim_sim_type("my_custom_3d", MyCustomSim)
    """

    @abstractmethod
    def build(
        self,
        model: object,
        F: object,
        solver_cfg: FestimSolverConfig,
        mat_def: object,
        helpers: FestimBuildHelpers,
    ) -> dict:
        """Configure ``model`` in-place.

        Called *before* ``model.initialise()`` and ``model.run()``.

        Args:
            model:      A blank ``F.Simulation()`` instance to configure.
            F:          The ``festim`` module (avoids redundant re-imports inside strategies).
            solver_cfg: Typed solver config from the flowsheet unit's ``solver_config`` dict.
            mat_def:    :class:`~processforge.types.MaterialDef` for the structural material.
            helpers:    Shared :class:`FestimBuildHelpers` (mesh/BC/trap builders).

        Returns:
            Dict of named scalar results for ``SimulationResult.scalars``.
            Return ``{}`` when no scalars apply.
        """


# Module-level registry: sim_type string → strategy class
_SIM_TYPE_REGISTRY: dict = {}


def register_festim_sim_type(name: str, strategy_cls: type) -> None:
    """Register a Festim simulation type by name.

    Args:
        name:         The ``sim_type`` string used in the flowsheet JSON.
        strategy_cls: A subclass of :class:`FestimSimStrategy`.
    """
    _SIM_TYPE_REGISTRY[name] = strategy_cls


# ---------------------------------------------------------------------------
# Built-in strategies
# ---------------------------------------------------------------------------

class _Heat1DStrategy(FestimSimStrategy):
    """1D steady or transient heat transfer.

    Corresponds to ``simple_1D.py``:
    MeshFromVertices → Material(thermal_cond) → HeatTransferProblem → BCs + sources.
    """

    def build(self, model, F, solver_cfg, mat_def, helpers):
        model.mesh = helpers.build_mesh_1d(F, solver_cfg.mesh or FestimMesh1D())
        model.materials = helpers.build_material(F, mat_def)
        model.T = F.HeatTransferProblem(transient=solver_cfg.transient)
        model.boundary_conditions = helpers.build_bcs(F, solver_cfg.boundary_conditions)
        model.sources = helpers.build_sources(F, solver_cfg.sources)
        model.settings = F.Settings(
            transient=solver_cfg.transient,
            absolute_tolerance=1e-10,
            relative_tolerance=1e-10,
        )
        return {}


class _Heat2DStrategy(FestimSimStrategy):
    """2D coupled heat transfer + hydrogen transport.

    Corresponds to ``simple_2D.py``:
    UnitSquareMesh → Material(thermal_cond, D_0, E_D) → traps →
    HeatTransferProblem → BCs + sources → XDMF exports.

    Surface marker convention: left=1, top+bottom=2, right=3.
    """

    def build(self, model, F, solver_cfg, mat_def, helpers):
        mesh_cfg = solver_cfg.mesh or FestimMesh2D()
        fenics_mesh, vol_markers, surf_markers = helpers.build_mesh_2d(mesh_cfg)
        mat = helpers.build_material(F, mat_def)

        model.materials = F.Materials([mat])
        model.mesh = F.Mesh(
            mesh=fenics_mesh,
            volume_markers=vol_markers,
            surface_markers=surf_markers,
        )

        if solver_cfg.traps:
            model.traps = helpers.build_traps(F, solver_cfg.traps, mat)

        model.T = F.HeatTransferProblem(transient=solver_cfg.transient)
        model.boundary_conditions = helpers.build_bcs(F, solver_cfg.boundary_conditions)
        model.sources = helpers.build_sources(F, solver_cfg.sources)

        if solver_cfg.exports:
            model.exports = helpers.build_exports(
                F, solver_cfg.exports, solver_cfg.export_folder
            )

        model.settings = F.Settings(
            transient=solver_cfg.transient,
            absolute_tolerance=1e-9,
            relative_tolerance=1e-9,
        )
        return {}


class _Permeation1DStrategy(FestimSimStrategy):
    """1D transient hydrogen permeation.

    Corresponds to ``permeation.py``:
    MeshFromVertices → Material(D_0, E_D) → T=constant →
    SievertsBC (upstream) + DirichletBC=0 (downstream) → HydrogenFlux export.
    """

    def build(self, model, F, solver_cfg, mat_def, helpers):
        model.mesh = helpers.build_mesh_1d(
            F, solver_cfg.mesh or FestimMesh1D(length=3e-4, num_vertices=1001)
        )
        model.materials = helpers.build_material(F, mat_def)
        model.T = solver_cfg.T_fixed if solver_cfg.T_fixed is not None else 500.0
        model.boundary_conditions = helpers.build_bcs(F, solver_cfg.boundary_conditions)

        derived = F.DerivedQuantities([F.HydrogenFlux(surface=2)], show_units=True)
        model.exports = [derived]

        model.settings = F.Settings(
            absolute_tolerance=1e-2,
            relative_tolerance=1e-10,
            final_time=100.0,
        )
        return {}


# Register built-ins at module load
register_festim_sim_type("heat_1d", _Heat1DStrategy)
register_festim_sim_type("heat_2d", _Heat2DStrategy)
register_festim_sim_type("permeation_1d", _Permeation1DStrategy)


# ---------------------------------------------------------------------------
# Material validation constants
# ---------------------------------------------------------------------------

_TRANSPORT_REQUIRED = frozenset({"D_0", "E_D"})
_PERMEATION_REQUIRED = frozenset({"D_0", "E_D", "S_0", "E_S"})

_FESTIM_MATERIAL_KEYS = frozenset({
    "D_0", "E_D", "S_0", "E_S", "thermal_cond",
    "k_0", "E_k", "p_0", "E_p", "friendly_material_id",
})


# ---------------------------------------------------------------------------
# FestimProvider
# ---------------------------------------------------------------------------

class FestimProvider(AbstractProvider):
    """Hydrogen transport and heat transfer provider backed by FESTIM.

    Responsibilities
    ----------------
    * **Material registry** — ``initialize()`` loads material defs from the
      flowsheet and indexes them by ``friendly_material_id`` and name.
    * **Material validation** — ``validate_material()`` enforces Festim-specific
      required properties (D_0/E_D, and S_0/E_S for permeation).
    * **Simulation dispatch** — ``run_simulation()`` looks up the strategy for
      the requested ``sim_type`` from ``_SIM_TYPE_REGISTRY`` and runs it.
    * **Stream-unit fallback** — ``compute_unit()`` returns ``None`` so
      ``FestimMembrane`` units run their own FEM stepping.
    """

    def __init__(self):
        self._materials: dict = {}
        self._initialized: bool = False

    def initialize(self, provider_config: dict, flowsheet_config: dict) -> None:
        """Initialize the FESTIM provider.

        Checks that FESTIM is installed, then builds an internal
        :class:`~processforge.types.MaterialDef` registry from the flowsheet's
        global ``materials`` section.

        Only Festim-relevant fields (D_0, E_D, S_0, E_S, thermal_cond, trap
        parameters) are retained — OpenMC/other fields (nuclides, density_units)
        are excluded, keeping the provider footprint minimal.

        Provider-block ``"materials"`` overrides take priority over the global
        section.  Materials are indexed by both ``friendly_material_id`` (as str)
        and name for flexible lookup.
        """
        try:
            import festim  # noqa: F401
        except ImportError as exc:
            raise RuntimeError(
                "FESTIM is not installed. Install with: pip install festim"
            ) from exc

        from processforge.types import MaterialDef

        self._materials = {}
        for mat_name, mat_dict in flowsheet_config.get("materials", {}).items():
            filtered = {k: v for k, v in mat_dict.items() if k in _FESTIM_MATERIAL_KEYS}
            mat = MaterialDef.from_dict(filtered)
            fid = mat_dict.get("friendly_material_id")
            if fid is not None:
                self._materials[str(fid)] = mat
            self._materials[mat_name] = mat

        # Provider-block inline overrides take priority
        for key, props in provider_config.get("materials", {}).items():
            self._materials[str(key)] = MaterialDef.from_dict(props)

        self._initialized = True
        logger.info(
            f"FestimProvider initialized with {len(set(self._materials.values()))} material(s). "
            f"Registered sim_types: {sorted(_SIM_TYPE_REGISTRY)}"
        )

    def get_thermo_properties(self, stream: dict) -> dict:
        """Fall back to CoolProp for stream thermodynamics.

        FESTIM handles hydrogen transport; enthalpy, Cp, and K-values
        are delegated to CoolProp.
        """
        from processforge.thermo import get_Cp_molar, get_enthalpy_molar, get_K_values

        z, T, P = stream["z"], stream["T"], stream["P"]
        return {
            "H": get_enthalpy_molar(z, T, P),
            "Cp": get_Cp_molar(z, T, P),
            "K_values": get_K_values(list(z.keys()), T, P),
        }

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        """Return ``None`` for all FESTIM units.

        ``FestimMembrane`` drives its own FEM time-stepping in ``_run_impl`` /
        ``run_dynamic``; returning ``None`` lets the unit's own logic run.
        ``SolverUnit`` never calls ``compute_unit`` — it calls ``run_simulation``
        directly.
        """
        return None

    def teardown(self) -> None:
        """Release provider-level resources."""
        self._initialized = False

    @classmethod
    def validate_material(cls, mat_name: str, mat_def, unit_cfg) -> list:
        """Validate Festim-specific material properties.

        Rules
        -----
        * All Festim units:       ``D_0``, ``E_D`` required
        * Permeation units:       also ``S_0``, ``E_S`` required
          (FestimMembrane or sim_type="permeation_1d")
        * ``thermal_cond`` string: must be syntactically valid Python in ``T``

        Returns:
            List of error strings (empty = valid).
        """
        errors = []
        needs_permeation = (
            unit_cfg.type == "FestimMembrane"
            or unit_cfg.sim_type == "permeation_1d"
        )
        required = _PERMEATION_REQUIRED if needs_permeation else _TRANSPORT_REQUIRED

        missing = required - set(mat_def.extra.keys())
        if missing:
            errors.append(
                f"Material '{mat_name}' is missing required Festim properties: "
                f"{sorted(missing)}"
            )

        tc = mat_def.extra.get("thermal_cond")
        if isinstance(tc, str):
            try:
                compile(tc, "<thermal_cond>", "eval")
            except SyntaxError as e:
                errors.append(
                    f"Material '{mat_name}' has invalid thermal_cond expression "
                    f"'{tc}': {e}"
                )
        return errors

    def run_simulation(self, unit_config, inlet: dict):
        """Build and run a FESTIM simulation from a typed ``UnitConfig``.

        Looks up the strategy class from ``_SIM_TYPE_REGISTRY`` by ``sim_type``,
        parses ``solver_config`` into a :class:`FestimSolverConfig`, retrieves
        the material, then runs ``strategy.build(model, …)`` followed by
        ``model.initialise()`` and ``model.run()``.

        To add a new simulation type, call :func:`register_festim_sim_type` —
        this method never needs to change.
        """
        import festim as F
        from processforge.types import MaterialDef, SimulationResult

        sim_type = unit_config.sim_type
        strategy_cls = _SIM_TYPE_REGISTRY.get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"FestimProvider: unknown sim_type '{sim_type}'. "
                f"Register custom types with register_festim_sim_type(). "
                f"Built-in types: {sorted(_SIM_TYPE_REGISTRY)}"
            )

        solver_cfg = FestimSolverConfig.from_dict(unit_config.solver_config)
        mat_def = self._materials.get(
            str(unit_config.material),
            MaterialDef(friendly_material_id=-1),
        )

        model = F.Simulation()
        helpers = FestimBuildHelpers()
        logger.info(
            f"FestimProvider: running sim_type='{sim_type}' "
            f"with strategy {strategy_cls.__name__}"
        )
        scalars = strategy_cls().build(model, F, solver_cfg, mat_def, helpers)
        model.initialise()
        model.run()
        logger.info(f"FestimProvider: '{sim_type}' completed")

        return SimulationResult(
            status="completed",
            sim_type=sim_type,
            scalars=scalars or {},
        )

    def get_material_props(self, material_id: str) -> dict:
        """Return the ``extra`` dict for a material by id or name.

        Used by :class:`~processforge.units.festim_membrane.FestimMembrane`
        to retrieve D_0, E_D, S_0, E_S for its internal FESTIM model.

        Returns an empty dict and logs a warning if the material is not found.
        """
        mat = self._materials.get(str(material_id)) or self._materials.get(material_id)
        if mat is None:
            logger.warning(
                f"FestimProvider: material '{material_id}' not found. "
                "FestimMembrane will use its own defaults."
            )
            return {}
        return mat.extra


# Register the provider
register_provider("festim", FestimProvider)
