"""FESTIM provider for hydrogen transport finite element simulations.

Architecture (follows OpenMC pattern)
--------------------------------------
1. **Build helpers** — module-level functions that translate JSON config dicts
   into FESTIM objects by calling FESTIM constructors directly.

2. **Strategy registry** — each ``sim_type`` string maps to a strategy class
   that builds a FESTIM problem.  New simulation types are added by
   subclassing ``FestimSimStrategy`` and calling
   ``register_festim_sim_type(name, MyStrategy)``.

3. **FestimProvider** — ``AbstractProvider`` subclass.  Owns material
   registry, validation, and simulation dispatch.

Adding a new sim_type::

    from processforge.providers.festim_provider import (
        FestimSimStrategy, register_festim_sim_type,
    )

    class MyCustomSim(FestimSimStrategy):
        def build(self, F, unit_cfg, materials_map, helpers):
            ...
            return problem_instance

    register_festim_sim_type("my_custom_sim", MyCustomSim)
"""

from __future__ import annotations

import pathlib
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional

import numpy as np
from loguru import logger

from .base import AbstractProvider
from .registry import register_provider

if TYPE_CHECKING:
    from processforge.types import (
        FestimProviderConfig,
        FlowsheetConfig,
        SimulationResult,
        UnitConfig,
    )


# ---------------------------------------------------------------------------
# Build helpers — thin translation from JSON dicts to FESTIM constructors
# ---------------------------------------------------------------------------


def _build_mesh(F, mesh_cfg: dict):
    """Build ``F.Mesh1D`` from a mesh config dict with ``vertices`` array."""
    return F.Mesh1D(vertices=np.array(mesh_cfg["vertices"]))


def _build_material(F, mat_name: str, mat_def):
    """Build ``F.Material`` from a ``MaterialDef``'s ``extra`` dict."""
    extra = mat_def.extra or {}
    return F.Material(
        D_0=extra["D_0"],
        E_D=extra["E_D"],
        name=mat_name,
        K_S_0=extra.get("K_S_0"),
        E_K_S=extra.get("E_K_S"),
        solubility_law=extra.get("solubility_law", "none"),
    )


def _build_subdomains(F, sub_cfg: dict, materials_map: dict):
    """Build volume and surface subdomain lists from config.

    Returns:
        Tuple of ``(volumes, surfaces)`` — lists of FESTIM subdomain objects.
    """
    volumes = []
    for v in sub_cfg.get("volume", []):
        mat = materials_map[v["material"]]
        volumes.append(
            F.VolumeSubdomain1D(id=v["id"], borders=v["borders"], material=mat)
        )

    surfaces = []
    for s in sub_cfg.get("surface", []):
        surfaces.append(F.SurfaceSubdomain1D(id=s["id"], x=s["x"]))

    return volumes, surfaces


def _build_species(F, species_cfgs: list):
    """Build ``F.Species`` objects from config dicts."""
    return [F.Species(name=s["name"], mobile=s.get("mobile", True)) for s in species_cfgs]


def _build_bcs(F, bc_cfgs: list, species_map: dict, surfaces_map: dict):
    """Build boundary condition objects from config dicts."""
    bcs = []
    for cfg in bc_cfgs:
        sub = surfaces_map[cfg["subdomain_id"]]
        sp = species_map[cfg["species"]]
        bc_type = cfg["type"]

        if bc_type == "fixed_concentration":
            bcs.append(
                F.FixedConcentrationBC(subdomain=sub, value=cfg["value"], species=sp)
            )
        elif bc_type == "particle_flux":
            bcs.append(
                F.ParticleFluxBC(subdomain=sub, value=cfg["value"], species=sp)
            )
        elif bc_type == "sieverts":
            bcs.append(
                F.SievertsBC(
                    subdomain=sub,
                    S_0=cfg["S_0"],
                    E_S=cfg["E_S"],
                    pressure=cfg["pressure"],
                    species=sp,
                )
            )
        elif bc_type == "henrys":
            bcs.append(
                F.HenrysBC(
                    subdomain=sub,
                    H_0=cfg["H_0"],
                    E_H=cfg["E_H"],
                    pressure=cfg["pressure"],
                    species=sp,
                )
            )
        else:
            raise ValueError(f"Unknown FESTIM boundary condition type: '{bc_type}'")

    return bcs


def _build_exports(F, export_cfgs: list, species_map: dict, volumes_map: dict,
                   surfaces_map: dict, milestones: list):
    """Build export objects from config dicts.

    ``"times": "milestones"`` is resolved to the stepsize milestones list.
    """
    exports = []
    for cfg in export_cfgs:
        export_type = cfg["type"]
        times = milestones if cfg.get("times") == "milestones" else cfg.get("times")

        if export_type == "vtx_species":
            field = species_map[cfg["field"]]
            subdomain = volumes_map.get(cfg.get("volume_id"))
            exports.append(
                F.VTXSpeciesExport(
                    field=field, filename=cfg["filename"],
                    subdomain=subdomain, times=times,
                )
            )
        elif export_type == "vtx_temperature":
            exports.append(
                F.VTXTemperatureExport(filename=cfg["filename"], times=times)
            )
        elif export_type == "profile_1d":
            field = species_map[cfg["field"]]
            subdomain = volumes_map.get(cfg.get("volume_id"))
            exports.append(
                F.Profile1DExport(field=field, subdomain=subdomain, times=times)
            )
        elif export_type == "surface_flux":
            field = species_map[cfg["field"]]
            surface = surfaces_map[cfg["surface_id"]]
            exports.append(
                F.SurfaceFlux(field=field, surface=surface, filename=cfg.get("filename"))
            )
        elif export_type == "total_volume":
            field = species_map[cfg["field"]]
            volume = volumes_map[cfg["volume_id"]]
            exports.append(
                F.TotalVolume(field=field, volume=volume, filename=cfg.get("filename"))
            )
        else:
            raise ValueError(f"Unknown FESTIM export type: '{export_type}'")

    return exports


# ---------------------------------------------------------------------------
# Strategy registry
# ---------------------------------------------------------------------------

_SIM_TYPE_REGISTRY: dict = {}


class FestimSimStrategy(ABC):
    """Base class for FESTIM simulation type strategies."""

    @abstractmethod
    def build(self, F, unit_cfg: dict, materials_map: dict, helpers):
        """Build and return a FESTIM problem object.

        Args:
            F: The ``festim`` module.
            unit_cfg: The full unit config dict (mesh, species, subdomains, etc.).
            materials_map: Dict mapping material names to ``F.Material`` objects.
            helpers: ``FestimBuildHelpers`` instance (reserved for future use).

        Returns:
            A FESTIM problem instance (e.g. ``F.HydrogenTransportProblem``).
        """


def register_festim_sim_type(name: str, strategy_cls: type) -> None:
    """Register a FESTIM simulation type strategy."""
    _SIM_TYPE_REGISTRY[name] = strategy_cls


def get_registered_sim_types() -> dict[str, type]:
    """Return a copy of the registered FESTIM sim type registry."""
    return dict(_SIM_TYPE_REGISTRY)


# ---------------------------------------------------------------------------
# Built-in strategies
# ---------------------------------------------------------------------------


class _HydrogenTransportStrategy(FestimSimStrategy):
    """Strategy for ``hydrogen_transport`` sim type."""

    def build(self, F, unit_cfg: dict, materials_map: dict, helpers):
        mesh = _build_mesh(F, unit_cfg["mesh"])

        volumes, surfaces = _build_subdomains(F, unit_cfg["subdomains"], materials_map)
        subdomains_list = volumes + surfaces

        species_list = _build_species(F, unit_cfg["species"])
        species_map = {s.name: s for s in species_list}
        surfaces_map = {s.id: s for s in surfaces}
        volumes_map = {v.id: v for v in volumes}

        ss_cfg = unit_cfg.get("solver_config", {}).get("stepsize", {})
        stepsize = F.Stepsize(
            initial_value=ss_cfg["initial_value"],
            growth_factor=ss_cfg.get("growth_factor"),
            cutback_factor=ss_cfg.get("cutback_factor"),
            target_nb_iterations=ss_cfg.get("target_nb_iterations"),
            max_stepsize=ss_cfg.get("max_stepsize"),
            milestones=ss_cfg.get("milestones"),
        )
        milestones = stepsize.milestones or []

        bcs = _build_bcs(F, unit_cfg.get("boundary_conditions", []), species_map, surfaces_map)
        exports = _build_exports(
            F, unit_cfg.get("exports", []), species_map, volumes_map, surfaces_map, milestones
        )

        sc = unit_cfg.get("solver_config", {})
        settings = F.Settings(
            atol=sc["atol"],
            rtol=sc["rtol"],
            final_time=sc.get("final_time"),
            transient=sc.get("transient", True),
            element_degree=sc.get("element_degree", 1),
        )
        settings.stepsize = stepsize

        model = F.HydrogenTransportProblem(
            mesh=mesh,
            subdomains=subdomains_list,
            species=species_list,
            temperature=sc.get("temperature", 300),
            boundary_conditions=bcs,
            settings=settings,
            exports=exports,
        )

        return model


register_festim_sim_type("hydrogen_transport", _HydrogenTransportStrategy)


# ---------------------------------------------------------------------------
# Provider
# ---------------------------------------------------------------------------


class FestimProvider(AbstractProvider):
    """Hydrogen transport FEM provider backed by FESTIM.

    Responsibilities
    ----------------
    * **Material registry** — ``initialize()`` loads all material defs from the
      flowsheet and indexes them by name.
    * **Material validation** — ``validate_material()`` checks that ``D_0`` and
      ``E_D`` are present in ``MaterialDef.extra``.
    * **Simulation dispatch** — ``run_simulation()`` looks up the strategy from
      ``_SIM_TYPE_REGISTRY`` by ``sim_type``, builds the FESTIM problem,
      calls ``initialise()`` and ``run()``, then reports results.
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
                "FESTIM is not installed. Install with: conda install -c conda-forge festim"
            ) from exc

        self._provider_output_dir = provider_config.output_dir
        self._materials = {}
        for mat_name, mat_def in flowsheet_config.materials.items():
            self._materials[mat_name] = mat_def

        self._initialized = True
        logger.info(
            f"FestimProvider initialized with {len(self._materials)} material(s). "
            f"Registered sim_types: {sorted(_SIM_TYPE_REGISTRY)}"
        )

    def get_thermo_properties(self, stream: dict) -> dict:
        raise NotImplementedError("FestimProvider does not support stream thermodynamics.")

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
        * ``extra.D_0`` is required (diffusion pre-exponential).
        * ``extra.E_D`` is required (diffusion activation energy).

        Returns:
            List of error strings (empty = valid).
        """
        errors = []
        extra = mat_def.extra or {}

        if "D_0" not in extra:
            errors.append(
                f"Material '{mat_name}' is missing 'extra.D_0' "
                "(diffusion pre-exponential, required for FESTIM)."
            )
        if "E_D" not in extra:
            errors.append(
                f"Material '{mat_name}' is missing 'extra.E_D' "
                "(diffusion activation energy, required for FESTIM)."
            )

        return errors

    def run_simulation(self, unit_config: "UnitConfig", inlet: dict) -> "SimulationResult":
        """Build and run a FESTIM simulation from a typed ``UnitConfig``.

        Execution flow
        --------------
        1. Look up strategy from ``_SIM_TYPE_REGISTRY``
        2. Build ``F.Material`` objects for all registry materials
        3. Call ``strategy.build()`` → problem instance
        4. Create output directory
        5. Call ``model.initialise()`` then ``model.run()``
        6. Extract results from exports
        7. Return ``SimulationResult``
        """
        import festim as F

        from processforge.types import SimulationResult

        sim_type = unit_config.sim_type
        strategy_cls = _SIM_TYPE_REGISTRY.get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"FestimProvider: unknown sim_type '{sim_type}'. "
                f"Register with register_festim_sim_type(). "
                f"Built-in types: {sorted(_SIM_TYPE_REGISTRY)}"
            )

        materials_map = {}
        for mat_name, mdef in self._materials.items():
            materials_map[mat_name] = _build_material(F, mat_name, mdef)

        logger.info(
            f"FestimProvider: running sim_type='{sim_type}' "
            f"with strategy {strategy_cls.__name__}"
        )

        full_cfg = {**unit_config.extra, "solver_config": unit_config.solver_config or {}}
        model = strategy_cls().build(F, full_cfg, materials_map, None)

        run_dir = pathlib.Path(self._provider_output_dir)
        run_dir.mkdir(parents=True, exist_ok=True)

        model.initialise()
        model.run()

        logger.info(f"FestimProvider: '{sim_type}' completed in '{run_dir}'")

        scalars, metadata = self._extract_results(model, run_dir)

        return SimulationResult(
            status="completed",
            sim_type=sim_type,
            scalars=scalars,
            metadata=metadata,
        )

    def _extract_results(self, model, run_dir: pathlib.Path) -> tuple:
        """Extract scalar results from the FESTIM model's derived quantities.

        Reads CSV output files from DerivedQuantity exports and extracts
        the final value.  Returns ``(scalars_dict, metadata_dict)``.
        """
        scalars: dict = {}
        metadata: dict = {"run_dir": str(run_dir.resolve())}

        for export in getattr(model, "exports", []):
            if not hasattr(export, "filename") or export.filename is None:
                continue

            filename = str(export.filename)
            csv_path = run_dir / filename

            if not csv_path.exists():
                # FESTIM writes derived quantities to CWD, which may differ
                # from run_dir if the model was built before chdir.
                # Try the current working directory as fallback.
                csv_path = pathlib.Path(filename)

            if csv_path.exists() and csv_path.suffix == ".csv":
                try:
                    import csv as _csv

                    with open(csv_path, newline="") as f:
                        reader = _csv.reader(f)
                        header = next(reader, None)
                        last_row = None
                        for row in reader:
                            last_row = row
                        if last_row and header:
                            # DerivedQuantity CSVs have columns: t(s), quantity
                            # Use the export class name + filename as key
                            key = pathlib.Path(filename).stem
                            scalars[key] = float(last_row[-1])
                            metadata[f"csv_{key}"] = str(csv_path.resolve())
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        f"FestimProvider: could not parse CSV '{csv_path}': {exc}"
                    )
            elif csv_path.exists():
                metadata[f"file_{pathlib.Path(filename).stem}"] = str(csv_path.resolve())

        return scalars, metadata


# Register the provider
register_provider("festim", FestimProvider)
