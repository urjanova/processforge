"""OpenMC simulation-type strategies.

Each subclass of :class:`OpenMCSimStrategy` encapsulates the model setup logic
for one ``sim_type`` value.  Strategies are registered via
:func:`register_openmc_sim_type` and looked up by
:func:`get_registered_sim_types`.
"""
from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING, Optional

from processforge.providers._sim_strategy import (
    SimStrategy,
    get_registered_sim_types as _get_registered,
    register_sim_type as _register,
)
from processforge.schemas.openmc.openmc_model import (
    OpenMCSetting,
    PointSourceGeometryConfig,
    ReactorCoreGeometryConfig,
)

if TYPE_CHECKING:
    from processforge.providers.openmc.build_helpers import OpenMCBuildHelpers


ENGINE = "openmc"


class OpenMCSimStrategy(SimStrategy):
    """Base class for an OpenMC simulation type.

    Subclasses declare ``config_model`` and implement ``build``.  The provider
    validates ``geometry_config`` before calling ``build``.
    """

    sim_type: str = ""
    run_mode: str = "eigenvalue"

    @abstractmethod
    def build(
        self,
        openmc,
        solver_cfg: OpenMCSetting,
        geometry_cfg,
        materials_map: dict,
        helpers: "OpenMCBuildHelpers",
    ) -> tuple:
        """Set up the OpenMC model and return ``(materials, geometry, settings, tallies)``."""


def register_openmc_sim_type(name: str, strategy_cls: type) -> None:
    """Register an OpenMC simulation type by name."""
    _register(ENGINE, name, strategy_cls)


def get_registered_sim_types() -> dict[str, type]:
    """Return a copy of the OpenMC sim_type → strategy mapping."""
    return _get_registered(ENGINE)


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
    from processforge.schemas.openmc.openmc_model import MeshTallyConfig

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


class _PointSourceSphereStrategy(OpenMCSimStrategy):
    """Shared build logic for CSG-sphere point-source simulations."""

    config_model = PointSourceGeometryConfig
    sim_type: str = "point_source_csg"
    run_mode: str = "eigenvalue"

    def build(
        self,
        openmc,
        solver_cfg: OpenMCSetting,
        geometry_cfg: PointSourceGeometryConfig,
        materials_map: dict,
        helpers: "OpenMCBuildHelpers",
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
    """Fixed-source point-source approximation using a simple CSG sphere geometry."""

    sim_type = "fixed_source_point"
    run_mode = "fixed source"


class _EigenvalueCSGStrategy(_PointSourceSphereStrategy):
    """Eigenvalue (criticality) simulation using a simple CSG sphere geometry."""

    sim_type = "eigenvalue_csg"
    run_mode = "eigenvalue"


class _ReactorCoreStrategy(OpenMCSimStrategy):
    """Approximate reactor-core simulation from a ``reactor_core`` geometry block."""

    config_model = ReactorCoreGeometryConfig
    sim_type: str = "reactor_core"
    run_mode: str = "eigenvalue"

    def build(
        self,
        openmc,
        solver_cfg: OpenMCSetting,
        geometry_cfg: ReactorCoreGeometryConfig,
        materials_map: dict,
        helpers: "OpenMCBuildHelpers",
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
