"""OpenMC model-building utilities shared by simulation strategies."""
from __future__ import annotations

from processforge.schemas.openmc.openmc_model import (
    MeshTallyConfig,
    OpenMCSetting,
    ReactorCoreGeometryConfig,
    RunMode,
    SourcePoint,
)


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
