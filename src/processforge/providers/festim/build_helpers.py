"""FESTIM model-building utilities shared by simulation strategies."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
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

if TYPE_CHECKING:
    from processforge.types import MaterialDef


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
    """Build a ``(x, t) -> UFL expr`` callable for a Gaussian particle source."""

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


class FestimBuildHelpers:
    """Shared FESTIM object-construction utilities injected into sim strategies."""

    # Optional FESTIM material properties read from ``extra``. FESTIM accepts
    # ``None`` for all of these *except* ``solubility_law`` (its setter raises on
    # ``None``), so that one needs a real default.
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
        vertices = mesh_cfg.vertices
        if vertices is None:
            blocks = [
                np.linspace(seg.start, seg.stop, seg.num) for seg in mesh_cfg.segments
            ]
            vertices = np.concatenate(blocks)
        return festim.Mesh1D(vertices)

    @staticmethod
    def build_species(festim, model: FestimModel) -> tuple:
        """Construct ``festim.Species`` / ``festim.ImplicitSpecies`` objects."""
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
        """Construct volume + surface subdomains."""
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
                kwargs.update(milestones=ss.milestones)
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
