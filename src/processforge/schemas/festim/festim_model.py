"""Typed Pydantic models for the FESTIM ``solver_config`` block.

These models are the single source of truth used by both flowsheet validation
and the runtime ``FestimProvider``.  They mirror the constructor signatures of
the FESTIM 2.0 API (see ``flowsheets/festim/FESTIM-main/src/festim``)::

    Mesh1D, Material, Species, ImplicitSpecies, VolumeSubdomain1D,
    SurfaceSubdomain1D, Reaction, ParticleSource,
    FixedConcentrationBC, ParticleFluxBC, SievertsBC, HenrysBC,
    Settings, Stepsize, SurfaceFlux, TotalVolume, Profile1DExport

Expression-like FESTIM arguments (source values, temperatures, BC values,
``max_stepsize``) cannot be serialised as Python callables; they are expressed
as structured profile objects (``GaussianProfile``, ``RampProfile``,
``AfterTProfile``) which the provider converts back into UFL expressions /
callables at build time.
"""

from enum import Enum
from typing import Annotated, List, Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator


class SolubilityLaw(str, Enum):
    """Solubility law of a material (mirrors ``festim.SolubilityLaw``)."""

    SIEVERT = "sievert"
    HENRY = "henry"
    NONE = "none"


# ---------------------------------------------------------------------------
# Expression-like value profiles
# ---------------------------------------------------------------------------


class GaussianProfile(BaseModel):
    """Gaussian spatial source profile.

    Builds the UFL expression ``amplitude / (width * sqrt(2*pi)) *
    exp(-0.5*((x - center)/width)^2)``, optionally gated to
    ``t <= active_until`` (``None`` = always on).
    """

    model_config = ConfigDict(extra="forbid")

    type: Literal["gaussian"] = "gaussian"
    amplitude: float = Field(description="Peak amplitude of the source (mol m-3 s-1).")
    center: float = Field(description="Gaussian centre position (m).")
    width: float = Field(gt=0, description="Gaussian standard deviation (m).")
    active_until: Optional[float] = Field(
        default=None,
        description="Time (s) until which the source is active; None = always on.",
    )


class RampProfile(BaseModel):
    """Linear ramp in time.

    Builds the callable ``T0`` for ``t < start_t`` else
    ``T0 + ramp_rate * (t - start_t)``.
    """

    model_config = ConfigDict(extra="forbid")

    type: Literal["ramp"] = "ramp"
    T0: float = Field(description="Value before the ramp starts.")
    start_t: float = Field(description="Time (s) at which the ramp starts.")
    ramp_rate: float = Field(description="Rate of change per second.")


class AfterTProfile(BaseModel):
    """Piecewise-constant profile active after a time.

    Builds the callable ``value`` for ``t > after_t`` else ``None`` — used for
    FESTIM's callable ``max_stepsize``.
    """

    model_config = ConfigDict(extra="forbid")

    type: Literal["after_t"] = "after_t"
    after_t: float = Field(description="Time (s) after which ``value`` applies.")
    value: float = Field(description="Value to return after ``after_t``.")


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------


class LinspaceSegment(BaseModel):
    """A single ``np.linspace`` mesh segment (mirrors the FESTIM tutorial mesh)."""

    model_config = ConfigDict(extra="forbid")

    start: float = Field(description="First vertex coordinate (m).")
    stop: float = Field(description="Last vertex coordinate (m).")
    num: int = Field(gt=1, description="Number of vertices in the segment.")


class MeshConfig(BaseModel):
    """1-D mesh definition.

    Exactly one of ``vertices`` or ``segments`` must be given.

    * ``vertices`` — explicit coordinates, either a flat list or a list of
      blocks (each block a disconnected segment of the mesh).
    * ``segments`` — ``np.linspace`` segments concatenated in order (this is
      how the FESTIM TDS tutorial builds its refined mesh).
    """

    model_config = ConfigDict(extra="forbid")

    vertices: Optional[Union[List[float], List[List[float]]]] = Field(
        default=None,
        description="Explicit vertex x-coordinates (m); flat list or list of blocks.",
    )
    segments: Optional[List[LinspaceSegment]] = Field(
        default=None,
        description="np.linspace segments concatenated in order.",
    )

    @model_validator(mode="after")
    def _validate_exactly_one(self) -> "MeshConfig":
        if (self.vertices is None) == (self.segments is None):
            raise ValueError("mesh must define exactly one of 'vertices' or 'segments'")
        if self.vertices is not None:
            if len(self.vertices) < 2:
                raise ValueError("mesh.vertices must contain at least 2 points")
            if isinstance(self.vertices[0], list) and len(self.vertices[0]) < 2:
                raise ValueError(
                    "each block of mesh.vertices must contain at least 2 points"
                )
        return self


# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------


class MaterialConfig(BaseModel):
    """FESTIM material properties carried in the flowsheet material ``extra``.

    Mirrors ``festim.Material(D_0, E_D, K_S_0, E_K_S, ...)``.  ``D_0`` and
    ``E_D`` are required; ``K_S_0``/``E_K_S`` must be given together.
    """

    model_config = ConfigDict(extra="forbid")

    name: Optional[str] = Field(default=None, description="Material name.")
    D_0: float = Field(gt=0, description="Diffusion pre-exponential factor (m2 s-1).")
    E_D: float = Field(ge=0, description="Diffusion activation energy (eV).")
    K_S_0: Optional[float] = Field(
        default=None,
        description="Solubility pre-exponential factor (H m-3 Pa-0.5).",
    )
    E_K_S: Optional[float] = Field(
        default=None, description="Solubility activation energy (eV)."
    )
    thermal_conductivity: Optional[float] = Field(
        default=None, description="Thermal conductivity (W m-1 K-1)."
    )
    density: Optional[float] = Field(default=None, description="Density (kg m-3).")
    heat_capacity: Optional[float] = Field(
        default=None, description="Heat capacity (J kg-1 K-1)."
    )
    solubility_law: SolubilityLaw = Field(
        default=SolubilityLaw.NONE,
        description="Solubility law: 'sievert', 'henry' or 'none'.",
    )

    @model_validator(mode="after")
    def _validate_solubility_pairing(self) -> "MaterialConfig":
        has_K = self.K_S_0 is not None
        has_E = self.E_K_S is not None
        if has_K != has_E:
            raise ValueError("'K_S_0' and 'E_K_S' must be provided together")
        return self


# ---------------------------------------------------------------------------
# Species
# ---------------------------------------------------------------------------


class SpeciesConfig(BaseModel):
    """Mobile or immobile hydrogen species (mirrors ``festim.Species``)."""

    model_config = ConfigDict(extra="forbid")

    name: str = Field(description="Unique species name.")
    mobile: bool = Field(default=True, description="Whether the species diffuses.")


class ImplicitSpeciesConfig(BaseModel):
    """Implicit species: ``c = n - sum(others)`` (mirrors ``festim.ImplicitSpecies``).

    Used to model empty trap sites: ``n`` is the total trap-site density and
    ``others`` the trapped species whose concentration is subtracted.
    """

    model_config = ConfigDict(extra="forbid")

    name: str = Field(description="Unique name of the implicit species.")
    n: float = Field(description="Total concentration of the implicit species.")
    others: List[str] = Field(
        default_factory=list,
        description="Names of species subtracted from ``n``.",
    )


# ---------------------------------------------------------------------------
# Subdomains
# ---------------------------------------------------------------------------


class VolumeSubdomainConfig(BaseModel):
    """1-D volume subdomain (mirrors ``festim.VolumeSubdomain1D``)."""

    model_config = ConfigDict(extra="forbid")

    id: int = Field(
        description="Unique, non-zero volume subdomain id.",
    )
    borders: List[float] = Field(
        min_length=2,
        max_length=2,
        description="[x0, x1] borders of the volume (m).",
    )
    material: str = Field(
        description="Name of a material defined in the flowsheet materials registry.",
    )


class SurfaceSubdomainConfig(BaseModel):
    """1-D surface subdomain (mirrors ``festim.SurfaceSubdomain1D``)."""

    model_config = ConfigDict(extra="forbid")

    id: int = Field(description="Unique surface subdomain id.")
    x: float = Field(description="x-coordinate of the surface (m).")


class SubdomainsConfig(BaseModel):
    """Volume and surface subdomains of the problem."""

    model_config = ConfigDict(extra="forbid")

    volume: List[VolumeSubdomainConfig] = Field(
        description="Volume subdomains.",
    )
    surface: List[SurfaceSubdomainConfig] = Field(
        description="Surface subdomains.",
    )


# ---------------------------------------------------------------------------
# Reactions
# ---------------------------------------------------------------------------


class ArrheniusReactionConfig(BaseModel):
    """Arrhenius trapping/detrapping reaction (mirrors ``festim.ArrheniusReaction``).

    ``reactant``/``product`` reference species or implicit-species names.
    ``p_0`` and ``E_p`` are required exactly when ``product`` is non-empty.
    """

    model_config = ConfigDict(extra="forbid")

    type: Literal["arrhenius"] = "arrhenius"
    reactant: List[str] = Field(
        min_length=1,
        description="Reactant species/implicit-species names.",
    )
    product: List[str] = Field(
        default_factory=list,
        description="Product species names (empty for irreversible reactions).",
    )
    k_0: float = Field(gt=0, description="Forward rate pre-exponential factor.")
    E_k: float = Field(ge=0, description="Forward rate activation energy (eV).")
    p_0: Optional[float] = Field(
        default=None, description="Backward rate pre-exponential factor."
    )
    E_p: Optional[float] = Field(
        default=None, description="Backward rate activation energy (eV)."
    )
    volume: int = Field(description="Volume subdomain id where the reaction occurs.")

    @model_validator(mode="after")
    def _validate_product_rates(self) -> "ArrheniusReactionConfig":
        has_product = len(self.product) > 0
        if has_product and (self.p_0 is None or self.E_p is None):
            raise ValueError("'p_0' and 'E_p' are required when 'product' is non-empty")
        if not has_product and (self.p_0 is not None or self.E_p is not None):
            raise ValueError("'p_0'/'E_p' must be omitted when there is no product")
        return self


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------


class SourceConfig(BaseModel):
    """Volumetric particle source (mirrors ``festim.ParticleSource``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["particle"] = "particle"
    species: str = Field(description="Species name the source feeds.")
    volume: int = Field(description="Volume subdomain id of the source.")
    value: Union[float, GaussianProfile] = Field(
        description="Constant source rate or a structured Gaussian profile.",
    )


# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------


class FixedConcentrationBCConfig(BaseModel):
    """``c = value`` on a surface (mirrors ``festim.FixedConcentrationBC``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["fixed_concentration"] = "fixed_concentration"
    subdomain: int = Field(description="Surface subdomain id.")
    species: str = Field(description="Species name.")
    value: Union[float, RampProfile] = Field(
        description="Concentration (mol m-3) or time-ramped value.",
    )


class ParticleFluxBCConfig(BaseModel):
    """Flux boundary condition (mirrors ``festim.ParticleFluxBC``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["particle_flux"] = "particle_flux"
    subdomain: int = Field(description="Surface subdomain id.")
    species: str = Field(description="Species name.")
    value: Union[float, RampProfile] = Field(
        description="Flux (mol m-2 s-1) or time-ramped value.",
    )


class SievertsBCConfig(BaseModel):
    """Sieverts' law boundary condition (mirrors ``festim.SievertsBC``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["sieverts"] = "sieverts"
    subdomain: int = Field(description="Surface subdomain id.")
    species: str = Field(description="Species name.")
    S_0: float = Field(
        gt=0, description="Sieverts constant pre-exponential factor (H m-3 Pa-0.5)."
    )
    E_S: float = Field(ge=0, description="Sieverts constant activation energy (eV).")
    pressure: float = Field(description="Gas pressure at the boundary (Pa).")


class HenrysBCConfig(BaseModel):
    """Henry's law boundary condition (mirrors ``festim.HenrysBC``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["henrys"] = "henrys"
    subdomain: int = Field(description="Surface subdomain id.")
    species: str = Field(description="Species name.")
    H_0: float = Field(
        gt=0, description="Henry's constant pre-exponential factor (H m-3 Pa-1)."
    )
    E_H: float = Field(ge=0, description="Henry's constant activation energy (eV).")
    pressure: float = Field(description="Gas pressure at the boundary (Pa).")


BoundaryConditionConfig = Annotated[
    Union[
        FixedConcentrationBCConfig,
        ParticleFluxBCConfig,
        SievertsBCConfig,
        HenrysBCConfig,
    ],
    Field(discriminator="type"),
]


# ---------------------------------------------------------------------------
# Temperature
# ---------------------------------------------------------------------------


class ConstantTemperatureConfig(BaseModel):
    """Constant temperature in K."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["constant"] = "constant"
    value: float = Field(description="Temperature (K).")


class RampTemperatureConfig(BaseModel):
    """Time-ramped temperature (mirrors the TDS tutorial's piecewise ramp)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["ramp"] = "ramp"
    T0: float = Field(description="Temperature before the ramp (K).")
    start_t: float = Field(description="Time (s) at which the ramp starts.")
    ramp_rate: float = Field(description="Heating rate (K s-1).")


TemperatureConfig = Annotated[
    Union[ConstantTemperatureConfig, RampTemperatureConfig],
    Field(discriminator="type"),
]


# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------


class SurfaceFluxExportConfig(BaseModel):
    """Surface flux derived quantity (mirrors ``festim.SurfaceFlux``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["surface_flux"] = "surface_flux"
    field: str = Field(description="Species name to integrate the flux of.")
    surface: int = Field(description="Surface subdomain id.")
    filename: Optional[str] = Field(
        default=None, description="Optional CSV/TXT output filename."
    )


class TotalVolumeExportConfig(BaseModel):
    """Volume integral of a field (mirrors ``festim.TotalVolume``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["total_volume"] = "total_volume"
    field: str = Field(description="Species name to integrate.")
    volume: int = Field(description="Volume subdomain id.")
    filename: Optional[str] = Field(
        default=None, description="Optional CSV/TXT output filename."
    )


class Profile1DExportConfig(BaseModel):
    """1-D concentration profile snapshots (mirrors ``festim.Profile1DExport``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["profile_1d"] = "profile_1d"
    field: str = Field(description="Species name to profile.")
    volume: int = Field(description="Volume subdomain id.")
    times: List[float] = Field(
        description="Times (s) at which profiles are exported.",
    )


ExportConfig = Annotated[
    Union[SurfaceFluxExportConfig, TotalVolumeExportConfig, Profile1DExportConfig],
    Field(discriminator="type"),
]


# ---------------------------------------------------------------------------
# Stepsize / settings
# ---------------------------------------------------------------------------


class StepsizeConfig(BaseModel):
    """Adaptive time-step control (mirrors ``festim.Stepsize``)."""

    model_config = ConfigDict(extra="forbid")

    initial_value: float = Field(gt=0, description="Initial time step (s).")
    growth_factor: Optional[float] = Field(
        default=None,
        description="Step growth factor (>= 1) for fast-converging iterations.",
    )
    cutback_factor: Optional[float] = Field(
        default=None,
        description="Step cut-back factor (<= 1) for slow iterations.",
    )
    target_nb_iterations: Optional[int] = Field(
        default=None, description="Newton iterations the adaptivity targets."
    )
    max_stepsize: Optional[Union[float, AfterTProfile]] = Field(
        default=None,
        description="Maximum step (s), constant or active after a given time.",
    )
    milestones: Optional[List[float]] = Field(
        default=None,
        description="Times (s) the simulation must pass through exactly.",
    )

    @model_validator(mode="after")
    def _validate_adaptive(self) -> "StepsizeConfig":
        if self.milestones and (
            self.growth_factor is None or self.cutback_factor is None
        ):
            raise ValueError(
                "'milestones' require 'growth_factor' and 'cutback_factor' "
                "(adaptive stepsize)"
            )
        if self.growth_factor is not None and self.growth_factor < 1:
            raise ValueError("'growth_factor' must be >= 1")
        if self.cutback_factor is not None and self.cutback_factor > 1:
            raise ValueError("'cutback_factor' must be <= 1")
        if isinstance(self.max_stepsize, float):
            if self.max_stepsize < self.initial_value:
                raise ValueError("'max_stepsize' must not be less than 'initial_value'")
        return self


class InitialConcentrationConfig(BaseModel):
    """Initial concentration of a species in a volume (mirrors
    ``festim.InitialConcentration``)."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["concentration"] = "concentration"
    value: float = Field(description="Initial concentration (mol m-3).")
    species: str = Field(description="Species name.")
    volume: int = Field(description="Volume subdomain id.")


# ---------------------------------------------------------------------------
# Top-level solver-config model
# ---------------------------------------------------------------------------


class FestimModel(BaseModel):
    """Typed representation of the ``solver_config`` block of a FESTIM
    ``SolverUnit``.

    Parsed automatically from the opaque ``solver_config`` JSON dict — this is
    the single source of truth used by both flowsheet validation and the
    runtime :class:`~processforge.providers.festim_provider.FestimProvider`.
    """

    model_config = ConfigDict(extra="forbid")

    mesh: MeshConfig = Field(description="1-D mesh definition.")
    subdomains: SubdomainsConfig = Field(description="Volume and surface subdomains.")
    species: List[SpeciesConfig] = Field(description="Mobile/immobile species.")
    implicit_species: List[ImplicitSpeciesConfig] = Field(
        default_factory=list,
        description="Implicit species (e.g. empty trap sites).",
    )
    reactions: List[ArrheniusReactionConfig] = Field(
        default_factory=list,
        description="Trapping/detrapping reactions.",
    )
    sources: List[SourceConfig] = Field(
        default_factory=list,
        description="Volumetric particle sources.",
    )
    initial_conditions: List[InitialConcentrationConfig] = Field(
        default_factory=list,
        description="Initial concentration conditions.",
    )
    boundary_conditions: List[BoundaryConditionConfig] = Field(
        default_factory=list,
        description="Boundary conditions.",
    )
    temperature: Union[float, TemperatureConfig] = Field(
        description="Temperature (K): a constant float or a structured profile.",
    )
    exports: List[ExportConfig] = Field(
        default_factory=list,
        description="Derived quantities / profile exports.",
    )
    atol: float = Field(gt=0, description="Absolute solver tolerance.")
    rtol: float = Field(gt=0, description="Relative solver tolerance.")
    max_iterations: int = Field(
        default=30, gt=0, description="Maximum Newton iterations."
    )
    transient: bool = Field(
        default=True, description="Whether the simulation is time-dependent."
    )
    final_time: Optional[float] = Field(
        default=None, description="End time of a transient simulation (s)."
    )
    element_degree: int = Field(default=1, ge=1, description="Finite-element degree.")
    stepsize: Optional[StepsizeConfig] = Field(
        default=None, description="Time-step control for transient runs."
    )
    convergence_criterion: Literal["residual", "incremental"] = Field(
        default="residual",
        description="Newton convergence criterion.",
    )

    # -- cross-field validation ------------------------------------------

    @model_validator(mode="after")
    def _validate_references(self) -> "FestimModel":
        species_names = {s.name for s in self.species}
        if len(species_names) != len(self.species):
            raise ValueError(
                "species names must be unique; got " f"{[s.name for s in self.species]}"
            )
        if not species_names:
            raise ValueError("at least one species must be defined")

        implicit_names = {s.name for s in self.implicit_species}
        if implicit_names & species_names:
            raise ValueError(
                "implicit_species names must not collide with species names; "
                f"collisions: {sorted(implicit_names & species_names)}"
            )

        all_names = species_names | implicit_names

        volume_ids = {v.id for v in self.subdomains.volume}
        if len(volume_ids) != len(self.subdomains.volume):
            raise ValueError("volume subdomain ids must be unique")
        if not volume_ids:
            raise ValueError("at least one volume subdomain must be defined")
        if any(vid == 0 for vid in volume_ids):
            raise ValueError("volume subdomain ids cannot be 0")

        surface_ids = {s.id for s in self.subdomains.surface}
        if len(surface_ids) != len(self.subdomains.surface):
            raise ValueError("surface subdomain ids must be unique")
        if not surface_ids:
            raise ValueError("at least one surface subdomain must be defined")

        if not any(s.mobile for s in self.species):
            raise ValueError("at least one mobile species is required")

        # species references
        for source in self.sources:
            if source.species not in species_names:
                raise ValueError(
                    f"source for species '{source.species}' references an "
                    f"undefined species; available: {sorted(species_names)}"
                )
            if source.volume not in volume_ids:
                raise ValueError(
                    f"source volume={source.volume} is not a defined volume "
                    f"subdomain; available: {sorted(volume_ids)}"
                )

        for ic in self.initial_conditions:
            if ic.species not in species_names:
                raise ValueError(
                    f"initial condition species '{ic.species}' is undefined; "
                    f"available: {sorted(species_names)}"
                )
            if ic.volume not in volume_ids:
                raise ValueError(
                    f"initial condition volume={ic.volume} is not a defined "
                    f"volume subdomain; available: {sorted(volume_ids)}"
                )

        for bc in self.boundary_conditions:
            if bc.species not in species_names:
                raise ValueError(
                    f"boundary condition species '{bc.species}' is undefined; "
                    f"available: {sorted(species_names)}"
                )
            if bc.subdomain not in surface_ids:
                raise ValueError(
                    f"boundary condition subdomain={bc.subdomain} is not a "
                    f"defined surface subdomain; available: {sorted(surface_ids)}"
                )

        for reaction in self.reactions:
            for name in reaction.reactant + reaction.product:
                if name not in all_names:
                    raise ValueError(
                        f"reaction references species '{name}' which is neither "
                        f"a species nor an implicit species; available: "
                        f"{sorted(all_names)}"
                    )
            if reaction.volume not in volume_ids:
                raise ValueError(
                    f"reaction volume={reaction.volume} is not a defined volume "
                    f"subdomain; available: {sorted(volume_ids)}"
                )

        for export in self.exports:
            if export.field not in species_names:
                raise ValueError(
                    f"export field '{export.field}' is undefined; available: "
                    f"{sorted(species_names)}"
                )
            if hasattr(export, "volume") and export.volume not in volume_ids:
                raise ValueError(
                    f"export volume={export.volume} is not a defined volume "
                    f"subdomain; available: {sorted(volume_ids)}"
                )
            if hasattr(export, "surface") and export.surface not in surface_ids:
                raise ValueError(
                    f"export surface={export.surface} is not a defined surface "
                    f"subdomain; available: {sorted(surface_ids)}"
                )

        # transient consistency
        if self.transient:
            if self.final_time is None:
                raise ValueError("'final_time' is required for a transient simulation")
            if self.stepsize is None:
                raise ValueError("'stepsize' is required for a transient simulation")

        return self
