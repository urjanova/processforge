"""Pydantic models for Geant4 provider configuration and settings."""

from __future__ import annotations

from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator


class SourcePoint(BaseModel):
    """Point source location and optional energy."""

    x: float = Field(default=0.0, description="X coordinate in cm.")
    y: float = Field(default=0.0, description="Y coordinate in cm.")
    z: float = Field(default=0.0, description="Z coordinate in cm.")
    energy_MeV: Optional[float] = Field(
        default=None, description="Source energy in MeV."
    )


class SourceSettings(BaseModel):
    """Particle source configuration for Geant4 simulations."""

    particle: str = Field(
        default="neutron", description="Particle type (neutron, gamma, etc.)."
    )
    energy_MeV: float = Field(default=2.0, description="Source energy in MeV.")
    distribution: Literal["isotropic", "beam", "cosine"] = Field(
        default="isotropic", description="Angular distribution."
    )
    position: Optional[SourcePoint] = Field(
        default=None, description="Source position."
    )


class Geant4Setting(BaseModel):
    """Solver settings for Geant4 simulations."""

    physics_list: str = Field(
        default="FTFP_BERT",
        description="Geant4 physics list (e.g., FTFP_BERT, QGSP_BERT, Shielding).",
    )
    events: int = Field(
        default=1000,
        description="Number of events to simulate.",
    )
    seed: Optional[int] = Field(
        default=None,
        description="Random seed for reproducibility.",
    )
    tracking_cutoff: Optional[float] = Field(
        default=None,
        description="Tracking cutoff energy in MeV.",
    )


class ReactorCoreGeometryConfig(BaseModel):
    """Approximate cylindrical reactor-core geometry for shielding simulations.

    Nested cylinders: a salt ``core`` surrounded by optional ``reflector``,
    ``vessel``, ``gap``, and outer ``structure`` shells. Uses the declared
    flowsheet materials so all of them participate in the simulation.
    """

    model_config = ConfigDict(extra="forbid")

    type: str = Field(
        default="reactor_core", description="Geometry kind discriminator."
    )
    core_radius: float = Field(description="Core (salt) radius in cm.")
    core_height: float = Field(description="Core (salt) height in cm.")
    reflector_thickness: float = Field(
        default=0.0, description="Reflector shell thickness in cm."
    )
    vessel_thickness: float = Field(
        default=0.0, description="Vessel shell thickness in cm."
    )
    gap_thickness: float = Field(default=0.0, description="Gap shell thickness in cm.")
    structure_thickness: float = Field(
        default=0.0, description="Outer structure shell thickness in cm."
    )
    core_material: str = Field(description="Material name for the core.")
    reflector_material: Optional[str] = Field(
        default=None, description="Material name for the reflector shell."
    )
    vessel_material: Optional[str] = Field(
        default=None, description="Material name for the vessel shell."
    )
    gap_material: Optional[str] = Field(
        default=None, description="Material name for the gap shell."
    )
    structure_material: Optional[str] = Field(
        default=None, description="Material name for the outer structure shell."
    )
    source_point: Optional[SourcePoint] = Field(
        default=None, description="Point source location (and optional energy)."
    )
    source: Optional[SourceSettings] = Field(
        default=None, description="Particle source configuration."
    )

    @model_validator(mode="after")
    def _check_materials_ctx(self, info) -> "ReactorCoreGeometryConfig":
        materials = (info.context or {}).get("materials")
        if not materials:
            return self
        refs = {
            "core_material": self.core_material,
            "reflector_material": self.reflector_material,
            "vessel_material": self.vessel_material,
            "gap_material": self.gap_material,
            "structure_material": self.structure_material,
        }
        for label, name in refs.items():
            if name is not None and name not in materials:
                raise ValueError(
                    f"{label}='{name}' is not a declared flowsheet material. "
                    f"Available: {sorted(materials)}"
                )
        return self
