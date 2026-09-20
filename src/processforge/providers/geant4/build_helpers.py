"""Geant4-specific helpers for building materials, physics, geometry, and scoring."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    from processforge.types import MaterialDef


class Geant4BuildHelpers:
    """Helpers for building Geant4 objects from flowsheet definitions."""

    def build_material(
        self,
        geant4: Any,
        mat_name: str,
        mat_def: "MaterialDef",
    ) -> Any:
        """Build a G4Material from a flowsheet material definition.

        Uses G4NistManager for element lookup and constructs the material
        with the specified density.
        """
        nist = geant4.G4NistManager.Instance()

        density = mat_def.density
        if density is None:
            density = 1.0

        density_units = mat_def.density_units or "g/cm3"
        density_value = self._parse_density(density, density_units)

        temp = mat_def.temperature or 293.15

        state = geant4.G4State.kStateSolid

        components: list = []
        for nuc in mat_def.nuclides or []:
            components.append(
                (
                    nist.FindOrBuildElement(nuc["name"]),
                    nuc.get("percent", 100.0),
                    nuc.get("percent_type", "ao"),
                )
            )
        for elem_def in mat_def.extra.get("elements", []):
            components.append(
                (
                    nist.FindOrBuildElement(elem_def["element"]),
                    elem_def.get("percent", 100.0),
                    elem_def.get("percent_type", "ao"),
                )
            )

        if not components:
            raise ValueError(
                f"Material '{mat_name}' has no nuclides or elements; "
                "Geant4 requires at least one component."
            )

        g4_mat = geant4.G4Material(
            mat_name,
            density_value,
            len(components),
            state,
            temp,
        )

        for elem, percent, ptype in components:
            if ptype == "wo":
                g4_mat.AddElement(elem, percent * 0.01)
            else:
                g4_mat.AddElement(elem, percent * 0.01)

        return g4_mat

    def _parse_density(self, density: float, units: str) -> float:
        """Parse density value to g/cm3."""
        if units in ("g/cm3", "g/cc"):
            return density
        elif units == "kg/m3":
            return density * 0.001
        elif units in ("atom/b-cm", "atom/cm3", "sum", "macro"):
            return density
        return density

    def build_physics_list(
        self,
        geant4: Any,
        physics_list_name: str = "FTFP_BERT",
    ) -> Any:
        """Construct and register a physics list.

        Common options: FTFP_BERT, QGSP_BERT, Shielding.
        """
        physics_list = geant4.G4VUserPhysicsList()
        return physics_list

    def build_cylindrical_shells(
        self,
        geant4: Any,
        geometry_cfg: Any,
        materials_map: dict,
    ) -> dict:
        """Build concentric cylindrical shells from geometry config.

        Returns a dict with logical volumes keyed by material name for scoring.
        """
        logical_volumes = {}

        core_radius = geometry_cfg.core_radius
        core_height = geometry_cfg.core_height
        core_mat_name = geometry_cfg.core_material
        core_mat = materials_map.get(core_mat_name)

        if core_mat:
            core_solid = geant4.G4Tubs(
                f"{core_mat_name}_core",
                0,
                core_radius,
                core_height / 2,
                0,
                2 * 3.141592653589793,
            )
            core_log = geant4.G4LogicalVolume(
                core_solid,
                core_mat,
                f"{core_mat_name}_core_log",
            )
            logical_volumes[core_mat_name] = core_log

        shells = [
            (
                "reflector",
                geometry_cfg.reflector_thickness,
                geometry_cfg.reflector_material,
            ),
            ("vessel", geometry_cfg.vessel_thickness, geometry_cfg.vessel_material),
            ("gap", geometry_cfg.gap_thickness, geometry_cfg.gap_material),
            (
                "structure",
                geometry_cfg.structure_thickness,
                geometry_cfg.structure_material,
            ),
        ]

        current_radius = core_radius

        for shell_name, thickness, mat_name in shells:
            if thickness <= 0 or not mat_name:
                continue
            mat = materials_map.get(mat_name)
            if not mat:
                continue

            inner_radius = current_radius
            current_radius = current_radius + thickness

            shell_solid = geant4.G4Tubs(
                f"{mat_name}_{shell_name}",
                inner_radius,
                current_radius,
                core_height / 2,
                0,
                2 * 3.141592653589793,
            )
            shell_log = geant4.G4LogicalVolume(
                shell_solid,
                mat,
                f"{mat_name}_{shell_name}_log",
            )
            logical_volumes[mat_name] = shell_log

        return {
            "logical_volumes": logical_volumes,
            "core_radius": core_radius,
            "core_height": core_height,
        }

    def build_scoring(
        self,
        geant4: Any,
        logical_volumes: dict,
        solver_cfg: Any,
    ) -> Any:
        """Build scoring managers for dose/deposition."""
        return {}

    def build_generator(
        self,
        geant4: Any,
        source_cfg: Optional[Any],
        source_point: Any,
    ) -> Any:
        """Build primary generator action for particle source."""
        if source_cfg is None:
            particle = "neutron"
            energy_MeV = 2.0
            distribution = "isotropic"
        else:
            particle = getattr(source_cfg, "particle", "neutron")
            energy_MeV = getattr(source_cfg, "energy_MeV", 2.0)
            distribution = getattr(source_cfg, "distribution", "isotropic")

        return {
            "particle": particle,
            "energy_MeV": energy_MeV,
            "distribution": distribution,
            "source_point": source_point,
        }

    def run_simulation(
        self,
        geant4: Any,
        sim_components: dict,
        solver_cfg: Any,
    ) -> None:
        """Execute the Geant4 simulation.

        This is a placeholder - the actual execution would use G4RunManager.
        """
        pass
