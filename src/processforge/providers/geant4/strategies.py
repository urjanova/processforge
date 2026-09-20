"""Geant4 simulation-type strategies.

Each subclass of :class:`Geant4SimStrategy` encapsulates the model setup logic
for one ``sim_type`` value.  Strategies are registered via
:func:`register_geant4_sim_type` and looked up by
:func:`get_registered_sim_types`.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING

from processforge.providers._sim_strategy import (
    SimStrategy,
    get_registered_sim_types as _get_registered,
    register_sim_type as _register,
)

if TYPE_CHECKING:
    from processforge.providers.geant4.build_helpers import Geant4BuildHelpers
    from processforge.schemas.geant4.geant4_model import Geant4Setting


ENGINE = "geant4"


class Geant4SimStrategy(SimStrategy):
    """Base class for a Geant4 simulation type.

    Subclasses declare ``config_model`` (optional) and implement ``build``.
    """

    sim_type: str = ""

    @abstractmethod
    def build(
        self,
        geant4,
        solver_cfg: "Geant4Setting",
        geometry_cfg,
        materials_map: dict,
        helpers: "Geant4BuildHelpers",
    ) -> tuple:
        """Set up the Geant4 model and return simulation components."""


def register_geant4_sim_type(name: str, strategy_cls: type) -> None:
    """Register a Geant4 simulation type by name."""
    _register(ENGINE, name, strategy_cls)


def get_registered_sim_types() -> dict[str, type]:
    """Return a copy of the Geant4 sim_type → strategy mapping."""
    return _get_registered(ENGINE)


class ShieldingAttenuationStrategy(Geant4SimStrategy):
    """Fixed-source particle transport through layered shielding.

    Geometry: concentric cylindrical shells (reactor_core type) or slab stack.
    Source: isotropic particle gun at specified energy.
    Scoring: per-layer energy deposition + transmission fraction.
    """

    sim_type: str = "shielding_attenuation"

    def build(
        self,
        geant4,
        solver_cfg: "Geant4Setting",
        geometry_cfg,
        materials_map: dict,
        helpers: "Geant4BuildHelpers",
    ) -> tuple:
        """Build Geant4 simulation components for shielding attenuation.

        Returns:
            Tuple of (detector, run_manager, physics_list, scoring) ready for execution.
        """
        from processforge.schemas.geant4.geant4_model import ReactorCoreGeometryConfig

        if not isinstance(geometry_cfg, ReactorCoreGeometryConfig):
            raise ValueError(
                f"shielding_attenuation requires ReactorCoreGeometryConfig, "
                f"got {type(geometry_cfg).__name__}"
            )

        # Build physics list
        physics_list_name = getattr(solver_cfg, "physics_list", "FTFP_BERT")
        physics_list = helpers.build_physics_list(geant4, physics_list_name)

        # Build geometry with concentric shells
        geometry_data = helpers.build_cylindrical_shells(
            geant4,
            geometry_cfg,
            materials_map,
        )

        # Build sensitive detector for scoring
        scoring = helpers.build_scoring(
            geant4,
            geometry_data["logical_volumes"],
            solver_cfg,
        )

        # Build primary generator action
        source_cfg = getattr(geometry_cfg, "source", None)
        generator = helpers.build_generator(
            geant4,
            source_cfg,
            geometry_cfg.source_point,
        )

        return {
            "physics_list": physics_list,
            "geometry": geometry_data,
            "scoring": scoring,
            "generator": generator,
        }


# Auto-register built-in strategies
register_geant4_sim_type("shielding_attenuation", ShieldingAttenuationStrategy)
