"""FESTIM simulation-type strategies."""
from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING

from processforge.providers._sim_strategy import (
    SimStrategy,
    get_registered_sim_types as _get_registered,
    register_sim_type as _register,
)
from processforge.schemas.festim.festim_model import FestimModel

if TYPE_CHECKING:
    from processforge.providers.festim.build_helpers import FestimBuildHelpers


ENGINE = "festim"


class FestimSimStrategy(SimStrategy):
    """Base class for a FESTIM simulation type."""

    sim_type: str = ""

    @abstractmethod
    def build(
        self,
        festim,
        festim_model: FestimModel,
        materials_map: dict,
        helpers: "FestimBuildHelpers",
    ) -> dict:
        """Set up the FESTIM model and return ``HydrogenTransportProblem`` kwargs."""


def register_festim_sim_type(name: str, strategy_cls: type) -> None:
    """Register a FESTIM simulation type by name."""
    _register(ENGINE, name, strategy_cls)


def get_registered_sim_types() -> dict[str, type]:
    """Return a copy of the FESTIM sim_type → strategy mapping."""
    return _get_registered(ENGINE)


class _TDSTrappingStrategy(FestimSimStrategy):
    """1-D hydrogen transport with trapping/detrapping (TDS-type simulation)."""

    sim_type: str = "hydrogen_transport_tds"

    def build(
        self,
        festim,
        festim_model: FestimModel,
        materials_map: dict,
        helpers: "FestimBuildHelpers",
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
