"""Shared simulation-type strategy registry for engine providers.

Engine providers such as OpenMC and FESTIM support multiple ``sim_type`` values.
Each simulation type is implemented as a strategy class that builds the
engine-specific model for that type.  This module provides the shared
infrastructure so new engines can reuse the same registration pattern.

Example::

    from processforge.providers._sim_strategy import SimStrategy, register_sim_type

    class MyEngineStrategy(SimStrategy):
        config_model = MyGeometryConfig

        def build(self, engine_module, solver_cfg, geometry_cfg, materials_map):
            ...
            return built_model

    register_sim_type("my_engine", "my_sim", MyEngineStrategy)
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class SimStrategy(ABC):
    """Base class for an engine-specific simulation type.

    Subclasses declare ``config_model`` — the Pydantic model validating the
    per-strategy ``geometry_config`` block of the flowsheet — and implement
    ``build`` to construct the engine model.
    """

    #: Pydantic model validating the unit's ``geometry_config`` block.  ``None``
    #: means the strategy takes no geometry config.
    config_model: type | None = None

    @abstractmethod
    def build(self, **kwargs: Any) -> Any:
        """Construct the engine model and return engine-specific objects."""


# Engine name -> sim_type name -> strategy class
_SIM_TYPE_REGISTRIES: dict[str, dict[str, type]] = {}


def register_sim_type(engine: str, name: str, strategy_cls: type) -> None:
    """Register a simulation type for a given engine.

    Args:
        engine: Engine/provider name (e.g. ``"openmc"``, ``"festim"``).
        name: The ``sim_type`` string used in the flowsheet JSON.
        strategy_cls: A subclass of :class:`SimStrategy`.
    """
    _SIM_TYPE_REGISTRIES.setdefault(engine, {})[name] = strategy_cls


def get_registered_sim_types(engine: str) -> dict[str, type]:
    """Return a copy of the sim_type → strategy mapping for *engine*.

    Returns:
        A dict mapping ``sim_type`` strings to strategy classes.
    """
    return dict(_SIM_TYPE_REGISTRIES.get(engine, {}))


def get_strategy(engine: str, sim_type: str) -> type | None:
    """Return the strategy class for (*engine*, *sim_type*), or ``None``."""
    return _SIM_TYPE_REGISTRIES.get(engine, {}).get(sim_type)
