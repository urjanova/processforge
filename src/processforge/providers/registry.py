"""Provider registry — maps provider type names to their classes.

Mirrors the ``_BACKENDS`` pattern in ``processforge.eo.solver``.

Optional providers (Cantera, Modelica) self-register by calling
``register_provider()`` at module import time.  They are imported lazily
in ``manager.py`` only when the user declares them in the flowsheet JSON,
so installing the extras is never required for base functionality.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .base import AbstractProvider

# Registry is seeded with the built-in CoolProp provider at bottom of file.
_PROVIDERS: dict[str, type] = {}


def get_provider_class(provider_type: str) -> type:
    """Return the provider class for the given type string.

    Raises:
        ValueError: If the type is not registered.
    """
    cls = _PROVIDERS.get(provider_type)
    if cls is None:
        raise ValueError(
            f"Unknown provider type '{provider_type}'. "
            f"Registered types: {sorted(_PROVIDERS)}"
        )
    return cls


def register_provider(name: str, cls: type) -> None:
    """Register a provider class under *name*.

    Called by each optional provider module so the registry stays current
    without hard imports at the top level.
    """
    _PROVIDERS[name] = cls


# Seed the registry with the always-available CoolProp provider.
# Import is deferred to a function-level import to avoid circular deps.
def _seed_registry() -> None:
    from .coolprop_provider import CoolPropProvider
    register_provider("coolprop", CoolPropProvider)


_seed_registry()
