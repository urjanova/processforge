"""Provider registry — maps provider type names to their classes.

Mirrors the ``_BACKENDS`` pattern in ``processforge.eo.solver``.

Optional providers (Cantera, Modelica) self-register by calling
``register_provider()`` at module import time.  They are imported lazily
only when the user declares them in the flowsheet JSON, so installing the
extras is never required for base functionality.
"""
from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .base import AbstractProvider

# Registry is seeded with the built-in CoolProp provider at bottom of file.
_PROVIDERS: dict[str, type] = {}

# Types that are always available (no optional dependency).
_BUILTIN_TYPES: frozenset[str] = frozenset({"coolprop"})


def get_provider_class(provider_type: str) -> type:
    """Return the provider class for the given type string.

    If the type is not yet registered, the registry attempts a lazy import
    of ``processforge.providers.{provider_type}_provider`` (convention-based)
    before failing.  This keeps import-time cost zero for providers the
    flowsheet never references.

    Raises:
        ValueError: If the type is not registered and cannot be imported.
    """
    cls = _PROVIDERS.get(provider_type)
    if cls is not None:
        return cls

    # Not yet registered — try the convention-based lazy import.
    if provider_type not in _BUILTIN_TYPES:
        module_name = f"processforge.providers.{provider_type}_provider"
        try:
            importlib.import_module(module_name)
        except ModuleNotFoundError as exc:
            raise ValueError(
                f"No provider module found for type '{provider_type}'. "
                f"Expected module '{module_name}' to exist and self-register."
            ) from exc

        cls = _PROVIDERS.get(provider_type)
        if cls is not None:
            return cls

    raise ValueError(
        f"Unknown provider type '{provider_type}'. "
        f"Registered types: {sorted(_PROVIDERS)}"
    )


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
