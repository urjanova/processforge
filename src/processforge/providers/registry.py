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

# Canonical catalog of all supported providers.
# ``list_providers()`` reads this to report what exists without importing anything.
_PROVIDER_CATALOG: dict[str, dict[str, str | None]] = {
    "coolprop": {
        "module": "processforge.providers.coolprop_provider",
        "class": "CoolPropProvider",
        "optional_dep": None,
        "description": "Thermodynamic properties via CoolProp (built-in)",
    },
    "cantera": {
        "module": "processforge.providers.cantera_provider",
        "class": "CanteraProvider",
        "optional_dep": "cantera",
        "description": "Thermochemistry and reactor kinetics via Cantera",
    },
    "modelica": {
        "module": "processforge.providers.modelica_provider",
        "class": "ModelicaProvider",
        "optional_dep": "modelica",
        "description": "FMU-based simulation via OpenModelica",
    },
    "openmc": {
        "module": "processforge.providers.openmc_provider",
        "class": "OpenMCProvider",
        "optional_dep": "openmc",
        "description": "Neutronics simulation via OpenMC",
    },
}


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


def list_providers() -> list[dict[str, object]]:
    """Return metadata for every supported provider.

    Each entry contains:

    * ``type`` — provider type string (e.g. ``"cantera"``)
    * ``description`` — one-line summary
    * ``optional_dep`` — pip extras name if the provider requires an optional
      dependency, or ``None`` for built-in providers
    * ``installed`` — whether the provider module is importable right now
      (i.e. the optional dependency is installed)
    * ``registered`` — whether the provider class has already been loaded
      into the runtime registry

    This function does **not** import any provider modules, so it is safe
    to call at startup regardless of which extras are installed.
    """
    results: list[dict[str, object]] = []
    for type_name, info in _PROVIDER_CATALOG.items():
        module_name = info["module"]
        try:
            importlib.util.find_spec(module_name)
            installed = True
        except (ModuleNotFoundError, ValueError):
            installed = False

        results.append(
            {
                "type": type_name,
                "description": info["description"],
                "optional_dep": info["optional_dep"],
                "installed": installed,
                "registered": type_name in _PROVIDERS,
            }
        )
    return results


# Seed the registry with the always-available CoolProp provider.
# Import is deferred to a function-level import to avoid circular deps.
def _seed_registry() -> None:
    from .coolprop_provider import CoolPropProvider
    register_provider("coolprop", CoolPropProvider)


_seed_registry()
