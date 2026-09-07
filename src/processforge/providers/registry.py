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

from pydantic import BaseModel, ConfigDict, Field

if TYPE_CHECKING:
    from .base import AbstractProvider


class ProviderCatalogEntry(BaseModel):
    """Validated metadata entry for a provider in the canonical catalog."""

    module: str
    class_name: str = Field(alias="class")
    optional_dep: str | None
    description: str
    docker_image: str | None = None
    default_port: int | None = None

    model_config = ConfigDict(populate_by_name=True)


# Runtime registry of loaded provider classes.
_PROVIDERS: dict[str, type] = {}

# Types that are always available (no optional dependency).
_BUILTIN_TYPES: frozenset[str] = frozenset({"coolprop"})

# Whether the built-in CoolProp provider has been seeded into `_PROVIDERS`.
# Seeding is deferred until first use so importing `registry` never imports
# provider backends unless they are actually needed.
_SEEDED: bool = False

# Canonical catalog of all supported providers.
# ``list_providers()`` reads this to report what exists without importing anything.
_PROVIDER_CATALOG: dict[str, ProviderCatalogEntry] = {
    "coolprop": ProviderCatalogEntry(
        module="processforge.providers.coolprop_provider",
        class_name="CoolPropProvider",
        optional_dep=None,
        description="Thermodynamic properties via CoolProp (built-in)",
    ),
    "cantera": ProviderCatalogEntry(
        module="processforge.providers.cantera_provider",
        class_name="CanteraProvider",
        optional_dep="cantera",
        description="Thermochemistry and reactor kinetics via Cantera",
    ),
    "modelica": ProviderCatalogEntry(
        module="processforge.providers.modelica_provider",
        class_name="ModelicaProvider",
        optional_dep="modelica",
        description="FMU-based simulation via OpenModelica",
    ),
    "openmc": ProviderCatalogEntry(
        module="processforge.providers.openmc",
        class_name="OpenMCProvider",
        optional_dep="openmc",
        description="Neutronics simulation via OpenMC",
        docker_image="ghcr.io/urjanova/processforge-openmc:latest",
        default_port=9001,
    ),
    "festim": ProviderCatalogEntry(
        module="processforge.providers.festim",
        class_name="FestimProvider",
        optional_dep=None,
        description="Hydrogen transport FEM via FESTIM (Docker service)",
        docker_image="ghcr.io/urjanova/processforge-festim:latest",
        default_port=9002,
    ),
}


def _ensure_seeded() -> None:
    """Seed the registry with the always-available CoolProp provider.

    This is a one-time, lazy operation. It is safe to call repeatedly.
    """
    global _SEEDED
    if _SEEDED:
        return
    # Set the flag BEFORE registering to prevent recursion: register_provider
    # also calls _ensure_seeded().
    _SEEDED = True
    from .coolprop_provider import CoolPropProvider

    _PROVIDERS["coolprop"] = CoolPropProvider


def get_provider_class(provider_type: str) -> type:
    """Return the provider class for the given type string.

    If the type is not yet registered, the registry attempts a lazy import
    of ``processforge.providers.{provider_type}_provider`` (convention-based)
    before failing.  This keeps import-time cost zero for providers the
    flowsheet never references.

    Raises:
        ValueError: If the type is not registered and cannot be imported.
    """
    _ensure_seeded()
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
    _ensure_seeded()
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
    _ensure_seeded()
    results: list[dict[str, object]] = []
    for type_name, info in _PROVIDER_CATALOG.items():
        module_name = info.module
        try:
            importlib.util.find_spec(module_name)
            installed = True
        except (ModuleNotFoundError, ValueError):
            installed = False

        results.append(
            {
                "type": type_name,
                "description": info.description,
                "optional_dep": info.optional_dep,
                "docker_image": info.docker_image,
                "default_port": info.default_port,
                "installed": installed,
                "registered": type_name in _PROVIDERS,
            }
        )
    return results


def get_provider_docker_image(provider_type: str) -> str | None:
    """Return the Docker image for a containerized provider, or None."""
    info = _PROVIDER_CATALOG.get(provider_type)
    if info is None:
        raise ValueError(f"Unknown provider type '{provider_type}'")
    return info.docker_image


def get_provider_default_port(provider_type: str) -> int | None:
    """Return the default port for a containerized provider, or None."""
    info = _PROVIDER_CATALOG.get(provider_type)
    if info is None:
        raise ValueError(f"Unknown provider type '{provider_type}'")
    return info.default_port


def is_containerized(provider_type: str) -> bool:
    """Return True if the provider runs in a Docker container.

    Returns False for unknown provider types rather than raising, so callers
    can cleanly fall through to the local/import path.
    """
    try:
        return get_provider_docker_image(provider_type) is not None
    except ValueError:
        return False
