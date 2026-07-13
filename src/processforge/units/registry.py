"""Unit type registry for processforge.

Unit classes self-register via :func:`register_unit` at module import time.
Both ``Flowsheet`` and ``EOFlowsheet`` use the registry as the single source
of truth instead of maintaining hardcoded type maps.
"""
from __future__ import annotations

_UNIT_TYPES: dict[str, type] = {}


def register_unit(name: str, cls: type) -> None:
    """Register a unit class under the given *name*."""
    _UNIT_TYPES[name] = cls


def get_unit_class(name: str) -> type | None:
    """Look up a unit class by name.  Returns ``None`` if not found."""
    return _UNIT_TYPES.get(name)


def get_all_unit_types() -> dict[str, type]:
    """Return a shallow copy of the full registry."""
    return dict(_UNIT_TYPES)
