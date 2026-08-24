"""Multi-physics coupling resolver for Processforge flowsheets.

This module implements *parameter-reference coupling* between heterogeneous
solver units (OpenMC neutronics, FESTIM thermal, CoolProp thermo, …).  Each
solver unit may declare an ``inputs`` block mapping a dotted configuration
path to a reference to another unit's (or stream's) output field::

    "inputs": {
        "solver_config.heat_source": {"ref": "openmc_solver.power", "as_unit": "W"},
        "solver_config.temperature": {"ref": "thermal.temperature", "reduce": "mean", "as_unit": "K"}
    }

A single shared :class:`ParameterStore` collects every unit/stream output as a
unit-bearing :class:`~processforge.quantity.Quantity` (the same ``Quantity``
type every engine already emits — see Phase 5A).  The :class:`CouplingResolver`
then, for each consumer unit and each declared input:

1. looks the referenced ``<unit>.<field>`` quantity up in the store,
2. reduces array-valued fields to a scalar (``mean`` / ``sum`` / ``max`` /
   ``min`` — the natural choice for uniform coupling),
3. converts units through the process-wide ``pint`` registry, and
4. injects the bare magnitude into the consumer's config at the dotted path.

Because solver units run in topological order, a single forward pass resolves
acyclic couplings; cyclic couplings (e.g. neutronics power ↔ thermal
temperature) converge via the iterative loop in :mod:`processforge.flowsheet`.
"""

from __future__ import annotations

import re
from typing import Any, Optional

from loguru import logger

from .quantity import Quantity, ureg

__all__ = [
    "CouplingError",
    "ParameterStore",
    "CouplingResolver",
    "deep_merge",
]


class CouplingError(ValueError):
    """Raised when a coupling reference cannot be resolved."""


_REF_RE = re.compile(r"^(?P<src>[A-Za-z0-9_]+)\.(?P<field>[A-Za-z0-9_]+)$")


def _reduce_value(value: Any, mode: str) -> float:
    """Collapse an (optionally array) value to a scalar for uniform coupling."""
    if isinstance(value, (list, tuple)):
        arr = [float(x) for x in value]
        if not arr:
            return 0.0
        if mode == "sum":
            return sum(arr)
        if mode == "max":
            return max(arr)
        if mode == "min":
            return min(arr)
        return sum(arr) / len(arr)
    if value is None:
        return 0.0
    return float(value)


def _convert(scalar: float, from_unit: str, to_unit: Optional[str]) -> float:
    """Convert *scalar* from *from_unit* to *to_unit* via the shared registry."""
    if not to_unit or not from_unit:
        return float(scalar)
    try:
        return float(ureg.Quantity(float(scalar), from_unit).to(to_unit).magnitude)
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            f"Coupling unit conversion {from_unit} -> {to_unit} failed ({exc}); "
            f"using raw magnitude {scalar}."
        )
        return float(scalar)


def _deep_set(d: dict, dotted_path: str, value: Any) -> None:
    """Set *value* at *dotted_path* inside nested dict *d*, creating maps."""
    parts = dotted_path.split(".")
    cur = d
    for p in parts[:-1]:
        nxt = cur.get(p)
        if not isinstance(nxt, dict):
            nxt = {}
            cur[p] = nxt
        cur = nxt
    cur[parts[-1]] = value


def deep_merge(base: dict, overrides: dict) -> dict:
    """Recursively merge *overrides* into a copy of *base*."""
    out = dict(base)
    for k, v in (overrides or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


class ParameterStore:
    """Registry of unit/stream output quantities, keyed by ``<src>.<field>``."""

    def __init__(self) -> None:
        self._quantities: dict[str, Quantity] = {}

    def register(self, src: str, field: str, quantity: Quantity) -> None:
        self._quantities[f"{src}.{field}"] = quantity

    def register_unit(self, unit_name: str, engine_output: Any) -> None:
        """Register every field of an :class:`EngineOutput`."""
        for field in getattr(engine_output, "fields", []) or []:
            self.register(unit_name, field.name, field.quantity)

    def register_stream(
        self,
        stream_name: str,
        stream_dict: dict,
        thermo_fields: Optional[dict] = None,
    ) -> None:
        """Register stream scalars (T in K, P in Pa) and any thermo fields."""
        if "T" in stream_dict:
            self.register(
                stream_name, "T", Quantity(value=float(stream_dict["T"]), unit="K")
            )
        if "P" in stream_dict:
            self.register(
                stream_name, "P", Quantity(value=float(stream_dict["P"]), unit="Pa")
            )
        for fname, q in (thermo_fields or {}).items():
            if isinstance(q, Quantity):
                self.register(stream_name, fname, q)

    def get(self, ref: str) -> Optional[Quantity]:
        return self._quantities.get(ref)

    def __contains__(self, ref: str) -> bool:
        return ref in self._quantities


class CouplingResolver:
    """Resolve a unit's ``inputs`` block against a :class:`ParameterStore`."""

    def __init__(self, store: ParameterStore) -> None:
        self.store = store

    def resolve(self, spec: dict) -> float:
        ref = spec.get("ref") if isinstance(spec, dict) else None
        if not ref:
            raise CouplingError("coupling spec missing 'ref'")
        if not _REF_RE.match(ref):
            raise CouplingError(
                f"invalid coupling ref '{ref}' (expected '<unit>.<field>')"
            )
        quantity = self.store.get(ref)
        if quantity is None:
            raise CouplingError(f"coupling ref '{ref}' not found in store")
        reduced = _reduce_value(quantity.value, spec.get("reduce", "mean"))
        return _convert(reduced, quantity.unit, spec.get("as_unit"))

    def build_overrides(self, inputs: dict):
        """Return ``(nested_overrides, injected_scalars)`` for a unit.

        *nested_overrides* mirrors the dotted input paths and is deep-merged
        into the unit config before the provider runs.  *injected_scalars*
        maps each local path to the resolved float, used for convergence
        checking by the flowsheet driver.
        """
        overrides: dict = {}
        injected: dict = {}
        for path, spec in (inputs or {}).items():
            value = self.resolve(spec)
            _deep_set(overrides, path, value)
            injected[path] = value
        return overrides, injected
