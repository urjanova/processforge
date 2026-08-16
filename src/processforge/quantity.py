"""Unit-aware quantities for Processforge simulation outputs.

A single :class:`Quantity` type is used everywhere engine outputs are
standardized (OpenMC, FESTIM, CoolProp, Cantera).  It wraps a `pint`
:meth:`UnitRegistry` so values carry real units with validation and
conversion, but stores magnitude as a plain ``list``/``float`` plus a unit
string so it round-trips cleanly through JSON / pydantic / Zarr.

The registry is exposed module-level so every provider tags outputs through
the *same* instance — there is exactly one source of units in the process.
"""
from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, field_serializer, field_validator
from pint import UnitRegistry

#: The single process-wide unit registry.
ureg = UnitRegistry()

# A few domain aliases that make engine outputs read naturally.
ureg.define("keff = []") if "keff" not in ureg else None
ureg.define("percent = 0.01 = %") if "percent" not in ureg else None


def _coerce_magnitude(v):
    """Normalize numpy / scalar magnitudes into JSON-friendly Python types."""
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.floating):
        return float(v)
    if isinstance(v, np.integer):
        return float(v)
    if isinstance(v, (list, tuple)):
        return [_coerce_magnitude(x) for x in v]
    return v


class Quantity(BaseModel):
    """A unit-aware value with optional uncertainty.

    ``value`` is a scalar or (possibly nested) list of floats; ``unit`` is a
    ``pint``-parsable unit string (``""`` means dimensionless).  ``std_dev``
    carries the same shape as ``value`` when present.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    value: float | list[float]
    unit: str = ""
    std_dev: float | list[float] | None = None

    @field_validator("value", "std_dev", mode="before")
    @classmethod
    def _coerce(cls, v):
        if v is None:
            return v
        return _coerce_magnitude(v)

    @field_serializer("value", "std_dev")
    def _serialize(self, v):
        return v

    # ------------------------------------------------------------------
    # pint interop
    # ------------------------------------------------------------------
    def to_pint(self):
        """Return the equivalent ``pint.Quantity`` (raises if dimensionless)."""
        if not self.unit:
            return np.asarray(self.value, dtype=float)
        return ureg.Quantity(np.asarray(self.value, dtype=float), self.unit)

    @classmethod
    def from_pint(cls, q, std_dev=None) -> "Quantity":
        """Build a :class:`Quantity` from a ``pint.Quantity``."""
        mag = q.magnitude
        unit = str(q.units)
        sd = None
        if std_dev is not None:
            sd = std_dev.magnitude if hasattr(std_dev, "magnitude") else std_dev
        return cls(value=mag, unit=unit, std_dev=sd)

    def to(self, units: str) -> "Quantity":
        """Convert to *units*, returning a new :class:`Quantity`."""
        if not self.unit:
            raise ValueError("Cannot convert a dimensionless Quantity.")
        conv = ureg.Quantity(np.asarray(self.value, dtype=float), self.unit).to(units)
        return Quantity(value=conv.magnitude.tolist(), unit=str(conv.units))

    def is_dimensionless(self) -> bool:
        if not self.unit:
            return True
        try:
            return ureg.Quantity(1, self.unit).dimensionless
        except Exception:
            return False

    def magnitude(self):
        return np.asarray(self.value, dtype=float)


def quantity(value, unit: str | None = None, std_dev=None) -> Quantity:
    """Convenience constructor for a :class:`Quantity`."""
    return Quantity(value=value, unit=unit or "", std_dev=std_dev)
