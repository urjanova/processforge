"""Stream unit consistency checker using pint (optional dependency)."""
from __future__ import annotations
from dataclasses import dataclass


@dataclass
class UnitMismatch:
    stream_name: str
    property_name: str
    declared_unit: str
    compatible: bool
    message: str


def check_unit_consistency(config: dict) -> list[UnitMismatch]:
    """
    Check dimensional validity of stream _units annotations.

    If pint is not installed, returns an empty list silently.
    If no streams carry _units annotations, returns an empty list.

    Tiers:
      1. Same dimensionality (kg/s ↔ lb/hr): info only, compatible=True
      2. Molar/mass ambiguity (kg/s ↔ mol/s): warning, compatible=True
      3. Incompatible dimensions (K ↔ Pa): error, compatible=False

    Args:
        config: Flowsheet config dict, streams may carry a '_units' sub-dict.

    Returns:
        List of UnitMismatch objects (empty if pint unavailable or no annotations).
    """
    try:
        import pint
    except ImportError:
        return []

    ureg = pint.UnitRegistry()

    # Expected dimensionalities for each stream property
    expected_dimensionality: dict[str, str] = {
        "T":        "[temperature]",
        "P":        "[mass] / [length] / [time] ** 2",
        "flowrate": None,  # accepts [mass]/[time] or [substance]/[time]
    }

    mismatches: list[UnitMismatch] = []

    for stream_name, stream_cfg in config.get("streams", {}).items():
        units_ann = stream_cfg.get("_units", {})
        if not units_ann:
            continue

        for prop, unit_str in units_ann.items():
            try:
                unit = ureg.Unit(unit_str)
            except pint.errors.UndefinedUnitError:
                mismatches.append(UnitMismatch(
                    stream_name=stream_name,
                    property_name=prop,
                    declared_unit=unit_str,
                    compatible=False,
                    message=f"Unrecognised unit '{unit_str}'",
                ))
                continue

            dim = unit.dimensionality

            if prop == "T":
                if "[temperature]" not in str(dim):
                    mismatches.append(UnitMismatch(
                        stream_name=stream_name,
                        property_name=prop,
                        declared_unit=unit_str,
                        compatible=False,
                        message=(
                            f"Stream '{stream_name}'.T annotated as '{unit_str}' "
                            f"but expected a temperature unit (e.g. K, degC)."
                        ),
                    ))

            elif prop == "P":
                expected = ureg.Unit("Pa").dimensionality
                if dim != expected:
                    mismatches.append(UnitMismatch(
                        stream_name=stream_name,
                        property_name=prop,
                        declared_unit=unit_str,
                        compatible=False,
                        message=(
                            f"Stream '{stream_name}'.P annotated as '{unit_str}' "
                            f"but expected a pressure unit (e.g. Pa, bar, psi)."
                        ),
                    ))

            elif prop == "flowrate":
                mass_dim = ureg.Unit("kg/s").dimensionality
                molar_dim = ureg.Unit("mol/s").dimensionality
                if dim == mass_dim:
                    mismatches.append(UnitMismatch(
                        stream_name=stream_name,
                        property_name=prop,
                        declared_unit=unit_str,
                        compatible=True,
                        message=(
                            f"Stream '{stream_name}'.flowrate annotated as '{unit_str}' "
                            f"(mass flow). Solver expects mol/s — conversion requires MW data."
                        ),
                    ))
                elif dim == molar_dim:
                    pass  # mol/s or equivalent — fine
                else:
                    mismatches.append(UnitMismatch(
                        stream_name=stream_name,
                        property_name=prop,
                        declared_unit=unit_str,
                        compatible=False,
                        message=(
                            f"Stream '{stream_name}'.flowrate annotated as '{unit_str}' "
                            f"which is neither a mass flow nor a molar flow unit."
                        ),
                    ))

    return mismatches


def strip_units_annotations(config: dict) -> None:
    """Remove _units keys from all stream dicts in-place."""
    for stream_cfg in config.get("streams", {}).values():
        stream_cfg.pop("_units", None)
