"""Degrees of Freedom (DOF) analysis for processforge flowsheets."""
from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class UnitDOFResult:
    unit_name: str
    unit_type: str
    n_variables: int
    n_equations: int
    dof: int
    status: str          # "determined" | "under-specified" | "over-specified"
    issues: list[str] = field(default_factory=list)


@dataclass
class SystemDOFReport:
    n_components: int
    component_names: list[str]
    feed_stream_specs: int
    total_variables: int
    total_equations: int
    system_dof: int
    per_unit: list[UnitDOFResult] = field(default_factory=list)


# Required parameters per unit type. Each entry is a list of options:
#   - a string means the key must be present
#   - a tuple means at least one of the keys must be present
_UNIT_REQUIRED: dict[str, list] = {
    "Pump":            ["deltaP"],
    "Valve":           ["pressure_ratio"],
    "Strainer":        ["deltaP"],
    "Pipes":           ["delta_p", "diameter"],
    "Heater":          [("duty", "outlet_T")],    # at least one of these
    "Flash":           ["P"],
    "Tank":            ["outlet_flow", "P"],
    "CSTR":            ["residence_time"],
    "PFR":             ["residence_time"],
    "IdealGasReactor": ["residence_time"],
}

# Units that produce two outlet streams (affects variables count)
_TWO_OUTLET_UNITS = frozenset({"Flash"})


def _check_required(unit_cfg: dict, unit_type: str, sim_mode: str) -> list[str]:
    """Return a list of issue strings for any missing required params."""
    issues = []
    required = list(_UNIT_REQUIRED.get(unit_type, []))

    # Dynamic mode: Tank needs initial conditions too
    if unit_type == "Tank" and sim_mode == "dynamic":
        required = required + ["initial_n", "initial_T"]

    for req in required:
        if isinstance(req, tuple):
            # At least one of these keys must be present
            if not any(k in unit_cfg for k in req):
                issues.append(
                    f"missing at least one of: {', '.join(req)}"
                )
        else:
            if req not in unit_cfg:
                issues.append(f"missing required param '{req}'")

    # Cantera-backed reactors: provider must be declared
    if unit_type in ("CSTR", "PFR", "IdealGasReactor"):
        provider_ref = unit_cfg.get("provider")
        if provider_ref is not None:
            # The provider declaration check is handled by validate_flowsheet_dict;
            # here we just note it for the DOF report.
            pass

    return issues


def analyze_dof(config: dict) -> SystemDOFReport:
    """
    Compute per-unit and system-level DOF for a validated flowsheet config.

    Args:
        config: Validated flowsheet config dict (post _resolve_material_mix_streams).

    Returns:
        SystemDOFReport with per-unit results and system totals.
    """
    sim_mode = config.get("simulation", {}).get("mode", "steady")

    # Determine component set from feed streams
    component_names: list[str] = []
    for stream_cfg in config.get("streams", {}).values():
        for comp in stream_cfg.get("z", {}).keys():
            if comp not in component_names:
                component_names.append(comp)
    component_names = sorted(component_names)
    n_c = len(component_names)
    vars_per_stream = 3 + n_c  # T, P, F, z_i × n_c

    feed_stream_specs = len(config.get("streams", {})) * vars_per_stream

    per_unit: list[UnitDOFResult] = []
    total_variables = 0
    total_equations = 0

    for unit_name, unit_cfg in config.get("units", {}).items():
        unit_type = unit_cfg.get("type", "")
        n_outlets = 2 if unit_type in _TWO_OUTLET_UNITS else 1
        variables = n_outlets * vars_per_stream
        equations = n_outlets * vars_per_stream  # each unit fully determines its outlets

        issues = _check_required(unit_cfg, unit_type, sim_mode)

        # Each missing required param = one DOF penalty
        dof_penalty = len(issues)
        dof = dof_penalty  # 0 = determined, >0 = under-specified

        if dof == 0:
            status = "determined"
        elif dof > 0:
            status = "under-specified"
        else:
            status = "over-specified"

        per_unit.append(UnitDOFResult(
            unit_name=unit_name,
            unit_type=unit_type,
            n_variables=variables,
            n_equations=equations,
            dof=dof,
            status=status,
            issues=issues,
        ))

        total_variables += variables
        total_equations += equations

    # System DOF = total missing required parameters across all units.
    # Feed streams are boundary conditions (fully given), not variables.
    # Each unit fully determines its outlets when all required params are present,
    # so the system is determined iff every unit has DOF=0.
    system_dof = sum(r.dof for r in per_unit)

    return SystemDOFReport(
        n_components=n_c,
        component_names=component_names,
        feed_stream_specs=feed_stream_specs,
        total_variables=total_variables,
        total_equations=total_equations,
        system_dof=system_dof,
        per_unit=per_unit,
    )
