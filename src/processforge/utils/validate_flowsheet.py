"""Flowsheet configuration validation: schema and graph-based connectivity checks.

Two API layers:

* :func:`_validate_config_impl` — programmatic entry point. Raises
  ``jsonschema.ValidationError`` on schema violations and ``ValueError`` on
  semantic/connectivity violations. The runner wraps these in its own
  :class:`~processforge.runner.FlowsheetValidationError`.
* :func:`validate_flowsheet` / :func:`validate_flowsheet_dict` — CLI-facing
  wrappers. These are the only places in this module that raise
  ``SystemExit``; all underlying checks raise plain exceptions so they stay
  library-friendly.

Provider-specific rules (material requirements, ``sim_type``/``solver_config``
checks, …) live on the provider classes via the ``validate_material`` and
``validate_unit_config`` hooks — this module contains zero provider-specific
logic, so adding a new provider never requires changes here.
"""
from __future__ import annotations

import copy
import json
import re

from jsonschema import ValidationError, validate
from loguru import logger

from .._schema import load_flowsheet_schema


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Unit types that run standalone FEM/neutronics simulations — no mandatory
# inlet or outlet stream.  They are exempt from stream connectivity rules and
# their outputs count as valid sources for downstream units.
_SOLVER_UNIT_TYPES = frozenset({"SolverUnit"})

# Unit types treated as valid stream sources for downstream units, alongside
# SolverUnit outputs (they transform streams like Pipes units do).
_PIPE_LIKE_UNIT_TYPES = frozenset({"Pipes", "CSTR", "PFR", "IdealGasReactor"})

# Unit types whose inlets must come from a feed stream or a Pipes / reactor /
# SolverUnit output (enforced by _check_pipe_linkage).
_LINKAGE_EXEMPT_TYPES = _PIPE_LIKE_UNIT_TYPES | _SOLVER_UNIT_TYPES

_COUPLE_REF_RE = re.compile(r"^[A-Za-z0-9_]+\.[A-Za-z0-9_]+$")
_COUPLE_REDUCE = {"mean", "sum", "max", "min"}

# Base schema for ``simulation.convergence_signal``. The ``provider`` property
# is added per-config by _check_convergence_signal, restricted to the provider
# names actually declared in the flowsheet.
_CONVERGENCE_SIGNAL_SCHEMA: dict = {
    "type": "object",
    "required": ["signal_key", "target"],
    "properties": {
        "signal_key": {"type": "string"},
        "target": {"type": "number"},
        "tolerance": {"type": "number", "minimum": 0, "maximum": 1},
    },
}


def _get_valid_density_units() -> frozenset:
    """Return the set of valid OpenMC material density units.

    The canonical list lives in the OpenMC provider module so it is defined in
    exactly one place; it is imported lazily (importing the provider module is
    safe at validation time — it does not import the ``openmc`` library). A
    local fallback keeps the check working if the provider module is unavailable.
    """
    try:
        from processforge.providers.openmc.provider import _VALID_DENSITY_UNITS
    except ImportError:
        return frozenset({
            "g/cm3", "g/cc", "kg/m3", "atom/b-cm", "atom/cm3", "sum", "macro"
        })
    return _VALID_DENSITY_UNITS


def _duplicates(values) -> list:
    """Return the sorted list of values that appear more than once."""
    seen: set = set()
    dupes: set = set()
    for value in values:
        if value in seen:
            dupes.add(value)
        seen.add(value)
    return sorted(dupes)


# ---------------------------------------------------------------------------
# Stream accessors
# ---------------------------------------------------------------------------

def _get_unit_inlets(unit_config: dict) -> list[str]:
    """Return inlet stream name(s) as a list, handling single string or list."""
    inlet = unit_config.get("in")
    if isinstance(inlet, list):
        return inlet
    return [inlet] if inlet else []


def _get_unit_outlets(unit_config: dict) -> list[str]:
    """Return outlet stream name(s) as a list, handling string or list for 'out'."""
    outlets = []
    out = unit_config.get("out")
    if isinstance(out, list):
        outlets.extend(out)
    elif out is not None:
        outlets.append(out)
    for key in ("retentate_out", "permeate_out"):
        val = unit_config.get(key)
        if val is not None:
            outlets.append(val)
    return outlets


def _collect_produced_streams(config: dict) -> set:
    """Return all stream names produced by any unit (out, retentate_out, permeate_out)."""
    produced = set()
    for unit in config.get("units", {}).values():
        produced.update(_get_unit_outlets(unit))
    return produced


# ---------------------------------------------------------------------------
# Core validation pipeline
# ---------------------------------------------------------------------------

def _validate_config_impl(config: dict, source_name: str) -> dict:
    """
    Validate a flowsheet configuration dict.

    The checks run as an explicit pipeline; ordering matters — e.g.
    ``_check_material_semantics`` must run before
    ``_resolve_material_mix_streams`` (mutual-exclusivity of ``z`` and
    ``material_mix``), and provider-reference checks must run before the
    provider dispatch resolves provider classes.

    Returns:
        The validated (and material-mix-resolved) config dict. The input dict
        is mutated in place by ``_resolve_material_mix_streams``; callers that
        must not observe mutation should pass a copy (as
        :func:`validate_flowsheet_dict` does).

    Raises:
        ValidationError: On JSON schema violations.
        ValueError:      On semantic/connectivity violations.
    """
    schema = load_flowsheet_schema()
    validate(instance=config, schema=schema)

    for check in (
        check_stream_connectivity,
        _check_material_semantics,
        _check_provider_references,
        _check_provider_unit_configs,
        _check_coupling_refs,
    ):
        check(config)

    config = _resolve_material_mix_streams(config)
    _check_convergence_signal(config)

    logger.info(f"✅ Flowsheet '{source_name}' is valid.")
    return config


def _validate_or_exit(config: dict, source_name: str) -> dict:
    """Run validation, logging failures and converting them to ``SystemExit(1)``.

    The only place in this module where ``SystemExit`` is raised — keeps the
    underlying checks library-friendly (plain exceptions).
    """
    try:
        return _validate_config_impl(config, source_name=source_name)
    except ValidationError as e:
        logger.error(f"❌ Validation failed for '{source_name}':")
        logger.error(f"   → {e.message}")
        if e.path:
            logger.error(f"   Path: {' → '.join(map(str, e.path))}")
        raise SystemExit(1) from e
    except ValueError as e:
        logger.error(f"❌ Semantic validation failed for '{source_name}':")
        logger.error(f"   → {e}")
        raise SystemExit(1) from e


def validate_flowsheet(config_path: str) -> dict:
    """
    Validate a flowsheet JSON file against the schema and connectivity rules.

    Args:
        config_path: Path to the flowsheet configuration JSON file.

    Returns:
        dict: The loaded and validated configuration dictionary.

    Raises:
        SystemExit: If the file cannot be read/parsed or validation fails.
    """
    try:
        with open(config_path, "r", encoding="utf-8") as f:
            config = json.load(f)
    except OSError as exc:
        logger.error(f"❌ Could not read flowsheet '{config_path}': {exc}")
        raise SystemExit(1) from exc
    except json.JSONDecodeError as exc:
        logger.error(f"❌ Flowsheet '{config_path}' is not valid JSON: {exc}")
        raise SystemExit(1) from exc

    return _validate_or_exit(config, source_name=config_path)


def validate_flowsheet_dict(config: dict, source_name: str = "<dict>") -> dict:
    """
    Validate an already-loaded flowsheet config dict.

    The input dict is deep-copied first, so validation-time expansion
    (``material_mix`` → ``z``) does not leak back into the caller's dict; the
    validated copy is returned.

    Args:
        config:      Flowsheet configuration dict.
        source_name: Optional label used in log/error messages.

    Returns:
        dict: The validated (and material-mix-resolved) config copy.

    Raises:
        SystemExit: If schema validation or connectivity checks fail.
    """
    return _validate_or_exit(copy.deepcopy(config), source_name=source_name)


# ---------------------------------------------------------------------------
# Stream connectivity
# ---------------------------------------------------------------------------

def check_stream_connectivity(config: dict) -> None:
    """
    Validate stream connectivity in a flowsheet configuration.

    Streams are valid if declared in 'streams' (feeds) or produced as the
    'out' of any unit. Recycle streams are auto-detected from graph topology.

    SolverUnit types (FEM/neutronics) are exempt — they have no mandatory
    inlet or outlet stream.

    Raises:
        ValueError: If any connectivity check fails (all problems reported).
    """
    defined_streams = set(config["streams"].keys())
    produced_streams = _collect_produced_streams(config)
    consumed_streams: set = set()
    errors: list[str] = []

    for name, unit in config["units"].items():
        if unit.get("type", "") in _SOLVER_UNIT_TYPES:
            # Standalone simulation units — no stream connectivity required
            continue
        inlets = _get_unit_inlets(unit)
        has_outlet = bool(_get_unit_outlets(unit))
        if not inlets or not has_outlet:
            errors.append(f"❌ Unit '{name}' missing 'in' or 'out' field.")
        consumed_streams.update(inlets)

    errors.extend(_check_inlet_sources(config, defined_streams, produced_streams))
    _check_unused_outlets(config, consumed_streams, defined_streams)
    errors.extend(_check_unreachable_units(config, defined_streams))
    errors.extend(_check_pipe_linkage(config))

    if errors:
        raise ValueError("\n".join(errors))


def _check_inlet_sources(
    config: dict, defined_streams: set, produced_streams: set
) -> list[str]:
    """Check every inlet has a valid source: feed stream or unit output."""
    valid_sources = defined_streams | produced_streams
    errors: list[str] = []
    for name, unit in config["units"].items():
        for inlet in _get_unit_inlets(unit):
            if inlet not in valid_sources:
                errors.append(
                    f"❌ Unit '{name}' inlet '{inlet}' is not a feed stream "
                    f"or the output of any unit."
                )
    return errors


def _check_unused_outlets(
    config: dict, consumed_streams: set, defined_streams: set
) -> None:
    """Warn about outlet streams not consumed by any downstream unit."""
    all_known_consumers = consumed_streams | defined_streams
    for name, unit in config["units"].items():
        for outlet in _get_unit_outlets(unit):
            if outlet not in all_known_consumers:
                logger.warning(
                    f"⚠️  Outlet '{outlet}' from unit '{name}' is not consumed "
                    f"by any downstream unit."
                )


def _check_unreachable_units(config: dict, defined_streams: set) -> list[str]:
    """
    Detect units disconnected from any valid boundary.

    Connectivity is checked on the undirected unit graph, so a feed-less
    closed loop is not treated as orphaned merely because it has no entry
    point. A connected group of units is valid when it contains at least one
    boundary:

    - a unit consuming a declared feed stream;
    - a SolverUnit (no inlet stream dependency);
    - in dynamic mode, a Tank carrying initial state.

    SolverUnit types are excluded from the orphan check — they have no inlet
    stream dependency.
    """
    units = config["units"]
    sim_mode = config.get("simulation", {}).get("mode", "steady")

    # stream -> producing unit
    producers: dict = {}
    for unit_name, ucfg in units.items():
        for outlet in _get_unit_outlets(ucfg):
            producers[outlet] = unit_name

    # Undirected unit adjacency: producer <-> consumer.
    adj: dict = {name: set() for name in units}
    for unit_name, ucfg in units.items():
        for inlet in _get_unit_inlets(ucfg):
            producer = producers.get(inlet)
            if producer is not None and producer != unit_name:
                adj[unit_name].add(producer)
                adj[producer].add(unit_name)

    # Valid boundaries seed the traversal.
    reachable = {
        name for name, ucfg in units.items()
        if ucfg.get("type") in _SOLVER_UNIT_TYPES
    }
    for name, ucfg in units.items():
        if any(i in defined_streams for i in _get_unit_inlets(ucfg)):
            reachable.add(name)
    if sim_mode == "dynamic":
        for name, ucfg in units.items():
            if ucfg.get("type") == "Tank" and ucfg.get("initial_n") is not None:
                reachable.add(name)

    stack = list(reachable)
    while stack:
        current = stack.pop()
        for neighbour in adj[current]:
            if neighbour not in reachable:
                reachable.add(neighbour)
                stack.append(neighbour)

    unreachable = [
        name for name, ucfg in units.items()
        if ucfg.get("type") not in _SOLVER_UNIT_TYPES
        and name not in reachable
    ]
    if unreachable:
        return [
            (
                f"❌ Unreachable units detected: {', '.join(unreachable)}. "
                f"Each connected group of units needs a feed stream or, in dynamic "
                f"mode, a Tank with initial conditions."
            )
        ]
    return []


def _check_pipe_linkage(config: dict) -> list[str]:
    """
    Enforce that non-pipe units receive inputs only from feeds or pipe/reactor outputs.

    Reactor unit types (CSTR, PFR, IdealGasReactor) are treated as valid
    sources alongside Pipes units, since they similarly transform streams.
    SolverUnit types are also valid sources and are exempt from the linkage
    requirement themselves.
    """
    valid_inputs = set(config["streams"].keys())
    for u in config["units"].values():
        if u.get("type") in _LINKAGE_EXEMPT_TYPES:
            valid_inputs.update(_get_unit_outlets(u))

    errors: list[str] = []
    for name, unit in config["units"].items():
        if unit.get("type") in _LINKAGE_EXEMPT_TYPES:
            continue
        for inlet in _get_unit_inlets(unit):
            if inlet not in valid_inputs:
                errors.append(
                    f"❌ Unit '{name}' input '{inlet}' must come from a "
                    f"feed stream or a Pipes / reactor / SolverUnit output."
                )
    return errors


# ---------------------------------------------------------------------------
# Provider references
# ---------------------------------------------------------------------------

def _require_declared_ref(
    errors: list[str],
    declared: set,
    owner_kind: str,
    owner_name: str,
    field: str,
    ref: str | None,
    required: bool,
) -> None:
    """Append a validation error if a provider reference is missing or unknown.

    Args:
        errors:     Accumulator to append messages to.
        declared:   Declared provider names (keys of the ``providers`` block).
        owner_kind: ``"Unit"`` or ``"Stream"`` — used in messages.
        owner_name: Name of the unit/stream owning the reference.
        field:      Reference field name (``"provider"`` or ``"thermo_provider"``).
        ref:        The referenced provider name (may be ``None``).
        required:   Whether an absent reference is itself an error.
    """
    if ref is None:
        if required:
            errors.append(
                f"❌ {owner_kind} '{owner_name}' omits '{field}' but 'coolprop' is not "
                f"declared in 'providers'. Add it to 'providers' or set '{field}' "
                f"explicitly on {owner_kind.lower()} '{owner_name}'."
            )
    elif ref not in declared:
        errors.append(
            f"❌ {owner_kind} '{owner_name}' references {field} '{ref}', which is "
            f"not declared in 'providers'. Declared providers: {sorted(declared)}"
        )


def _check_provider_references(config: dict) -> None:
    """Validate that provider references in the flowsheet are consistent.

    Checks:
    1. Every unit's ``"provider"`` key names a declared entry in ``providers``.
    2. Every unit's ``"thermo_provider"`` key names a declared entry in ``providers``.
    3. Every stream's ``"thermo_provider"`` key names a declared entry in ``providers``.
    4. ``"default_provider"`` (if set) names a declared entry in ``providers``.

    ``coolprop`` is the built-in default thermo provider: when it is not
    declared, every unit/stream must set its provider references explicitly;
    when it is declared, the references are optional but must resolve if present.
    """
    declared = set(config.get("providers", {}).keys())
    coolprop_declared = "coolprop" in declared
    required = not coolprop_declared
    errors: list[str] = []

    default = config.get("default_provider")
    if default is not None and default not in declared:
        errors.append(
            f"❌ 'default_provider' references '{default}', which is not "
            f"declared in the 'providers' block. "
            f"Declared providers: {sorted(declared)}"
        )

    for unit_name, unit_cfg in config.get("units", {}).items():
        _require_declared_ref(
            errors, declared, "Unit", unit_name,
            "provider", unit_cfg.get("provider"), required,
        )
        _require_declared_ref(
            errors, declared, "Unit", unit_name,
            "thermo_provider", unit_cfg.get("thermo_provider"), required,
        )
    for stream_name, stream_cfg in config.get("streams", {}).items():
        _require_declared_ref(
            errors, declared, "Stream", stream_name,
            "thermo_provider", stream_cfg.get("thermo_provider"), required,
        )

    if errors:
        raise ValueError("\n".join(errors))


# ---------------------------------------------------------------------------
# Material semantics
# ---------------------------------------------------------------------------

def _check_material_semantics(config: dict) -> None:
    """
    Validate material definitions and enforce that every unit has a valid material.

    Checks (only when a materials section is present, unless noted):
    1. density_units values are valid OpenMC units.
    2. material id values are unique.
    3. material_mix id values are unique.
    4. material_mixes component names resolve to entries in materials.
    5. material_mixes fraction sum equals 1.0 (only when all fractions are provided).
    6. Stream z-keys resolve to defined materials.
    7. Stream z and material_mix are mutually exclusive.
    8. Stream material_mix references resolve to a valid material_mix id.
    9. Every unit references a valid material id (always checked; SolverUnit
       without a material field is exempt).
    """
    materials = config.get("materials", {})
    valid_ids = {mat_def["id"] for mat_def in materials.values()}
    errors: list[str] = []

    if materials:
        all_material_names = set(materials.keys())
        material_mixes = config.get("material_mixes", {})
        valid_density_units = _get_valid_density_units()

        # 1. density_units validity — only checked when the field is present
        for mat_name, mat_def in materials.items():
            units = mat_def.get("density_units")
            if units is not None and units not in valid_density_units:
                errors.append(
                    f"❌ Material '{mat_name}' has invalid density_units '{units}'. "
                    f"Must be one of: {sorted(valid_density_units)}"
                )

        # 2. material id uniqueness
        ids = [mat_def["id"] for mat_def in materials.values()]
        dupes = _duplicates(ids)
        if dupes:
            errors.append(f"❌ Duplicate id values found: {dupes}")

        # 3. material_mix id uniqueness
        mix_ids = [mix_def["id"] for mix_def in material_mixes.values()]
        dupes = _duplicates(mix_ids)
        if dupes:
            errors.append(f"❌ Duplicate material_mix id values found: {dupes}")

        # 4 & 5. Mix component name resolution and fraction sum
        for mix_name, mix_def in material_mixes.items():
            components = mix_def.get("components", [])
            fractions = []
            for component in components:
                comp_name = component["name"]
                if comp_name not in all_material_names:
                    errors.append(
                        f"❌ material_mixes['{mix_name}'] references '{comp_name}', "
                        f"which is not defined in materials."
                    )
                if "fraction" in component:
                    fractions.append(component["fraction"])

            if len(fractions) == len(components) and fractions:
                frac_sum = sum(fractions)
                if abs(frac_sum - 1.0) > 1e-6:
                    errors.append(
                        f"❌ material_mixes['{mix_name}'] component fractions sum to "
                        f"{frac_sum:.8f}, expected 1.0."
                    )

        # 6. Stream z-keys resolve to materials
        for stream_name, stream_def in config.get("streams", {}).items():
            for comp_name in stream_def.get("z", {}):
                if comp_name not in all_material_names:
                    errors.append(
                        f"❌ Stream '{stream_name}' z-key '{comp_name}' is not "
                        f"defined in materials."
                    )

        # 7. z and material_mix are mutually exclusive on a stream
        for stream_name, stream_def in config.get("streams", {}).items():
            if "z" in stream_def and "material_mix" in stream_def:
                errors.append(
                    f"❌ Stream '{stream_name}' defines both 'z' and 'material_mix'. "
                    f"Use 'material_mix' to reference a predefined mix, or 'z' for an "
                    f"explicit composition — not both."
                )

        # 8. Stream material_mix references must resolve to a valid material_mix id
        valid_mix_ids = {mix_def["id"] for mix_def in material_mixes.values()}
        for stream_name, stream_def in config.get("streams", {}).items():
            mix_ref = stream_def.get("material_mix")
            if mix_ref is not None and mix_ref not in valid_mix_ids:
                errors.append(
                    f"❌ Stream '{stream_name}' references material_mix id {mix_ref}, "
                    f"which does not match any id in material_mixes."
                )

    # 9. Every unit must reference a valid id (SolverUnit without material field is exempt)
    hint = " (no materials section defined)" if not materials else ""
    for unit_name, unit_def in config.get("units", {}).items():
        mat_id = unit_def.get("material")
        if mat_id is None and unit_def.get("type") in _SOLVER_UNIT_TYPES:
            continue  # SolverUnit without material field — validated by provider
        if mat_id not in valid_ids:
            errors.append(
                f"❌ Unit '{unit_name}' references material id {mat_id}, "
                f"which does not match any id in materials{hint}."
            )

    if errors:
        raise ValueError("\n".join(errors))


# ---------------------------------------------------------------------------
# Provider-owned unit validation (dispatch)
# ---------------------------------------------------------------------------

def _check_provider_unit_configs(config: dict) -> None:
    """Run provider-owned validation hooks for every unit backed by a provider.

    For every unit that references a declared provider, resolves the provider
    class via the registry and calls:

    * ``cls.validate_material(mat_name, mat_def, unit_cfg)`` — provider-specific
      material requirements (OpenMC: nuclides and density; FESTIM: diffusivity …).
    * ``cls.validate_unit_config(unit_dict, config)`` — provider-specific unit
      configuration rules (OpenMC: sim_type/solver_config/geometry_config;
      FESTIM: the full problem model …).

    This function contains zero provider-specific logic. Adding a new provider
    never requires changes here; providers that are not installed are skipped
    (the runtime error will surface later).
    """
    from processforge.providers.registry import get_provider_class
    from processforge.types import MaterialDef, UnitConfig

    providers = config.get("providers", {})
    materials = config.get("materials", {})

    # Build id → (name, raw dict) lookup for material validation
    id_to_mat = {
        mat_dict["id"]: (mat_name, mat_dict)
        for mat_name, mat_dict in materials.items()
    }

    errors: list[str] = []
    for unit_name, unit_dict in config.get("units", {}).items():
        provider_ref = unit_dict.get("provider")
        if provider_ref not in providers:
            continue
        ptype = providers[provider_ref].get("type", "")
        try:
            provider_cls = get_provider_class(ptype)
        except (ValueError, ImportError):
            # Provider not installed — skip (runtime error will surface later)
            continue

        # 1. Provider-owned material validation.
        mat_id = unit_dict.get("material")
        if mat_id in id_to_mat:
            mat_name, mat_dict = id_to_mat[mat_id]
            unit_cfg = UnitConfig.from_dict(unit_dict)
            for msg in provider_cls.validate_material(
                mat_name, MaterialDef.from_dict(mat_dict), unit_cfg
            ):
                errors.append(f"❌ Unit '{unit_name}': {msg}")

        # 2. Provider-owned unit-config validation.
        for msg in provider_cls.validate_unit_config(unit_dict, config):
            errors.append(f"❌ Unit '{unit_name}': {msg}")

    if errors:
        raise ValueError("\n".join(errors))


# ---------------------------------------------------------------------------
# Coupling references (cross-engine parameter coupling)
# ---------------------------------------------------------------------------

def _check_coupling_refs(config: dict) -> None:
    """Validate ``inputs`` coupling declarations on every unit.

    Each entry maps a dotted config path to ``{"ref": "<unit>.<field>", ...}``.
    We check the ref format and that its source (a unit or a stream) exists.
    Cycles between coupled units are allowed — the flowsheet driver iterates
    until they converge.
    """
    units = config.get("units", {})
    # Declared feed streams and streams produced by any unit (out,
    # retentate_out, permeate_out) are valid coupling sources.
    streams = set(config.get("streams", {}).keys()) | _collect_produced_streams(config)
    unit_names = set(units.keys())
    errors: list[str] = []

    for unit_name, unit_cfg in units.items():
        inputs = unit_cfg.get("inputs")
        if not inputs:
            continue
        if not isinstance(inputs, dict):
            errors.append(
                f"❌ Unit '{unit_name}': 'inputs' must be an object of "
                f"dotted-path -> coupling spec."
            )
            continue
        for path, spec in inputs.items():
            if not isinstance(spec, dict) or "ref" not in spec:
                errors.append(
                    f"❌ Unit '{unit_name}': coupling '{path}' must be an object "
                    f"with a 'ref' key."
                )
                continue
            ref = spec["ref"]
            if not isinstance(ref, str) or not _COUPLE_REF_RE.match(ref):
                errors.append(
                    f"❌ Unit '{unit_name}': coupling '{path}' has invalid ref "
                    f"'{ref}' (expected '<unit>.<field>')."
                )
                continue
            src = ref.split(".")[0]
            if src not in unit_names and src not in streams:
                errors.append(
                    f"❌ Unit '{unit_name}': coupling '{path}' references unknown "
                    f"source '{src}' (not a declared unit or stream)."
                )
            reduce_mode = spec.get("reduce")
            if reduce_mode is not None and reduce_mode not in _COUPLE_REDUCE:
                errors.append(
                    f"❌ Unit '{unit_name}': coupling '{path}' has invalid "
                    f"reduce '{reduce_mode}' (expected one of {sorted(_COUPLE_REDUCE)})."
                )
            as_unit = spec.get("as_unit")
            if as_unit is not None and not isinstance(as_unit, str):
                errors.append(
                    f"❌ Unit '{unit_name}': coupling '{path}' 'as_unit' must be a string."
                )

    if errors:
        raise ValueError("\n".join(errors))


# ---------------------------------------------------------------------------
# Material-mix resolution
# ---------------------------------------------------------------------------

def _resolve_material_mix_streams(config: dict) -> dict:
    """Expand material_mix references on streams into z composition dicts.

    For every stream that carries a ``material_mix`` integer
    (a ``material_mix_id``), this function looks up the corresponding
    mix definition and writes its component fractions into ``stream["z"]``.
    This expansion happens after validation so that the solver receives a fully
    populated ``z`` dict regardless of whether the author wrote ``z`` explicitly
    or used a ``material_mix`` reference.

    The mutual-exclusivity rule (``z`` and ``material_mix`` cannot coexist) is
    enforced earlier in ``_check_material_semantics``, so by the time this
    function runs every stream has at most one of the two.

    Args:
        config: Validated flowsheet configuration dict (mutated in place and
            returned, so callers can write ``config = _resolve(config)``).

    Returns:
        The same config dict with ``z`` populated on any stream that had
        ``material_mix``.
    """
    material_mixes = config.get("material_mixes", {})
    if not material_mixes:
        return config

    id_to_mix = {
        mix_def["id"]: mix_def
        for mix_def in material_mixes.values()
    }

    for stream_name, stream_def in config.get("streams", {}).items():
        mix_ref = stream_def.get("material_mix")
        if mix_ref is None:
            continue
        mix = id_to_mix.get(mix_ref)  # existence already validated
        if mix is None:  # defensive — keep the failure message actionable
            raise ValueError(
                f"❌ Stream '{stream_name}' references material_mix id {mix_ref}, "
                f"which does not match any id in material_mixes."
            )
        stream_def["z"] = {
            comp["name"]: comp["fraction"]
            for comp in mix["components"]
            if "fraction" in comp
        }

    return config


# ---------------------------------------------------------------------------
# Convergence signal
# ---------------------------------------------------------------------------

def _check_convergence_signal(config: dict) -> None:
    """Validate the convergence_signal configuration in the simulation block.

    The JSON schema is built per-config so the ``provider`` field is restricted
    to the provider names actually declared in the flowsheet — the runtime
    homotopy driver compares ``convergence_signal.provider`` against declared
    provider *names*.

    When ``provider`` is omitted, at least one unit must reference a provider
    whose type is ``openmc`` (the only engine shipping a default convergence
    signal).
    """
    sim = config.get("simulation", {})
    conv_signal = sim.get("convergence_signal")
    if not conv_signal:
        return

    declared = set(config.get("providers", {}).keys())
    schema = copy.deepcopy(_CONVERGENCE_SIGNAL_SCHEMA)
    schema["properties"]["provider"] = {"type": "string"}

    validate(instance=conv_signal, schema=schema)

    provider_ref = conv_signal.get("provider")
    if provider_ref is not None:
        if provider_ref not in declared:
            raise ValueError(
                f"❌ convergence_signal references provider '{provider_ref}', which "
                f"is not declared in the 'providers' block. "
                f"Declared providers: {sorted(declared)}"
            )
        return

    providers_cfg = config.get("providers", {})
    openmc_provider_names = {
        name for name, pdef in providers_cfg.items()
        if pdef.get("type") == "openmc"
    }
    providers_in_use = {
        unit_cfg.get("provider")
        for unit_cfg in config.get("units", {}).values()
        if unit_cfg.get("provider") in openmc_provider_names
    }
    if not providers_in_use:
        raise ValueError(
            "❌ convergence_signal expects at least one SolverUnit with an "
            "OpenMC provider in the flowsheet when provider is not specified."
        )
