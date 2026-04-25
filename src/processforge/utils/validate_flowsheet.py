"""Flowsheet configuration validation: schema and graph-based connectivity checks."""
from dataclasses import dataclass
import json
from jsonschema import validate, ValidationError
from loguru import logger

from .._schema import load_flowsheet_schema

_VALID_DENSITY_UNITS = frozenset({
    "g/cm3", "g/cc", "kg/m3", "atom/b-cm", "atom/cm3", "sum", "macro"
})


@dataclass
class ConvergenceSignalConfig:
    """Configuration for provider-aware homotopy convergence signal."""

    signal_key: str
    target: float
    tolerance: float | None = None
    provider: str | None = None

    def to_jsonschema(self) -> dict:
        """Serialize to JSON schema for validation."""
        schema = {
            "type": "object",
            "required": ["signal_key", "target"],
            "properties": {
                "signal_key": {"type": "string"},
                "target": {"type": "number"},
                "tolerance": {"type": "number", "minimum": 0, "maximum": 1},
                "provider": {"type": "string", "enum": ["openmc", "festim"]},
            },
        }
        if self.tolerance is not None:
            schema["properties"]["tolerance"]["minimum"] = 0
            schema["properties"]["tolerance"]["maximum"] = 1
        return schema


def _validate_config_impl(config: dict, source_name: str) -> dict:
    """
    Validate a flowsheet configuration dict.

    Returns:
        The validated (and material-mix-resolved) config dict.

    Raises:
        ValidationError: On JSON schema violations.
        ValueError:      On semantic/connectivity violations.
    """
    schema = load_flowsheet_schema()
    validate(instance=config, schema=schema)
    logger.info(f"✅ Flowsheet '{source_name}' is valid.")
    check_stream_connectivity(config)
    _check_material_semantics(config)
    _check_provider_references(config)
    _check_provider_material_props(config)
    _check_openmc_unit_config(config)
    _resolve_material_mix_streams(config)
    _check_convergence_signal(config)
    return config


def validate_flowsheet(config_path):
    """
    Validate a flowsheet JSON file against the schema and connectivity rules.

    Args:
        config_path (str): Path to the flowsheet configuration JSON file.

    Returns:
        dict: The loaded and validated configuration dictionary.

    Raises:
        SystemExit: If schema validation or connectivity checks fail.
    """
    with open(config_path, "r", encoding="utf-8") as f:
        config = json.load(f)

    try:
        return _validate_config_impl(config, source_name=config_path)
    except ValidationError as e:
        logger.error(f"❌ Validation failed for '{config_path}':")
        logger.error(f"   → {e.message}")
        if e.path:
            logger.error(f"   Path: {' → '.join(map(str, e.path))}")
        raise SystemExit(1) from e
    except ValueError as e:
        logger.error(f"❌ Semantic validation failed for '{config_path}':")
        logger.error(f"   → {e}")
        raise SystemExit(1) from e


def validate_flowsheet_dict(config: dict, source_name: str = "<dict>") -> dict:
    """
    Validate an already-loaded flowsheet config dict.

    Args:
        config:      Flowsheet configuration dict.
        source_name: Optional label used in log/error messages.

    Returns:
        dict: The validated (and material-mix-resolved) config dict.

    Raises:
        SystemExit: If schema validation or connectivity checks fail.
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


def _get_unit_inlets(unit_config):
    """Return inlet stream name(s) as a list, handling single string or list."""
    inlet = unit_config.get("in")
    if isinstance(inlet, list):
        return inlet
    return [inlet] if inlet else []


def _get_unit_outlets(unit_config):
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


# Unit types that run standalone FEM/neutronics simulations — no mandatory
# inlet or outlet stream.  They are exempt from stream connectivity rules.
_SOLVER_UNIT_TYPES = frozenset({"SolverUnit"})


def _collect_produced_streams(config: dict) -> set:
    """Return all stream names produced by any unit (out, retentate_out, permeate_out)."""
    produced = set()
    for unit in config["units"].values():
        produced.update(_get_unit_outlets(unit))
    return produced


def check_stream_connectivity(config):
    """
    Validate stream connectivity in a flowsheet configuration.

    Streams are valid if declared in 'streams' (feeds) or produced as the
    'out' of any unit. Recycle streams are auto-detected from graph topology.

    SolverUnit types (FEM/neutronics) are exempt — they have no mandatory
    inlet or outlet stream.

    Raises:
        ValueError: If any connectivity check fails.

    Returns:
        bool: True if all checks pass.
    """
    defined_streams = set(config["streams"].keys())
    produced_streams = _collect_produced_streams(config)
    consumed_streams = set()

    for name, unit in config["units"].items():
        utype = unit.get("type", "")
        if utype in _SOLVER_UNIT_TYPES:
            # Standalone simulation units — no stream connectivity required
            continue
        inlets = _get_unit_inlets(unit)
        has_outlet = bool(_get_unit_outlets(unit))
        if not inlets or not has_outlet:
            raise ValueError(f"❌ Unit '{name}' missing 'in' or 'out' field.")
        consumed_streams.update(inlets)

    _check_inlet_sources(config, defined_streams, produced_streams)
    _check_unused_outlets(config, consumed_streams, defined_streams)
    _check_unreachable_units(config, defined_streams)
    _check_pipe_linkage(config)

    return True


def _check_inlet_sources(config, defined_streams, produced_streams):
    """Check every inlet has a valid source: feed stream or unit output."""
    valid_sources = defined_streams | produced_streams
    for name, unit in config["units"].items():
        for inlet in _get_unit_inlets(unit):
            if inlet not in valid_sources:
                raise ValueError(
                    f"❌ Unit '{name}' inlet '{inlet}' is not a feed stream "
                    f"or the output of any unit."
                )


def _check_unused_outlets(config, consumed_streams, defined_streams):
    """Warn about outlet streams not consumed by any downstream unit."""
    all_known_consumers = consumed_streams | defined_streams
    for name, unit in config["units"].items():
        for outlet in _get_unit_outlets(unit):
            if outlet not in all_known_consumers:
                logger.warning(
                    f"⚠️  Outlet '{outlet}' from unit '{name}' is not consumed "
                    f"by any downstream unit."
                )


def _check_unreachable_units(config, defined_streams):
    """
    Detect unreachable units by propagating reachability from feed streams.

    Runs N iterations (N = unit count) to cover all paths including cycles.
    SolverUnit types are excluded — they have no inlet stream dependency.
    """
    reachable = set(defined_streams)
    for _ in range(len(config["units"])):
        for ucfg in config["units"].values():
            if ucfg.get("type") in _SOLVER_UNIT_TYPES:
                continue
            if any(i in reachable for i in _get_unit_inlets(ucfg)):
                for outlet in _get_unit_outlets(ucfg):
                    reachable.add(outlet)

    unreachable = [
        u for u, cfg in config["units"].items()
        if cfg.get("type") not in _SOLVER_UNIT_TYPES
        and not any(o in reachable for o in _get_unit_outlets(cfg))
    ]
    if unreachable:
        raise ValueError(
            f"❌ Unreachable units detected: {', '.join(unreachable)}"
        )


_PIPE_LIKE_UNIT_TYPES = frozenset({"Pipes", "CSTR", "PFR", "IdealGasReactor"})

# Festim/solver unit types that produce streams valid as downstream inputs
_FESTIM_UNIT_TYPES = frozenset({"FestimMembrane", "SolverUnit"})


def _check_pipe_linkage(config):
    """
    Enforce that non-pipe units receive inputs only from feeds or pipe/reactor outputs.

    Reactor unit types (CSTR, PFR, IdealGasReactor) are treated as valid
    sources alongside Pipes units, since they similarly transform streams.
    Festim and SolverUnit types are also valid sources and are exempt from
    the linkage requirement themselves.
    """
    pipe_like_outputs = set()
    for u in config["units"].values():
        if u.get("type") in _PIPE_LIKE_UNIT_TYPES:
            pipe_like_outputs.update(_get_unit_outlets(u))

    festim_outputs = set()
    for u in config["units"].values():
        if u.get("type") in _FESTIM_UNIT_TYPES:
            festim_outputs.update(_get_unit_outlets(u))

    feed_streams = set(config["streams"].keys())
    valid_inputs = feed_streams | pipe_like_outputs | festim_outputs

    exempt_types = _PIPE_LIKE_UNIT_TYPES | _FESTIM_UNIT_TYPES
    for name, unit in config["units"].items():
        if unit.get("type") in exempt_types:
            continue
        for inlet in _get_unit_inlets(unit):
            if inlet not in valid_inputs:
                raise ValueError(
                    f"❌ Unit '{name}' input '{inlet}' must come from a "
                    f"feed stream or a Pipes / reactor / Festim unit output."
                )


def _check_provider_references(config: dict) -> None:
    """Validate that provider references in the flowsheet are consistent.

    Checks:
    1. Every unit's ``"provider"`` key names a declared entry in ``providers``.
    2. ``"default_provider"`` (if set) names a declared entry in ``providers``.
    """
    declared = set(config.get("providers", {}).keys())

    default = config.get("default_provider")
    if default is not None and default not in declared:
        raise ValueError(
            f"❌ 'default_provider' references '{default}', which is not "
            f"declared in the 'providers' block. "
            f"Declared providers: {sorted(declared)}"
        )

    for unit_name, unit_cfg in config.get("units", {}).items():
        ref = unit_cfg.get("provider")
        if ref is not None and ref not in declared:
            raise ValueError(
                f"❌ Unit '{unit_name}' references provider '{ref}', which is "
                f"not declared in the 'providers' block. "
                f"Declared providers: {sorted(declared)}"
            )


def _check_material_semantics(config):
    """
    Validate material definitions and enforce that every unit has a valid material.

    Checks:
    1. density_units values are valid OpenMC units (when materials section is present).
    2. id values are unique across all materials.
    3. material_mixes component names resolve to entries in materials.
    4. material_mixes fraction sum equals 1.0 (only when all fractions are provided).
    5. Stream z-keys resolve to defined materials (when materials section is present).
    6. Every unit must have a material field referencing a valid id.
    """
    materials = config.get("materials", {})
    valid_ids = {mat_def["id"] for mat_def in materials.values()}

    if materials:
        all_material_names = set(materials.keys())
        material_mixes = config.get("material_mixes", {})

        # 1. density_units validity — only checked when the field is present
        for mat_name, mat_def in materials.items():
            units = mat_def.get("density_units")
            if units is not None and units not in _VALID_DENSITY_UNITS:
                raise ValueError(
                    f"❌ Material '{mat_name}' has invalid density_units '{units}'. "
                    f"Must be one of: {sorted(_VALID_DENSITY_UNITS)}"
                )

        # 2. id uniqueness
        ids = [mat_def["id"] for mat_def in materials.values()]
        if len(ids) != len(set(ids)):
            seen, dupes = set(), set()
            for i in ids:
                if i in seen:
                    dupes.add(i)
                seen.add(i)
            raise ValueError(
                f"❌ Duplicate id values found: {sorted(dupes)}"
            )

        # 3. friendly_material_mix_id uniqueness
        mix_ids = [
            mix_def["friendly_material_mix_id"]
            for mix_def in material_mixes.values()
        ]
        if len(mix_ids) != len(set(mix_ids)):
            seen, dupes = set(), set()
            for i in mix_ids:
                if i in seen:
                    dupes.add(i)
                seen.add(i)
            raise ValueError(
                f"❌ Duplicate friendly_material_mix_id values found: {sorted(dupes)}"
            )

        # 4 & 5. Mix component name resolution and fraction sum
        for mix_name, mix_def in material_mixes.items():
            fractions = []
            for component in mix_def.get("components", []):
                comp_name = component["name"]
                if comp_name not in all_material_names:
                    raise ValueError(
                        f"❌ material_mixes['{mix_name}'] references '{comp_name}', "
                        f"which is not defined in materials."
                    )
                if "fraction" in component:
                    fractions.append(component["fraction"])

            components = mix_def.get("components", [])
            if len(fractions) == len(components) and fractions:
                frac_sum = sum(fractions)
                if abs(frac_sum - 1.0) > 1e-6:
                    raise ValueError(
                        f"❌ material_mixes['{mix_name}'] component fractions sum to "
                        f"{frac_sum:.8f}, expected 1.0."
                    )

        # 5. Stream z-keys resolve to materials
        for stream_name, stream_def in config.get("streams", {}).items():
            for comp_name in stream_def.get("z", {}):
                if comp_name not in all_material_names:
                    raise ValueError(
                        f"❌ Stream '{stream_name}' z-key '{comp_name}' is not defined in materials."
                    )

        # 5b. z and material_mix are mutually exclusive on a stream
        for stream_name, stream_def in config.get("streams", {}).items():
            if "z" in stream_def and "material_mix" in stream_def:
                raise ValueError(
                    f"❌ Stream '{stream_name}' defines both 'z' and 'material_mix'. "
                    f"Use 'material_mix' to reference a predefined mix, or 'z' for an "
                    f"explicit composition — not both."
                )

        # 6. Stream material_mix references must resolve to a valid friendly_material_mix_id
        valid_mix_ids = {
            mix_def["friendly_material_mix_id"]
            for mix_def in material_mixes.values()
        }
        for stream_name, stream_def in config.get("streams", {}).items():
            mix_ref = stream_def.get("material_mix")
            if mix_ref is not None and mix_ref not in valid_mix_ids:
                raise ValueError(
                    f"❌ Stream '{stream_name}' references material_mix id {mix_ref}, "
                    f"which does not match any friendly_material_mix_id in material_mixes."
                )

    # 7. Every unit must reference a valid id
    hint = " (no materials section defined)" if not materials else ""
    for unit_name, unit_def in config.get("units", {}).items():
        mat_id = unit_def.get("material")
        if mat_id not in valid_ids:
            raise ValueError(
                f"❌ Unit '{unit_name}' references material id {mat_id}, "
                f"which does not match any id in materials{hint}."
            )

    return True


def _check_provider_material_props(config: dict) -> None:
    """Delegate material validation to each provider's ``validate_material()`` classmethod.

    For every unit that references a declared provider, retrieves the provider
    class via the registry and calls ``cls.validate_material(mat_name, mat_def, unit_cfg)``.

    This function contains zero provider-specific logic — each provider
    enforces its own material requirements (Festim: D_0/E_D; OpenMC: nuclides; etc.).
    Adding a new provider never requires changes here.
    """
    from processforge.providers.manager import _maybe_import_provider
    from processforge.providers.registry import get_provider_class
    from processforge.types import MaterialDef, UnitConfig

    providers = config.get("providers", {})
    materials = config.get("materials", {})

    # Build id → (name, MaterialDef) lookup
    id_to_mat = {
        mat_dict["id"]: (mat_name, MaterialDef.from_dict(mat_dict))
        for mat_name, mat_dict in materials.items()
    }

    errors = []
    for unit_name, unit_dict in config.get("units", {}).items():
        provider_ref = unit_dict.get("provider")
        if provider_ref not in providers:
            continue
        ptype = providers[provider_ref].get("type", "")
        try:
            _maybe_import_provider(ptype)
            provider_cls = get_provider_class(ptype)
        except (ValueError, ImportError):
            # Provider not installed — skip (runtime error will surface later)
            continue

        mat_id = unit_dict.get("material")
        if mat_id not in id_to_mat:
            continue  # Already caught by _check_material_semantics

        mat_name, mat_def = id_to_mat[mat_id]
        unit_cfg = UnitConfig.from_dict(unit_dict)

        for msg in provider_cls.validate_material(mat_name, mat_def, unit_cfg):
            errors.append(f"❌ Unit '{unit_name}': {msg}")

    if errors:
        raise ValueError("\n".join(errors))


def _get_openmc_sim_type_registry() -> dict[str, type]:
    """Return a snapshot of registered OpenMC sim types from the provider module.

    Importing ``openmc_provider`` is safe at validation time because it does not
    import the ``openmc`` library at module level.
    """
    from processforge.providers.openmc_provider import get_registered_sim_types  # noqa: PLC0415
    return get_registered_sim_types()


def _check_openmc_unit_config(config: dict) -> None:
    """Validate OpenMC-specific ``SolverUnit`` fields when provider type is ``"openmc"``.

    Checks (only applied to units whose ``provider`` resolves to an OpenMC provider):

    1. ``sim_type`` is present.
    2. ``sim_type`` is registered in the provider's ``_SIM_TYPE_REGISTRY``
       (supports types added via :func:`~processforge.providers.openmc_provider.register_openmc_sim_type`).
    3. DAGMC sim types require ``dagmc_path`` in ``solver_config``.
    4. DAGMC sim types require ``source_box`` in ``solver_config``.
    5. ``inactive`` must be less than ``batches`` (when both are provided).
    6. Uses ``openmc_model.Material.model_validate()`` for deep Pydantic validation
       of each material referenced by an OpenMC unit.
    """
    providers = config.get("providers", {})
    openmc_provider_names = {
        name for name, pdef in providers.items()
        if pdef.get("type") == "openmc"
    }
    if not openmc_provider_names:
        return

    materials = config.get("materials", {})
    # Build id → (name, raw dict) lookup for material validation
    id_to_mat_raw = {
        mat_dict["id"]: (mat_name, mat_dict)
        for mat_name, mat_dict in materials.items()
    }

    registry = _get_openmc_sim_type_registry()
    errors = []

    for unit_name, unit_cfg in config.get("units", {}).items():
        if unit_cfg.get("provider") not in openmc_provider_names:
            continue

        sim_type = unit_cfg.get("sim_type")
        if not sim_type:
            errors.append(
                f"❌ Unit '{unit_name}' uses OpenMC provider but is missing 'sim_type'."
            )
            continue

        if sim_type not in registry:
            errors.append(
                f"❌ Unit '{unit_name}' has unknown OpenMC sim_type '{sim_type}'. "
                f"Registered types: {sorted(registry)}"
            )
            continue

        sc = unit_cfg.get("solver_config") or {}

        if sim_type.endswith("_dagmc") and not sc.get("dagmc_path"):
            errors.append(
                f"❌ Unit '{unit_name}' sim_type='{sim_type}' requires "
                f"'dagmc_path' in solver_config."
            )

        if sim_type.endswith("_dagmc") and not sc.get("source_box"):
            errors.append(
                f"❌ Unit '{unit_name}' sim_type='{sim_type}' requires "
                f"'source_box' in solver_config."
            )

        batches = sc.get("batches", 0)
        inactive = sc.get("inactive", 0)
        if batches > 0 and inactive >= batches:
            errors.append(
                f"❌ Unit '{unit_name}': solver_config 'inactive' ({inactive}) "
                f"must be less than 'batches' ({batches})."
            )

        # Deep Pydantic validation using the now-self-contained openmc_model.Material
        mat_id = unit_cfg.get("material")
        if mat_id in id_to_mat_raw:
            mat_name, mat_raw = id_to_mat_raw[mat_id]
            try:
                from processforge.schemas.openmc.openmc_model import Material as OpenMCMaterial
                # Build a dict compatible with the Pydantic Material model.
                # Pass through only raw values present in the source material so
                # missing required fields are surfaced by Pydantic validation
                # instead of being masked by local defaults.
                pydantic_dict = {
                    "name": mat_name,
                    "id": mat_raw["id"],
                    "nuclides": [
                        {
                            "nuclide": n["name"],
                            "fraction": n["percent"],
                            **(
                                {"percent_type": n["percent_type"]}
                                if "percent_type" in n else {}
                            ),
                        }
                        for n in mat_raw.get("nuclides", [])
                    ],
                    "temperature": mat_raw.get("temperature"),
                    "elements": mat_raw.get("elements"),
                }
                if "density" in mat_raw:
                    pydantic_dict["density"] = mat_raw["density"]
                if "density_units" in mat_raw:
                    pydantic_dict["density_units"] = mat_raw["density_units"]
                OpenMCMaterial.model_validate(pydantic_dict)
            except Exception as exc:  # noqa: BLE001
                errors.append(
                    f"❌ Unit '{unit_name}' material '{mat_name}' failed "
                    f"OpenMC schema validation: {exc}"
                )

    if errors:
        raise ValueError("\n".join(errors))


def _resolve_material_mix_streams(config: dict) -> dict:
    """Expand material_mix references on streams into z composition dicts.

    For every stream that carries a ``material_mix`` integer
    (a ``friendly_material_mix_id``), this function looks up the corresponding
    mix definition and writes its component fractions into ``stream["z"]``.
    This expansion happens after validation so that the solver receives a fully
    populated ``z`` dict regardless of whether the author wrote ``z`` explicitly
    or used a ``material_mix`` reference.

    The mutual-exclusivity rule (``z`` and ``material_mix`` cannot coexist) is
    enforced earlier in ``_check_material_semantics``, so by the time this
    function runs every stream has at most one of the two.

    Args:
        config: Validated flowsheet configuration dict (mutated in place).

    Returns:
        The same config dict with ``z`` populated on any stream that had
        ``material_mix``.
    """
    material_mixes = config.get("material_mixes", {})
    if not material_mixes:
        return config

    id_to_mix = {
        mix_def["friendly_material_mix_id"]: mix_def
        for mix_def in material_mixes.values()
    }

    for stream_def in config.get("streams", {}).values():
        mix_ref = stream_def.get("material_mix")
        if mix_ref is not None:
            mix = id_to_mix[mix_ref]  # existence already validated
            stream_def["z"] = {
                comp["name"]: comp["fraction"]
                for comp in mix["components"]
                if "fraction" in comp
            }

    return config


def _check_convergence_signal(config: dict) -> None:
    """Validate the convergence_signal configuration in simulation block."""
    from jsonschema import validate as jsonschema_validate

    sim = config.get("simulation", {})
    conv_signal = sim.get("convergence_signal")
    if not conv_signal:
        return

    schema = ConvergenceSignalConfig(
        signal_key="",  # placeholder for schema generation
        target=0.0,
    ).to_jsonschema()
    jsonschema_validate(instance=conv_signal, schema=schema)

    if conv_signal.get("provider") is None:
        providers_in_use = {
            unit_cfg.get("provider")
            for unit_cfg in config.get("units", {}).values()
            if unit_cfg.get("provider") in ("openmc", "festim")
        }
        if not providers_in_use:
            raise ValueError(
                "❌ convergence_signal expects at least one SolverUnit with provider "
                "'openmc' or 'festim' in the flowsheet when provider is not specified."
            )
