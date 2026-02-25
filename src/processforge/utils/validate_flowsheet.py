"""Flowsheet configuration validation: schema and graph-based connectivity checks."""
import json
from jsonschema import validate, ValidationError
from loguru import logger

from .._schema import load_flowsheet_schema


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
    schema = load_flowsheet_schema()

    with open(config_path, "r", encoding="utf-8") as f:
        config = json.load(f)

    try:
        validate(instance=config, schema=schema)
        logger.info(f"✅ Flowsheet '{config_path}' is valid.")
        check_stream_connectivity(config)
        return config
    except ValidationError as e:
        logger.error(f"❌ Validation failed for '{config_path}':")
        logger.error(f"   → {e.message}")
        if e.path:
            logger.error(f"   Path: {' → '.join(map(str, e.path))}")
        raise SystemExit(1) from e


def _get_unit_inlets(unit_config):
    """Return inlet stream name(s) as a list, handling single string or list."""
    inlet = unit_config.get("in")
    if isinstance(inlet, list):
        return inlet
    return [inlet] if inlet else []


def check_stream_connectivity(config):
    """
    Validate stream connectivity in a flowsheet configuration.

    Streams are valid if declared in 'streams' (feeds) or produced as the
    'out' of any unit. Recycle streams are auto-detected from graph topology.

    Raises:
        ValueError: If any connectivity check fails.

    Returns:
        bool: True if all checks pass.
    """
    defined_streams = set(config["streams"].keys())
    produced_streams = {unit["out"] for unit in config["units"].values()}
    consumed_streams = set()

    for name, unit in config["units"].items():
        inlets = _get_unit_inlets(unit)
        if not inlets or not unit.get("out"):
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
        outlet = unit["out"]
        if outlet not in all_known_consumers:
            logger.warning(
                f"⚠️  Outlet '{outlet}' from unit '{name}' is not consumed "
                f"by any downstream unit."
            )


def _check_unreachable_units(config, defined_streams):
    """
    Detect unreachable units by propagating reachability from feed streams.

    Runs N iterations (N = unit count) to cover all paths including cycles.
    """
    reachable = set(defined_streams)
    for _ in range(len(config["units"])):
        for ucfg in config["units"].values():
            if any(i in reachable for i in _get_unit_inlets(ucfg)):
                reachable.add(ucfg["out"])

    unreachable = [
        u for u, cfg in config["units"].items()
        if cfg["out"] not in reachable
    ]
    if unreachable:
        raise ValueError(
            f"❌ Unreachable units detected: {', '.join(unreachable)}"
        )


def _check_pipe_linkage(config):
    """
    Enforce that non-pipe units receive inputs only from feeds or pipe outputs.
    """
    pipe_outputs = {
        u["out"] for u in config["units"].values() if u["type"] == "Pipes"
    }
    feed_streams = set(config["streams"].keys())
    valid_inputs = feed_streams | pipe_outputs

    for name, unit in config["units"].items():
        if unit["type"] == "Pipes":
            continue
        for inlet in _get_unit_inlets(unit):
            if inlet not in valid_inputs:
                raise ValueError(
                    f"❌ Unit '{name}' input '{inlet}' must come from a "
                    f"feed stream or a Pipes unit output."
                )
