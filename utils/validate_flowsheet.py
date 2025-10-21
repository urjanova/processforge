import json
import os
from jsonschema import validate, ValidationError
from loguru import logger


def validate_flowsheet(config_path):
    """
    Validates a flowsheet configuration file against a predefined JSON schema and performs additional connectivity checks.
    This function loads the flowsheet schema from a fixed path relative to the script's location, validates the provided
    configuration file against it using the jsonschema library, and then runs custom connectivity checks on the streams.
    If validation succeeds, it logs a success message and returns the loaded configuration dictionary. If validation fails,
    it logs an error message with details and raises SystemExit to terminate the program.
    Args:
        config_path (str): The file path to the flowsheet configuration JSON file to be validated.
    Returns:
        dict: The loaded and validated configuration dictionary if validation and checks pass.
    Raises:
        SystemExit: If the configuration fails schema validation or connectivity checks, with an exit code of 1.
    Note:
        The schema file is expected to be located at '../schemas/flowsheet_schema.json' relative to the directory of this script.
        Additional connectivity checks are performed via the check_stream_connectivity function.
    """

    pwd = os.path.dirname(os.path.abspath(__file__))
    # go one level up
    schema_directory = os.path.dirname(pwd)
    schema_path = os.path.join(schema_directory, "schemas", "flowsheet_schema.json")

    with open(schema_path, "r") as f:
        schema = json.load(f)

    with open(config_path, "r") as f:
        config = json.load(f)

    try:
        validate(instance=config, schema=schema)
        logger.info(f"✅ Flowsheet '{config_path}' is valid.")
        # Additional connectivity checks
        check_stream_connectivity(config)
        return config
    except ValidationError as e:
        logger.error(f"❌ Validation failed for '{config_path}':")
        logger.error(f"   → {e.message}")
        if e.path:
            logger.error(f"   Path: {' → '.join(map(str, e.path))}")
        raise SystemExit(1)


def check_stream_connectivity(config):
    """
    Validates the stream connectivity in a flowsheet configuration.
    This function performs several checks on the provided configuration dictionary to ensure
    that streams are properly connected between units, feeds, and products. It verifies that
    all inlet streams have valid sources, detects unused outlets, identifies unreachable units,
    and enforces that all non-pipe units are linked through pipe units without direct connections.
    Args:
        config (dict): A dictionary containing the flowsheet configuration with keys:
            - "streams": A dict of feed/product streams.
            - "units": A dict of unit configurations, each with "type", "in", and "out" keys.
    Raises:
        ValueError: If any validation fails, such as missing fields, invalid sources,
            unreachable units, or direct connections between non-pipe units.
    Returns:
        bool: True if all validations pass.
    Notes:
        - Prints warnings for unused outlet streams.
        - Assumes "Pipes" type units are used to connect non-pipe units.
    """

    defined_streams = set(config["streams"].keys())
    produced_streams = set()
    consumed_streams = set()

    # Collect all unit input/output streams
    for name, unit in config["units"].items():
        in_stream = unit.get("in")
        out_stream = unit.get("out")

        if not in_stream or not out_stream:
            raise ValueError(f"❌ Unit '{name}' missing 'in' or 'out' field.")

        consumed_streams.add(in_stream)
        produced_streams.add(out_stream)

    _check_inlet_sources(config, defined_streams, produced_streams)
    _check_unused_outlets(config, consumed_streams, defined_streams)
    _check_unreachable_units(config, defined_streams)
    _check_pipe_linkage(config)

    return True


def _check_inlet_sources(config, defined_streams, produced_streams):
    """Check that every inlet stream has a valid source."""
    valid_sources = defined_streams.union(produced_streams)
    for name, unit in config["units"].items():
        inlet = unit["in"]
        if inlet not in valid_sources:
            raise ValueError(
                f"❌ Unit '{name}' inlet '{inlet}' not defined as a stream or previous output."
            )


def _check_unused_outlets(config, consumed_streams, defined_streams):
    """Check for unused outlets and print warnings."""
    all_known_consumers = consumed_streams.union(defined_streams)
    for name, unit in config["units"].items():
        outlet = unit["out"]
        if outlet not in all_known_consumers:
            print(
                f"⚠️  Warning: Outlet stream '{outlet}' from unit '{name}' is not consumed by any downstream unit or defined as product."
            )


def _check_unreachable_units(config, defined_streams):
    """Detect and raise error for unreachable units."""
    reachable = set(defined_streams)
    for _ in range(len(config["units"])):
        for u, ucfg in config["units"].items():
            if ucfg["in"] in reachable:
                reachable.add(ucfg["out"])

    unreachable_units = [
        u for u, cfg in config["units"].items() if cfg["out"] not in reachable
    ]

    if unreachable_units:
        raise ValueError(
            f"❌ Unreachable units detected: {', '.join(unreachable_units)}"
        )


def _check_pipe_linkage(config):
    """Check that units are linked through pipes and no direct connections between non-pipe units."""
    pipe_outputs = set()
    pipe_inputs = set()
    non_pipe_units = []

    for name, unit in config["units"].items():
        if unit["type"] == "Pipes":
            pipe_outputs.add(unit["out"])
            pipe_inputs.add(unit["in"])
        else:
            non_pipe_units.append((name, unit))

    # Non-pipe units' inputs should come from feeds or pipe outputs
    feed_streams = set(config["streams"].keys())
    valid_inputs = feed_streams.union(pipe_outputs)

    for name, unit in non_pipe_units:
        if unit["in"] not in valid_inputs:
            raise ValueError(
                f"❌ Unit '{name}' input '{unit['in']}' must come from a feed stream or a pipe output."
            )

    # Check for direct connections between non-pipe units
    non_pipe_inputs = {unit["in"] for _, unit in non_pipe_units}

    direct_connections = []
    for name, unit in non_pipe_units:
        if unit["out"] in non_pipe_inputs and unit["out"] not in feed_streams:
            # Find which unit consumes this output
            for consumer_name, consumer_unit in non_pipe_units:
                if consumer_name != name and consumer_unit["in"] == unit["out"]:
                    direct_connections.append((name, consumer_name, unit["out"]))

    if direct_connections:
        conn_str = ", ".join(
            f"'{a}' → '{b}' via '{s}'" for a, b, s in direct_connections
        )
        raise ValueError(
            f"❌ Direct connections between non-pipe units detected: {conn_str}. All units must be linked through pipes."
        )
