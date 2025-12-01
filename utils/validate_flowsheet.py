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


def _get_unit_inlets(unit_config):
    """Get inlet stream(s) for a unit, handling both single and multiple inlets."""
    inlet = unit_config.get("in")
    if isinstance(inlet, list):
        return inlet
    return [inlet] if inlet else []


def check_stream_connectivity(config):
    """
    Validates the stream connectivity in a flowsheet configuration.
    This function performs several checks on the provided configuration dictionary to ensure
    that streams are properly connected between units, feeds, and products. It verifies that
    all inlet streams have valid sources, detects unused outlets, identifies unreachable units,
    and enforces that all non-pipe units are linked through pipe units without direct connections.
    
    For flowsheets with recycles, it also validates that each recycle loop contains at least
    one Tank unit for stability.
    
    Args:
        config (dict): A dictionary containing the flowsheet configuration with keys:
            - "streams": A dict of feed/product streams.
            - "units": A dict of unit configurations, each with "type", "in", and "out" keys.
            - "recycle_streams" (optional): A dict of recycle stream initial guesses.
    Raises:
        ValueError: If any validation fails, such as missing fields, invalid sources,
            unreachable units, or direct connections between non-pipe units.
    Returns:
        bool: True if all validations pass.
    Notes:
        - Prints warnings for unused outlet streams.
        - Assumes "Pipes" type units are used to connect non-pipe units.
        - Recycle loops must contain at least one Tank unit.
    """
    defined_streams = set(config["streams"].keys())
    recycle_streams = set(config.get("recycle_streams", {}).keys())
    produced_streams = set()
    consumed_streams = set()

    # Collect all unit input/output streams
    for name, unit in config["units"].items():
        inlets = _get_unit_inlets(unit)
        out_stream = unit.get("out")

        if not inlets or not out_stream:
            raise ValueError(f"❌ Unit '{name}' missing 'in' or 'out' field.")

        for inlet in inlets:
            consumed_streams.add(inlet)
        produced_streams.add(out_stream)

    _check_inlet_sources(config, defined_streams, produced_streams, recycle_streams)
    _check_unused_outlets(config, consumed_streams, defined_streams)
    _check_unreachable_units(config, defined_streams, recycle_streams)
    _check_pipe_linkage(config, recycle_streams)
    _check_recycle_has_tank(config, recycle_streams)

    return True


def _check_inlet_sources(config, defined_streams, produced_streams, recycle_streams):
    """Check that every inlet stream has a valid source."""
    valid_sources = defined_streams.union(produced_streams).union(recycle_streams)
    for name, unit in config["units"].items():
        inlets = _get_unit_inlets(unit)
        for inlet in inlets:
            if inlet not in valid_sources:
                raise ValueError(
                    f"❌ Unit '{name}' inlet '{inlet}' not defined as a stream, recycle stream, or previous output."
                )


def _check_unused_outlets(config, consumed_streams, defined_streams):
    """Check for unused outlets and print warnings."""
    recycle_streams = set(config.get("recycle_streams", {}).keys())
    all_known_consumers = consumed_streams.union(defined_streams).union(recycle_streams)
    for name, unit in config["units"].items():
        outlet = unit["out"]
        if outlet not in all_known_consumers:
            print(
                f"⚠️  Warning: Outlet stream '{outlet}' from unit '{name}' is not consumed by any downstream unit or defined as product."
            )


def _check_unreachable_units(config, defined_streams, recycle_streams):
    """Detect and raise error for unreachable units."""
    # Start with feed streams and recycle streams as reachable
    reachable = set(defined_streams).union(recycle_streams)
    
    # Propagate reachability
    for _ in range(len(config["units"])):
        for u, ucfg in config["units"].items():
            inlets = _get_unit_inlets(ucfg)
            if any(inlet in reachable for inlet in inlets):
                reachable.add(ucfg["out"])

    unreachable_units = [
        u for u, cfg in config["units"].items() if cfg["out"] not in reachable
    ]

    if unreachable_units:
        raise ValueError(
            f"❌ Unreachable units detected: {', '.join(unreachable_units)}"
        )


def _check_pipe_linkage(config, recycle_streams):
    """Check that units are linked through pipes and no direct connections between non-pipe units."""
    pipe_outputs = set()
    pipe_inputs = set()
    non_pipe_units = []

    for name, unit in config["units"].items():
        if unit["type"] == "Pipes":
            pipe_outputs.add(unit["out"])
            inlets = _get_unit_inlets(unit)
            for inlet in inlets:
                pipe_inputs.add(inlet)
        else:
            non_pipe_units.append((name, unit))

    # Non-pipe units' inputs should come from feeds, recycle streams, or pipe outputs
    feed_streams = set(config["streams"].keys())
    valid_inputs = feed_streams.union(pipe_outputs).union(recycle_streams)

    for name, unit in non_pipe_units:
        inlets = _get_unit_inlets(unit)
        for inlet in inlets:
            if inlet not in valid_inputs:
                raise ValueError(
                    f"❌ Unit '{name}' input '{inlet}' must come from a feed stream, recycle stream, or a pipe output."
                )

    # Check for direct connections between non-pipe units
    non_pipe_inputs = set()
    for _, unit in non_pipe_units:
        inlets = _get_unit_inlets(unit)
        non_pipe_inputs.update(inlets)

    direct_connections = []
    for name, unit in non_pipe_units:
        if unit["out"] in non_pipe_inputs and unit["out"] not in feed_streams and unit["out"] not in recycle_streams:
            # Find which unit consumes this output
            for consumer_name, consumer_unit in non_pipe_units:
                consumer_inlets = _get_unit_inlets(consumer_unit)
                if consumer_name != name and unit["out"] in consumer_inlets:
                    direct_connections.append((name, consumer_name, unit["out"]))

    if direct_connections:
        conn_str = ", ".join(
            f"'{a}' → '{b}' via '{s}'" for a, b, s in direct_connections
        )
        raise ValueError(
            f"❌ Direct connections between non-pipe units detected: {conn_str}. All units must be linked through pipes."
        )


def _check_recycle_has_tank(config, recycle_streams):
    """
    Validate that each recycle loop contains at least one Tank unit.
    Tanks provide the necessary holdup for recycle stability.
    """
    if not recycle_streams:
        return
    
    # Build a graph to detect cycles and check for Tanks
    stream_producers = {}
    stream_consumers = {}
    
    for unit_name, unit_config in config["units"].items():
        out_stream = unit_config["out"]
        stream_producers[out_stream] = unit_name
        
        inlets = _get_unit_inlets(unit_config)
        for inlet in inlets:
            if inlet not in stream_consumers:
                stream_consumers[inlet] = []
            stream_consumers[inlet].append(unit_name)
    
    # For each recycle stream, trace the cycle and check for Tank
    for recycle_stream in recycle_streams:
        # Find the unit that produces this recycle stream
        if recycle_stream not in stream_producers:
            continue
        
        producer = stream_producers[recycle_stream]
        
        # Find the unit that consumes this recycle stream
        if recycle_stream not in stream_consumers:
            continue
        
        consumers = stream_consumers[recycle_stream]
        
        # Trace from consumer back to producer to find the cycle
        for consumer in consumers:
            cycle_units = _trace_cycle(config, consumer, producer, stream_producers)
            if cycle_units:
                # Check if any unit in the cycle is a Tank
                has_tank = any(
                    config["units"][u]["type"] == "Tank"
                    for u in cycle_units
                    if u in config["units"]
                )
                
                if not has_tank:
                    cycle_str = " -> ".join(cycle_units)
                    raise ValueError(
                        f"❌ Recycle loop must contain at least one Tank unit for stability. "
                        f"Cycle without Tank: {cycle_str}"
                    )


def _trace_cycle(config, start_unit, end_unit, stream_producers):
    """
    Trace the path from start_unit to end_unit through the flowsheet.
    Returns the list of units in the cycle, or empty list if no path found.
    """
    visited = set()
    path = []
    
    def dfs(unit):
        if unit in visited:
            return False
        if unit == end_unit:
            path.append(unit)
            return True
        
        visited.add(unit)
        path.append(unit)
        
        # Get the output stream of this unit
        out_stream = config["units"].get(unit, {}).get("out")
        if not out_stream:
            path.pop()
            return False
        
        # Find units that consume this output
        for next_unit, next_config in config["units"].items():
            next_inlets = _get_unit_inlets(next_config)
            if out_stream in next_inlets:
                if dfs(next_unit):
                    return True
        
        path.pop()
        return False
    
    if dfs(start_unit):
        return path
    return []
