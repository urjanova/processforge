import json
import os
from jsonschema import validate, ValidationError

def validate_flowsheet(config_path):
    """Validate a flowsheet JSON file against the schema."""
    schema_path = os.path.join(os.path.dirname(__file__), "..", "schemas", "flowsheet_schema..json")

    with open(schema_path, "r") as f:
        schema = json.load(f)

    with open(config_path, "r") as f:
        config = json.load(f)

    try:
        validate(instance=config, schema=schema)
        print(f"✅ Flowsheet '{config_path}' is valid.")
        return config
    except ValidationError as e:
        print(f"❌ Validation failed for '{config_path}':")
        print(f"   → {e.message}")
        if e.path:
            print(f"   Path: {' → '.join(map(str, e.path))}")
        raise SystemExit(1)

def check_stream_connectivity(config):
    """Ensure all unit inputs/outputs are connected correctly."""
    defined_streams = set(config["streams"].keys())  # feeds
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

    # 1️⃣ Check that every inlet stream has a source (feed or previous unit output)
    valid_sources = defined_streams.union(produced_streams)
    for name, unit in config["units"].items():
        inlet = unit["in"]
        if inlet not in valid_sources:
            raise ValueError(f"❌ Unit '{name}' inlet '{inlet}' not defined as a stream or previous output.")

    # 2️⃣ Check for unused outlets
    all_known_consumers = consumed_streams.union(defined_streams)
    for name, unit in config["units"].items():
        outlet = unit["out"]
        if outlet not in all_known_consumers:
            print(f"⚠️  Warning: Outlet stream '{outlet}' from unit '{name}' is not consumed by any downstream unit or defined as product.")

    # 3️⃣ Optionally, detect orphaned units (never reached)
    reachable = set(defined_streams)
    for _ in range(len(config["units"])):
        for u, ucfg in config["units"].items():
            if ucfg["in"] in reachable:
                reachable.add(ucfg["out"])

    all_outputs = set(u["out"] for u in config["units"].values())
    unreachable_units = [u for u, cfg in config["units"].items() if cfg["out"] not in reachable]

    if unreachable_units:
        raise ValueError(f"❌ Unreachable units detected: {', '.join(unreachable_units)}")

    return True