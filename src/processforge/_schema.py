"""Internal helper to load the bundled JSON schema."""

import json
from importlib.resources import files


def load_flowsheet_schema() -> dict:
    """Load and return the flowsheet JSON schema bundled with the package."""
    schema_text = (
        files("processforge.schemas")
        .joinpath("flowsheet_schema.json")
        .read_text(encoding="utf-8")
    )
    return json.loads(schema_text)
