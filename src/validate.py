import json
import jsonschema
from jsonschema import validate
from loguru import logger

def validate_flowsheet(config_path, schema_path="schemas/flowsheet_schema.json"):
    """
    Validates a flowsheet JSON against the defined schema.
    Raises jsonschema.ValidationError if invalid.
    """
    with open(config_path, "r") as f:
        config = json.load(f)
    with open(schema_path, "r") as f:
        schema = json.load(f)

    try:
        validate(instance=config, schema=schema)
        logger.info(f"Flowsheet '{config_path}' validated successfully.")
        return config
    except jsonschema.exceptions.ValidationError as err:
        logger.error(f"Validation error in '{config_path}': {err.message}")
        raise
