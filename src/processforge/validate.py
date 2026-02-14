import json
import jsonschema
from jsonschema import validate
from loguru import logger

from ._schema import load_flowsheet_schema


def validate_flowsheet(config_path, schema=None):
    """
    Validates a flowsheet JSON against the defined schema.
    Raises jsonschema.ValidationError if invalid.
    """
    if schema is None:
        schema = load_flowsheet_schema()

    with open(config_path, "r") as f:
        config = json.load(f)

    try:
        validate(instance=config, schema=schema)
        logger.info(f"Flowsheet '{config_path}' validated successfully.")
        return config
    except jsonschema.exceptions.ValidationError as err:
        logger.error(f"Validation error in '{config_path}': {err.message}")
        raise
