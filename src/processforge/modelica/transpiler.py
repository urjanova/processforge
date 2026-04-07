"""Transpile a Processforge flowsheet JSON to a self-contained Modelica .mo file."""
from __future__ import annotations

import os
import re

from loguru import logger

from ..fmu.builder import _analyze_config, _get_slave_class_name
from ..utils.validate_flowsheet import validate_flowsheet
from .mo_writer import build_model_source


def transpile(
    config_path: str,
    output_dir: str = "outputs",
    output_path: str | None = None,
) -> str:
    """Translate a Processforge flowsheet JSON to a Modelica ``.mo`` file.

    Args:
        config_path:  Path to the flowsheet JSON.
        output_dir:   Directory for the generated ``.mo`` (used when
                      *output_path* is ``None``).
        output_path:  Explicit destination path for the ``.mo`` file.
                      When given, *output_dir* is ignored.

    Returns:
        Absolute path to the written ``.mo`` file.

    Raises:
        SystemExit:  If the flowsheet fails schema or connectivity validation.
        IOError:     If the output file cannot be written.
    """
    config = validate_flowsheet(config_path)

    model_name = _derive_model_name(config, config_path)
    interface = _analyze_config(config)

    source = build_model_source(model_name, config, interface)

    if output_path is None:
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{model_name}.mo")

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(source)

    abs_path = os.path.abspath(output_path)
    logger.info(f"Modelica model written to: {abs_path}")
    return abs_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _derive_model_name(config: dict, config_path: str) -> str:
    """Return a valid Modelica identifier for the top-level model.

    Uses the same logic as ``_get_slave_class_name`` from the FMU builder.
    Modelica identifiers must start with a letter and contain only
    letters, digits, and underscores.
    """
    name = config.get("metadata", {}).get("name", "")
    if not name:
        name = os.path.splitext(os.path.basename(config_path))[0]

    name = re.sub(r"[^a-zA-Z0-9]", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    if name and name[0].isdigit():
        name = "PF_" + name
    if not name:
        name = "ProcessforgeModel"

    # PascalCase
    name = "".join(part.capitalize() for part in name.split("_"))
    return name or "ProcessforgeModel"
