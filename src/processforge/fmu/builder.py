"""Build a PythonFMU co-simulation FMU from a Processforge flowsheet config."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile

from ..utils.validate_flowsheet import validate_flowsheet
from ..utils.topology import TOPOLOGY_KEYS, get_outlets
from ._fmi_vars import _sanitize_name
from .slave_template import render_slave_source


def build_fmu(
    config_path: str,
    output_dir: str = "outputs",
    backend: str = "scipy",
) -> str:
    """Build an FMI 2.0 co-simulation FMU from a Processforge flowsheet JSON.

    Args:
        config_path: Path to the flowsheet JSON file.
        output_dir:  Directory where the ``.fmu`` file will be written.
        backend:     EO solver backend for steady-state mode
                     (``"scipy"``, ``"pyomo"``, or ``"casadi"``).

    Returns:
        Absolute path to the generated ``.fmu`` file.

    Raises:
        RuntimeError: If ``pythonfmu`` is not installed or the build fails.
    """
    if not shutil.which("pythonfmu"):
        raise RuntimeError(
            "PythonFMU is not installed or not on PATH. "
            "Install it with:  pip install processforge[fmu]"
        )

    config = validate_flowsheet(config_path)
    interface = _analyze_config(config)
    slave_class_name = _get_slave_class_name(config, config_path)

    os.makedirs(output_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as staging_dir:
        # Place the config JSON in the staging directory so PythonFMU bundles
        # it into the FMU's resources/ directory.
        config_staging = os.path.join(staging_dir, "flowsheet_config.json")
        shutil.copy(config_path, config_staging)

        slave_source = render_slave_source(
            slave_class_name=slave_class_name,
            mode=interface["mode"],
            backend=backend,
            feed_streams=interface["feed_streams"],
            output_streams=interface["output_streams"],
            components=interface["components"],
            unit_params=interface["unit_params"],
            config=config,
        )

        slave_py = os.path.join(staging_dir, f"{slave_class_name}.py")
        with open(slave_py, "w") as f:
            f.write(slave_source)

        try:
            result = subprocess.run(
                [
                    "pythonfmu",
                    "build",
                    "-f", slave_py,
                    "-d", os.path.abspath(output_dir),
                    config_staging,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"pythonfmu build failed:\n{exc.stderr}"
            ) from exc

    fmu_path = os.path.abspath(
        os.path.join(output_dir, f"{slave_class_name}.fmu")
    )
    return fmu_path


# ---------------------------------------------------------------------------
# Config analysis
# ---------------------------------------------------------------------------

def _analyze_config(config: dict) -> dict:
    """Extract the FMI interface description from a validated config dict.

    Returns a dict with:
        feed_streams   – list of feed stream names (from config["streams"])
        output_streams – all unit outlet stream names (not feeds)
        components     – sorted list of component names
        unit_params    – {unit_name: {param_key: value}} for numeric params
        mode           – "steady" or "dynamic"
    """
    feed_streams: list[str] = list(config["streams"].keys())

    # Collect all outlet stream names from unit definitions
    all_outlet_streams: set[str] = set()
    for unit_cfg in config["units"].values():
        for outlet in get_outlets(unit_cfg):
            all_outlet_streams.add(outlet)

    # Output streams = outlets not already in feeds
    feed_set = set(feed_streams)
    output_streams: list[str] = [s for s in all_outlet_streams if s not in feed_set]
    # Preserve stable order (insertion order of config["units"])
    seen: set[str] = set()
    output_streams_ordered: list[str] = []
    for unit_cfg in config["units"].values():
        for outlet in get_outlets(unit_cfg):
            if outlet not in feed_set and outlet not in seen:
                output_streams_ordered.append(outlet)
                seen.add(outlet)

    # Component discovery — same logic as EOFlowsheet._collect_components()
    comp_set: set[str] = set()
    for stream in config["streams"].values():
        comp_set.update(stream.get("z", {}).keys())
    components = sorted(comp_set)

    # Unit parameters — numeric config keys, excluding topology keys
    unit_params: dict[str, dict] = {}
    for unit_name, unit_cfg in config["units"].items():
        params = {
            k: v
            for k, v in unit_cfg.items()
            if k not in TOPOLOGY_KEYS and isinstance(v, (int, float))
        }
        if params:
            unit_params[unit_name] = params

    # Simulation mode — force "dynamic" if any Tank is present
    mode = config.get("simulation", {}).get("mode", "steady")
    has_tank = any(
        cfg.get("type") == "Tank" for cfg in config["units"].values()
    )
    if has_tank:
        mode = "dynamic"

    return {
        "feed_streams": feed_streams,
        "output_streams": output_streams_ordered,
        "components": components,
        "unit_params": unit_params,
        "mode": mode,
    }

def _get_slave_class_name(config: dict, config_path: str) -> str:
    """Derive a valid Python class name for the FMU slave."""
    name = config.get("metadata", {}).get("name", "")
    if not name:
        name = os.path.splitext(os.path.basename(config_path))[0]

    # Convert to PascalCase-safe identifier
    name = re.sub(r"[^a-zA-Z0-9]", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    # Ensure starts with a letter
    if name and name[0].isdigit():
        name = "FMU_" + name
    if not name:
        name = "ProcessforgeFMU"

    # PascalCase: capitalise each word segment
    name = "".join(part.capitalize() for part in name.split("_"))
    return name or "ProcessforgeFMU"
