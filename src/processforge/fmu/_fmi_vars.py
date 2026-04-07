"""FMI variable spec generation for Processforge FMU export."""
from __future__ import annotations

import re


def _sanitize_name(s: str) -> str:
    """Replace characters invalid in Python identifiers with underscores."""
    return re.sub(r"[^a-zA-Z0-9_]", "_", s)


def get_fmi_variable_specs(
    feed_streams: list[str],
    output_streams: list[str],
    components: list[str],
    unit_params: dict[str, dict],
    config: dict,
    mode: str,
) -> list[dict]:
    """Return ordered list of FMI variable spec dicts for a flowsheet.

    Each dict contains:
        attr_name     – Python attribute name on the slave instance
        initial_value – float initial value
        causality     – "input" | "output" | "parameter"
        variability   – "continuous" | "fixed"
        description   – human-readable string
    """
    specs: list[dict] = []

    # --- Inputs: feed stream boundary conditions ---
    for stream_name in feed_streams:
        stream_cfg = config["streams"][stream_name]
        safe_s = _sanitize_name(stream_name)

        specs.append({
            "attr_name": f"feed_{safe_s}_T",
            "initial_value": float(stream_cfg.get("T", 298.15)),
            "causality": "input",
            "variability": "continuous",
            "description": f"Temperature of feed stream '{stream_name}' [K]",
        })
        specs.append({
            "attr_name": f"feed_{safe_s}_P",
            "initial_value": float(stream_cfg.get("P", 101325.0)),
            "causality": "input",
            "variability": "continuous",
            "description": f"Pressure of feed stream '{stream_name}' [Pa]",
        })
        specs.append({
            "attr_name": f"feed_{safe_s}_flowrate",
            "initial_value": float(stream_cfg.get("flowrate", 1.0)),
            "causality": "input",
            "variability": "continuous",
            "description": f"Molar flowrate of feed stream '{stream_name}' [mol/s]",
        })
        z_cfg = stream_cfg.get("z", {})
        for comp in components:
            safe_c = _sanitize_name(comp)
            specs.append({
                "attr_name": f"feed_{safe_s}_z_{safe_c}",
                "initial_value": float(z_cfg.get(comp, 0.0)),
                "causality": "input",
                "variability": "continuous",
                "description": f"Mole fraction of {comp} in feed stream '{stream_name}'",
            })

    # --- Outputs: calculated stream properties ---
    for stream_name in output_streams:
        safe_s = _sanitize_name(stream_name)

        specs.append({
            "attr_name": f"out_{safe_s}_T",
            "initial_value": 0.0,
            "causality": "output",
            "variability": "continuous",
            "description": f"Temperature of stream '{stream_name}' [K]",
        })
        specs.append({
            "attr_name": f"out_{safe_s}_P",
            "initial_value": 0.0,
            "causality": "output",
            "variability": "continuous",
            "description": f"Pressure of stream '{stream_name}' [Pa]",
        })
        specs.append({
            "attr_name": f"out_{safe_s}_flowrate",
            "initial_value": 0.0,
            "causality": "output",
            "variability": "continuous",
            "description": f"Molar flowrate of stream '{stream_name}' [mol/s]",
        })
        for comp in components:
            safe_c = _sanitize_name(comp)
            specs.append({
                "attr_name": f"out_{safe_s}_z_{safe_c}",
                "initial_value": 0.0,
                "causality": "output",
                "variability": "continuous",
                "description": f"Mole fraction of {comp} in stream '{stream_name}'",
            })

    # --- Parameters: unit design values (fixed for FMU lifetime) ---
    for unit_name, params in unit_params.items():
        safe_u = _sanitize_name(unit_name)
        for key, value in params.items():
            if not isinstance(value, (int, float)):
                continue
            safe_k = _sanitize_name(key)
            specs.append({
                "attr_name": f"param_{safe_u}_{safe_k}",
                "initial_value": float(value),
                "causality": "parameter",
                "variability": "fixed",
                "description": f"Parameter '{key}' of unit '{unit_name}'",
            })

    # --- Tank state outputs (dynamic mode only) ---
    if mode == "dynamic":
        for unit_name, unit_cfg in config["units"].items():
            if unit_cfg.get("type") != "Tank":
                continue
            safe_u = _sanitize_name(unit_name)
            initial_T = float(unit_cfg.get("initial_T", 298.15))
            specs.append({
                "attr_name": f"state_{safe_u}_T",
                "initial_value": initial_T,
                "causality": "output",
                "variability": "continuous",
                "description": f"Tank temperature in unit '{unit_name}' [K]",
            })
            initial_n = unit_cfg.get("initial_n", {})
            for comp in components:
                safe_c = _sanitize_name(comp)
                specs.append({
                    "attr_name": f"state_{safe_u}_n_{safe_c}",
                    "initial_value": float(initial_n.get(comp, 0.0)),
                    "causality": "output",
                    "variability": "continuous",
                    "description": f"Molar holdup of {comp} in tank '{unit_name}' [mol]",
                })

    return specs
