"""Generate a Mermaid.js flowchart diagram from a flowsheet config."""
from __future__ import annotations
import os


def generate_mermaid(config: dict, output_path: str | None = None) -> str:
    """
    Build a Mermaid flowchart LR diagram of the flowsheet topology.

    Node shapes:
      - Feed streams:        ([name<br>feed stream])  — stadium
      - Standard units:      [name<br>Type]            — rectangle
      - Flash (2 outlets):   {name<br>Flash}           — diamond

    Edge labels show the intermediate stream name.

    Args:
        config:      Validated flowsheet config dict.
        output_path: If given, write the .mmd source to this path.

    Returns:
        The Mermaid source string.
    """
    lines: list[str] = ["flowchart LR"]

    feed_streams = set(config.get("streams", {}).keys())
    units = config.get("units", {})

    # Declare feed stream nodes
    for stream_name in feed_streams:
        node_id = _node_id(stream_name)
        lines.append(f'    {node_id}(["📥 {stream_name}<br>feed stream"])')

    # Declare unit nodes
    for unit_name, unit_cfg in units.items():
        unit_type = unit_cfg.get("type", unit_name)
        node_id = _node_id(unit_name)
        if unit_type == "Flash":
            lines.append(f'    {node_id}{{"{unit_name}<br>{unit_type}"}}')
        else:
            lines.append(f'    {node_id}["{unit_name}<br>{unit_type}"]')

    lines.append("")  # blank line before edges

    # Build stream-producer map: stream_name → producer node_id
    # Feed streams are produced by their own node; unit outputs by unit node
    producer: dict[str, str] = {s: _node_id(s) for s in feed_streams}
    for unit_name, unit_cfg in units.items():
        producer[unit_cfg["out"]] = _node_id(unit_name)

    # Draw edges: for each unit, one edge per inlet
    for unit_name, unit_cfg in units.items():
        inlets = unit_cfg.get("in", [])
        if isinstance(inlets, str):
            inlets = [inlets]
        dest = _node_id(unit_name)
        for inlet_stream in inlets:
            src = producer.get(inlet_stream, _node_id(inlet_stream))
            lines.append(f'    {src} -->|"{inlet_stream}"| {dest}')

    mermaid_src = "\n".join(lines) + "\n"

    if output_path:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(mermaid_src)

    return mermaid_src


def _node_id(name: str) -> str:
    """Convert a stream/unit name to a valid Mermaid node identifier."""
    # Replace characters that break Mermaid syntax
    return name.replace("-", "_").replace(" ", "_").replace(".", "_")
