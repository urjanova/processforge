import graphviz
import os
from typing import Dict, Any
import sys
import json
from loguru import logger
def draw_flowsheet(
    schema: Dict[str, Any],
    output_directory: str = ".",
    output_filename: str = "flowsheet",
    file_format: str = "png",
    engine: str = "dot",
) -> str:
    """
    Generates a diagram from a flowsheet schema and saves it as an image.

    The schema is expected to be a dictionary with 'units' and 'streams'.
    - 'units': An object where keys are unit IDs and values are dictionaries
               with at least 'type', 'in', and 'out'.
    - 'streams': An object where keys are stream IDs and values are dictionaries
                 with properties like 'P', 'T', 'flowrate'.

    Args:
        schema (Dict[str, Any]): The flowsheet schema dictionary.
        output_directory (str): The directory to save the image in.
        output_filename (str): The name of the output file (without extension).
        file_format (str): The output file format (e.g., 'png', 'svg', 'pdf').
        engine (str): The layout engine to use (e.g., 'dot', 'neato', 'fdp').
                      'dot' is standard for directed graphs.

    Returns:
        str: The full path to the generated image file.

    Raises:
        KeyError: If the schema is missing 'units' or 'streams'.
        ImportError: If the 'graphviz' library is not installed.
    """


    dot = graphviz.Digraph(comment="Process Flowsheet")
    dot.attr(rankdir="LR", splines="ortho")  # Left-to-Right layout
    dot.attr("node", shape="box", style="rounded")

    units = schema.get("units", {})
    streams = schema.get("streams", {})

    # Add stream nodes
    for stream_id, stream in streams.items():
        label = f"{stream_id}\nP: {stream.get('P', 'N/A')}\nT: {stream.get('T', 'N/A')}\nFlowrate: {stream.get('flowrate', 'N/A')}"
        dot.node(stream_id, label=label, shape="ellipse", style="filled", fillcolor="lightblue")

    # Add nodes (units)
    for unit_id, unit in units.items():
        logger.info(f"Adding unit: {unit_id} with details: {unit}")
        unit_type = unit.get("type", "")
        label = f"{unit_id} ({unit_type})"
        dot.node(unit_id, label=label)

    # Build edges by matching unit 'out' to another unit 'in'
    unit_out_to_in = {}
    for unit_id, unit in units.items():
        out_stream = unit.get("out")
        if out_stream:
            unit_out_to_in[out_stream] = unit_id

    for unit_id, unit in units.items():
        in_stream = unit.get("in")
        if in_stream and in_stream in unit_out_to_in:
            source_unit = unit_out_to_in[in_stream]
            dot.edge(source_unit, unit_id)

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Render the graph
    output_path_without_ext = os.path.join(output_directory, output_filename)
    rendered_path = dot.render(
        output_path_without_ext, format=file_format, view=False, cleanup=True
    )

    return rendered_path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python flowsheet_diagram.py <flowsheet_schema.json>")
        sys.exit(1)

    with open(sys.argv[1], 'r') as f:
        flowsheet_schema = json.load(f)

    try:
        # Generate a PNG image
        png_file_path = draw_flowsheet(
            flowsheet_schema,
            output_directory="utils/diagrams",
            output_filename="process_flow_diagram",
            file_format="png",
        )
        print(f"Successfully generated PNG diagram: {png_file_path}")

        # Generate an SVG image
        svg_file_path = draw_flowsheet(
            flowsheet_schema,
            output_directory="utils/diagrams",
            output_filename="process_flow_diagram",
            file_format="svg",
        )
        print(f"Successfully generated SVG diagram: {svg_file_path}")

    except ImportError:
        print("Graphviz library not found.")
        print("Please install it using: pip install graphviz")
        print("You also need to install the Graphviz system package.")
        print("  - On Ubuntu/Debian: sudo apt-get install graphviz")
        print("  - On macOS (with Homebrew): brew install graphviz")
        print("  - On Windows (with Chocolatey): choco install graphviz")
    except Exception as e:
        print(f"An error occurred: {e}")