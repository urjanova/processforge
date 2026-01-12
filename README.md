# Process Forge

A Python-based process simulation framework for chemical process engineering applications.

## Table of Contents

- [Features](#features)
- [Available Unit Operations](#available-unit-operations)
- [Installation](#installation)
- [Usage](#usage)
- [Quick Start Examples](#quick-start-examples)
- [Project Structure](#project-structure)
- [Dependencies](#dependencies)
- [License](#license)

## Features

### Core Capabilities
- Steady-state and dynamic process simulations
- Thermodynamic property calculations using CoolProp
- Multiple unit operations for hydraulic and thermal systems
- Flowsheet modeling and solving with closed-loop (recycle) support

### Recycle Loop Support
- Automatic tear stream detection using graph analysis
- Wegstein convergence acceleration for faster solving
- Validation requirements ensure stability (Tank units required in loops)
- Configurable convergence tolerance and iteration limits

### Results & Visualization
- Export formats: CSV and JSON for easy data analysis
- Automatic plotting: temperature profiles and composition charts
- Graphviz flowsheet diagrams (PNG and SVG)
- Timeseries visualization for dynamic simulations

### Validation & Quality
- JSON schema validation for flowsheet configurations
- Connectivity checks (inlet sources, unused outlets, unreachable units)
- Comprehensive logging for debugging

## Available Unit Operations

| Unit Type | Mode | Description | Key Parameters |
|-----------|------|-------------|----------------|
| **Pump** | Steady-state | Adds pressure rise with efficiency losses | `deltaP`, `efficiency` |
| **Valve** | Steady-state | Isenthalpic pressure reduction | `pressure_ratio` |
| **Strainer** | Steady-state | Fixed pressure drop element | `deltaP` |
| **Pipes** | Steady-state & Dynamic | Laminar flow with friction losses | `delta_p`, `diameter` |
| **Tank** | Steady-state & Dynamic | Well-mixed molar tank | `inlet`, `outlet_flow`, `initial_n`, `initial_T`, `P`, `duty` |
| **Flash** | Steady-state | Isothermal flash separator | `P` |
| **Heater** | Steady-state | Temperature control unit | `duty`, `flowrate` |

## Installation

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd processforge
   ```

2. Install the package using uv:
   ```bash
   uv sync
   ```

   This will create a virtual environment (`.venv`) and install all dependencies.

   **Note:** uv is the recommended package manager for this project. If you prefer pip, you can use `pip install -r requirements.txt` for basic dependencies, but uv provides better dependency resolution and virtual environment management.

3. Activate the virtual environment:
   ```bash
   source .venv/bin/activate  # On Unix/macOS
   # or
   .venv\Scripts\activate     # On Windows
   ```

## Usage

Run simulations using the `processforge` command:

```bash
processforge flowsheets/closed-loop-chain.json
```

This will generate multiple output files in the `outputs/` directory:
- `*_results.json` - Simulation results in JSON format
- `*_timeseries.json` - Time-series data for dynamic simulations
- `*_timeseries.csv` - Tabular results with component compositions
- `temps_*.png` - Temperature profile plots
- `comps_*.png` - Composition charts for each stream

Example flowsheet files are available in the `flowsheets/` directory, including:
- `closed-loop-chain.json` - Hydraulic system with recycle loop
- `archive/example_flash.json` - Simple flash separator
- `archive/example_dynamic_tank.json` - Dynamic tank simulation

## Quick Start Examples

### Run a simulation
```bash
processforge flowsheets/closed-loop-chain.json
```

### Validate a flowsheet
```bash
python utils/validate_flowsheet.py flowsheets/closed-loop-chain.json
```

### Generate a flowsheet diagram
```bash
python utils/flowsheet_diagram.py flowsheets/closed-loop-chain.json
```

Diagrams are saved to `utils/diagrams/` as PNG and SVG files.

## Project Structure

```
processforge/
├── src/                          # Core source code
│   ├── flowsheet.py             # Flowsheet modeling with closed-loop handling
│   ├── thermo.py                # Thermodynamic calculations via CoolProp
│   ├── result.py                # Results export (CSV, JSON, plotting)
│   ├── simulate.py              # Main simulation CLI entry point
│   ├── solver.py                # Solver interface
│   └── units/                   # Unit operations
│       ├── pump.py              # Pump with efficiency
│       ├── valve.py             # Pressure-reducing valve
│       ├── strainer.py          # Pressure drop element
│       ├── pipes.py             # Pipe with friction losses
│       ├── tank.py              # Well-mixed tank (steady & dynamic)
│       ├── flash.py             # Isothermal flash separator
│       └── heater.py            # Temperature control heater
├── flowsheets/                   # Example flowsheet configurations
│   ├── closed-loop-chain.json   # Main example with recycle
│   └── archive/                 # Additional examples
├── utils/                        # Utilities
│   ├── validate_flowsheet.py    # JSON schema validation
│   ├── flowsheet_diagram.py     # Graphviz visualization
│   └── diagrams/                # Generated flowsheet diagrams
├── schemas/                      # JSON schema definitions
│   └── flowsheet_schema.json    # Flowsheet validation schema
├── outputs/                      # Simulation results
│   ├── *.json                   # Results data
│   ├── *.csv                    # Tabular data
│   └── *.png                    # Plots and charts
└── pyproject.toml               # Project configuration
```

## Dependencies

Core dependencies with their purposes:

- **numpy** - Numerical computing
- **scipy** - Scientific computing and ODE solvers
- **coolprop** - Thermodynamic property calculations
- **networkx** - Graph analysis for recycle loop detection
- **h5py** - Data storage
- **casadi** - Dynamic simulation and optimization
- **matplotlib** - Plotting and visualization
- **loguru** - Logging
- **jsonschema** - Configuration validation
- **graphviz** - Flowsheet diagram generation

## License

This project is proprietary software. See the [LICENSE](LICENSE) file for details.

For licensing inquiries, please [contact the development](https://forms.gle/wUweVnoSqA9VeD7m9) team.