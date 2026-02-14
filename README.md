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

### From PyPI

```bash
pip install processforge
```

### From source (development)

1. Clone the repository:
   ```bash
   git clone https://github.com/urjanova/processforge.git
   cd processforge
   ```

2. Install the package:
   ```bash
   pip install -e ".[dev]"
   ```

   Or using uv:
   ```bash
   uv sync
   ```

## Usage

### Command Line Interface

ProcessForge provides a CLI with three subcommands:

```bash
# Run a simulation
processforge run flowsheets/closed-loop-chain.json

# Validate a flowsheet configuration
processforge validate flowsheets/closed-loop-chain.json

# Generate a flowsheet diagram
processforge diagram flowsheets/closed-loop-chain.json
processforge diagram flowsheets/closed-loop-chain.json --format svg --output-dir diagrams/
```

Running a simulation generates output files in the `outputs/` directory:
- `*_results.json` - Simulation results in JSON format
- `*_timeseries.json` - Time-series data for dynamic simulations
- `*_timeseries.csv` - Tabular results with component compositions
- `*_validation.xlsx` - Validation report

### As a Python Module

```python
from processforge import Flowsheet, Pump, Tank, Pipes, validate_flowsheet

# Load and validate a flowsheet
config = validate_flowsheet("flowsheets/closed-loop-chain.json")

# Run simulation
fs = Flowsheet(config)
results = fs.run()

# Use individual unit operations
pump = Pump("my_pump", deltaP=1e5, efficiency=0.8)
outlet = pump.run({"T": 300, "P": 101325, "flowrate": 10, "z": {"Water": 1.0}})
```

## Quick Start Examples

### Run a simulation
```bash
processforge run flowsheets/closed-loop-chain.json
```

### Validate a flowsheet
```bash
processforge validate flowsheets/closed-loop-chain.json
```

### Generate a flowsheet diagram
```bash
processforge diagram flowsheets/closed-loop-chain.json
```

## Project Structure

```
processforge/
├── src/processforge/              # Core package
│   ├── __init__.py               # Public API
│   ├── flowsheet.py              # Flowsheet modeling with closed-loop handling
│   ├── thermo.py                 # Thermodynamic calculations via CoolProp
│   ├── result.py                 # Results export (CSV, JSON, Excel, plotting)
│   ├── simulate.py               # CLI entry point with subcommands
│   ├── solver.py                 # Solver interface
│   ├── validate.py               # Simple schema validation
│   ├── _schema.py                # Schema loader (importlib.resources)
│   ├── units/                    # Unit operations
│   │   ├── pump.py               # Pump with efficiency
│   │   ├── valve.py              # Pressure-reducing valve
│   │   ├── strainer.py           # Pressure drop element
│   │   ├── pipes.py              # Pipe with friction losses
│   │   ├── tank.py               # Well-mixed tank (steady & dynamic)
│   │   ├── flash.py              # Isothermal flash separator
│   │   └── heater.py             # Temperature control heater
│   ├── utils/                    # Utilities
│   │   ├── validate_flowsheet.py # Schema + connectivity validation
│   │   └── flowsheet_diagram.py  # Graphviz visualization
│   └── schemas/                  # Bundled JSON schemas
│       └── flowsheet_schema.json
├── flowsheets/                    # Example flowsheet configurations
│   ├── closed-loop-chain.json    # Main example with recycle
│   └── archive/                  # Additional examples
├── pyproject.toml                # Project configuration
└── MANIFEST.in                   # Source distribution manifest
```

## Dependencies

Core dependencies:

- **numpy** - Numerical computing
- **scipy** - Scientific computing and ODE solvers
- **coolprop** - Thermodynamic property calculations
- **matplotlib** - Plotting and visualization
- **loguru** - Logging
- **jsonschema** - Configuration validation
- **graphviz** - Flowsheet diagram generation
- **pandas** - Data manipulation
- **openpyxl** - Excel report generation

## License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

For licensing inquiries, please [contact the development](https://forms.gle/wUweVnoSqA9VeD7m9) team.
