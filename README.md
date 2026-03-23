# Process Forge


![processforge-logo](images/processforge-logo.jpg)

A Python-based process simulation framework for chemical process engineering applications.

## Table of Contents

- [Features](#features)
- [Available Unit Operations](#available-unit-operations)
- [Installation](#installation)
- [Usage](#usage)
- [Flowsheet Configuration](#flowsheet-configuration)
- [Quick Start Examples](#quick-start-examples)
- [Project Structure](#project-structure)
- [Dependencies](#dependencies)
- [License](#license)

## Features

### Core Capabilities
- Steady-state EO (equation-oriented) and dynamic process simulations
- Thermodynamic property calculations using CoolProp
- Multiple unit operations for hydraulic and thermal systems
- Flowsheet modeling with automatic recycle stream detection
- Modelica/OMPython bridge: transpile flowsheets to native `.mo` source and compile to Model Exchange FMU via OpenModelica

### Steady-State EO Solver
- All unit equations assembled into a global `F(x) = 0` system and solved simultaneously
- Newton-Raphson with Armijo backtracking line search (SciPy backend)
- Pluggable backends: **SciPy** (built-in), **Pyomo + IPOPT** (optional), **CasADi** (optional)
- Recycle loops handled natively — no tear stream initialisation or Wegstein iteration needed
- Recycle streams are auto-detected from the flowsheet graph topology

### Dynamic Simulation
- ODE time-marching for Tank-based dynamic flowsheets
- Sequential-modular propagation with Wegstein convergence for recycle loops
- Timeseries results for transient analysis

### Results & Visualization
- Export formats: Zarr stores containing both scalar and timeseries results for downstream analysis
- Optional plotting (use `--export-images`) for temperature profiles and composition charts
- Graphviz flowsheet diagrams (PNG and SVG)
- Timeseries visualization for dynamic simulations

### Validation & Quality
- JSON schema validation for flowsheet configurations
- Connectivity checks (inlet sources, unused outlets, unreachable units)
- Comprehensive logging for debugging

## Available Unit Operations

| Unit Type | Mode | Description | Key Parameters |
|-----------|------|-------------|----------------|
| **Pump** | Steady-state (EO) | Adds pressure rise with efficiency losses | `deltaP`, `efficiency` |
| **Valve** | Steady-state (EO) | Isenthalpic pressure reduction | `pressure_ratio` |
| **Strainer** | Steady-state (EO) | Fixed pressure drop element | `deltaP` |
| **Pipes** | Steady-state (EO) & Dynamic | Laminar flow with friction losses | `delta_p`, `diameter` |
| **Tank** | Dynamic only | Well-mixed molar tank (ODE) | `outlet_flow`, `initial_n`, `initial_T`, `P`, `duty` |
| **Flash** | Steady-state (EO, SciPy backend) | Isothermal flash separator | `P` |
| **Heater** | Steady-state (EO, SciPy backend) | Temperature control unit | `duty`, `flowrate` |

> **Note:** Flash and Heater use CoolProp internally and are supported on the SciPy backend only.
> Pump, Valve, Strainer, and Pipes support all three backends (SciPy, Pyomo, CasADi).

## Installation

### From PyPI

```bash
# pip
pip install processforge

# uv
uv add processforge
```

### With EO solver backends (optional)

```bash
# Pyomo + IPOPT backend
pip install "processforge[eo]"
uv add "processforge[eo]"

# Pyomo + IPOPT + CasADi backend
pip install "processforge[eo-casadi]"
uv add "processforge[eo-casadi]"

# Modelica transpiler + OMPython bridge (requires OpenModelica installed separately)
pip install "processforge[modelica]"
uv add "processforge[modelica]"
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
# Run a simulation (add --export-images to generate PNG plots)
processforge run flowsheets/closed-loop-chain.json [--export-images]

# Validate a flowsheet configuration
processforge validate flowsheets/closed-loop-chain.json

# Generate a flowsheet diagram
processforge diagram flowsheets/closed-loop-chain.json
processforge diagram flowsheets/closed-loop-chain.json --format svg --output-dir diagrams/

# Export flowsheet as Modelica .mo and compile to Model Exchange FMU via OMPython
processforge export-modelica flowsheets/my-flowsheet.json
processforge export-modelica flowsheets/my-flowsheet.json --output-dir modelica/ --no-compile
```

Running a simulation generates output files in the `outputs/` directory:
- `*_results.zarr` - Simulation results stored as a Zarr directory
- `*_validation.xlsx` - Validation report derived directly from the Zarr store
- PNG plots for temperatures and compositions when `--export-images` is supplied

### As a Python Module

```python
from processforge import EOFlowsheet, validate_flowsheet

# Load and validate a flowsheet
config = validate_flowsheet("flowsheets/my-flowsheet.json")

# Run steady-state EO simulation (default SciPy backend)
fs = EOFlowsheet(config, backend="scipy")
results = fs.run()
# results: {stream_name: {"T": ..., "P": ..., "flowrate": ..., "z": {...}}}

# Use Pyomo + IPOPT backend (requires: pip install "processforge[eo]")
fs = EOFlowsheet(config, backend="pyomo")
results = fs.run()

# Use CasADi backend (requires: pip install "processforge[eo-casadi]")
fs = EOFlowsheet(config, backend="casadi")
results = fs.run()

```

## Flowsheet Configuration

Flowsheets are defined as JSON files. The `simulation.mode` field controls which solver is used:

| `mode` | Solver | Use case |
|--------|--------|----------|
| `"steady"` (default) | EO — global Newton-Raphson | Steady-state without Tank units |
| `"dynamic"` | SM — ODE time-marching | Flowsheets containing Tank units |

```json
{
  "metadata": { "name": "My Flowsheet", "version": "2.0" },
  "materials": {
    "Water":   { "friendly_material_id": 1 },
    "Toluene": { "friendly_material_id": 2 },
    "Steel":   { "friendly_material_id": 3 }
  },
  "material_mixes": {
    "Water_Toluene_Mix": {
      "friendly_material_mix_id": 1,
      "percent_type": "ao",
      "components": [
        { "name": "Water",   "fraction": 0.8 },
        { "name": "Toluene", "fraction": 0.2 }
      ]
    }
  },
  "streams": {
    "feed": { "T": 298.15, "P": 101325, "flowrate": 1.0, "material_mix": 1 }
  },
  "units": {
    "pump_1":  { "type": "Pump",  "in": "feed",       "out": "after_pump", "deltaP": 200000, "efficiency": 0.75, "material": 3 },
    "valve_1": { "type": "Valve", "in": "after_pump",  "out": "product",   "pressure_ratio": 0.5,                "material": 3 }
  },
  "simulation": {
    "mode": "steady",
    "backend": "scipy"
  }
}
```

The optional `backend` key selects the EO solver backend (`"scipy"`, `"pyomo"`, or `"casadi"`). Defaults to `"scipy"`.

### Materials and composition

Stream composition can be defined in two ways:

**Explicit `z` dict** — list component mole fractions directly on the stream:
```json
"feed": { "T": 298.15, "P": 101325, "flowrate": 1.0, "z": { "Water": 0.8, "Toluene": 0.2 } }
```

**`material_mix` reference** — define a reusable mix in the top-level `material_mixes` section and reference it by `friendly_material_mix_id`. The validator automatically expands the reference into a `z` dict before simulation:
```json
"material_mixes": {
  "Water_Toluene_Mix": {
    "friendly_material_mix_id": 1,
    "percent_type": "ao",
    "components": [
      { "name": "Water",   "fraction": 0.8 },
      { "name": "Toluene", "fraction": 0.2 }
    ]
  }
},
"streams": {
  "feed": { "T": 298.15, "P": 101325, "flowrate": 1.0, "material_mix": 1 }
}
```

Rules:
- `z` and `material_mix` are mutually exclusive on a single stream.
- `friendly_material_mix_id` values must be unique across all mixes.
- Each component `name` in a mix must match a key in the top-level `materials` section.
- When all component fractions are provided they must sum to 1.0.

Every unit also requires a `material` integer field pointing to a `friendly_material_id` in the `materials` section (this identifies the structural material the unit is made of, separate from the fluid composition).

### Recycle streams

Recycle streams require no special configuration. Any stream produced as the `out` of one unit can be used as the `in` of any other unit — including upstream units. The EO solver resolves the full coupled system simultaneously.

```json
"units": {
  "tank_1": { "type": "Tank", "in": ["feed", "recycle"], "out": "after_tank", ... },
  "pipe_1": { "in": "after_tank", "out": "recycle", ... }
}
```

## Quick Start Examples

### Run a steady-state simulation
```bash
processforge run flowsheets/hydraulic-chain.json
```

### Run a dynamic simulation
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
│   ├── flowsheet.py              # Sequential-modular solver (dynamic mode)
│   ├── thermo.py                 # Thermodynamic calculations via CoolProp
│   ├── result.py                 # Results export (Zarr, Excel, plotting)
│   ├── simulate.py               # CLI entry point with subcommands
│   ├── solver.py                 # ODE solver interface (dynamic)
│   ├── validate.py               # Simple schema validation
│   ├── _schema.py                # Schema loader (importlib.resources)
│   ├── eo/                       # Equation-oriented (EO) steady-state solver
│   │   ├── flowsheet.py          # EOFlowsheet — build, warm-start, solve
│   │   ├── solver.py             # EOSolver — backend selector
│   │   ├── jacobian.py           # GlobalJacobianManager — F(x), J(x)
│   │   ├── stream_var.py         # StreamVar — per-stream variable container
│   │   ├── mixin.py              # EOUnitModelMixin — unit residual interface
│   │   ├── backends/             # Pluggable solver backends
│   │   │   ├── scipy_backend.py  # Newton-Raphson + Armijo (built-in)
│   │   │   ├── pyomo_backend.py  # Pyomo ConcreteModel + IPOPT (optional)
│   │   │   └── casadi_backend.py # CasADi SX + rootfinder (optional)
│   │   └── units/                # EO residual equations per unit type
│   │       ├── pump_eo.py
│   │       ├── valve_eo.py
│   │       ├── strainer_eo.py
│   │       ├── pipes_eo.py
│   │       ├── heater_eo.py
│   │       └── flash_eo.py
│   ├── units/                    # Unit operation implementations
│   │   ├── pump.py               # Pump with efficiency
│   │   ├── valve.py              # Pressure-reducing valve
│   │   ├── strainer.py           # Pressure drop element
│   │   ├── pipes.py              # Pipe with friction losses
│   │   ├── tank.py               # Well-mixed tank (dynamic ODE)
│   │   ├── flash.py              # Isothermal flash separator
│   │   └── heater.py             # Temperature control heater
│   ├── modelica/                 # OMPython bridge: .mo transpiler + omc runner
│   │   ├── transpiler.py         # transpile() — config → Modelica source
│   │   ├── mo_writer.py          # build_model_source() — string builder
│   │   ├── unit_equations.py     # Per-unit equation generators
│   │   └── omc_runner.py         # compile_modelica() — OMPython validate + FMU
│   ├── utils/                    # Utilities
│   │   ├── validate_flowsheet.py # Schema + connectivity validation
│   │   └── flowsheet_diagram.py  # Graphviz visualization
│   └── schemas/                  # Bundled JSON schemas
│       └── flowsheet_schema.json
├── flowsheets/                    # Example flowsheet configurations
│   ├── closed-loop-chain.json    # Dynamic recycle example
│   └── archive/                  # Additional examples
├── pyproject.toml                # Project configuration
└── MANIFEST.in                   # Source distribution manifest
```

## Dependencies

Core dependencies (always installed):

- **numpy** - Numerical computing
- **scipy** - Scientific computing, ODE solvers, sparse linear algebra
- **coolprop** - Thermodynamic property calculations
- **matplotlib** - Plotting and visualization
- **loguru** - Logging
- **jsonschema** - Configuration validation
- **graphviz** - Flowsheet diagram generation
- **pandas** - Data manipulation
- **openpyxl** - Excel report generation
- **zarr** - Chunked storage for simulation outputs

Optional EO solver backends:

- **pyomo** ≥ 6.7 — Pyomo + IPOPT backend (`pip install "processforge[eo]"`)
- **casadi** ≥ 3.6 — CasADi AD-based backend (`pip install "processforge[eo-casadi]"`)
- **OMPython** ≥ 1.4 — Python API to OpenModelica compiler (`pip install "processforge[modelica]"`). Requires [OpenModelica](https://openmodelica.org) installed on the system.


## Logo credit
Google Gemini / Nano Banana

## License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

For licensing inquiries, please [contact the development](https://forms.gle/wUweVnoSqA9VeD7m9) team.
