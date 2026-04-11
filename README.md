# Processforge


![processforge-logo](images/processforge-logo.jpg)

A Python-based process simulation framework for coupling different simulation engines.

## Table of Contents

- [Features](#features)
- [Available Unit Operations](#available-unit-operations)
- [Installation](#installation)
- [Usage](#usage)
- [Flowsheet Configuration](#flowsheet-configuration)
- [Quick Start Examples](#quick-start-examples)
- [Plan / Apply Workflow](#plan--apply-workflow-detail)
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

### Plan / Apply Workflow
- **`pf init`**: Initialises the `.processforge/` project directory and `outputs/` folder. Run once per project.
- **`pf plan`**: Validates the flowsheet (schema, DOF, Pint unit consistency), performs a structural diff against the last saved state (`+` added, `~` modified, `-` removed units), and generates a Mermaid diagram — all without running the solver.
- **`pf apply`**: Solves the flowsheet using the last converged state as a warm start. Falls back automatically to a step-wise homotopy/continuation solver if the direct Newton solve fails. Topology changes (added/removed units) trigger a cold start with a warning.
- **Snapshot Versioning**: Every successful `apply` creates a new numbered snapshot in `.pfstate/snapshots/`. Previous snapshots are never deleted, enabling rollback to any prior converged design.
- **Convergence Guardrails**: If both the direct solve and homotopy fail, the engine auto-reverts `latest` to the last good snapshot and writes a divergence debug report (`*_divergence.json`) with the final residual norm, drifted parameters, and solver statistics.
- **Dynamic t=0 from State**: `pf run` (dynamic mode) automatically loads the latest `.pfstate` converged values as the initial conditions for time-integration, replacing arbitrary feed defaults with a physically meaningful starting point.

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

Processforge provides a CLI with the following subcommands. Both `processforge` and the shorter alias `pf` are interchangeable.

#### Recommended workflow: init → plan → apply

```bash
# 1. Initialise project (run once per project)
pf init

# 2. Preview changes: structural diff, DOF analysis, unit consistency, Mermaid diagram
pf plan flowsheets/my-flowsheet.json

# 3. Apply changes: warm-start from last state, homotopy fallback, snapshot versioned
pf apply flowsheets/my-flowsheet.json
```

#### All commands

```bash
# Initialise .processforge/ directory and outputs/ folder
pf init [--path PATH]

# Preview changes against saved state (no solver run)
pf plan flowsheets/my-flowsheet.json [--no-diagram] [--output-dir diagrams/]

# Apply flowsheet using state-based warm start and homotopy fallback
pf apply flowsheets/my-flowsheet.json

# Run a simulation directly (steady-state or dynamic)
pf run flowsheets/my-flowsheet.json [--export-images]

# Generate a standalone flowsheet diagram
pf diagram flowsheets/my-flowsheet.json [--format svg] [--output-dir diagrams/]

# Export flowsheet as Modelica .mo and optionally compile to FMU via OMPython
pf export-modelica flowsheets/my-flowsheet.json [--output-dir modelica/] [--no-compile]

# Export flowsheet as FMI 2.0 co-simulation FMU
pf export-fmu flowsheets/my-flowsheet.json [--output-dir outputs/] [--backend scipy]
```

Running `pf apply` or `pf run` generates output files in the `outputs/` directory:
- `*_results.zarr` — Simulation results stored as a Zarr directory
- `*_validation.xlsx` — Validation report derived directly from the Zarr store
- `*.pfstate/` — Versioned state store: `snapshots/` directory with one Zarr group per successful apply, plus a `latest` pointer for rollback
- `*_divergence.json` — Written when both the direct solve and homotopy fail; contains drifted parameters, final residual norm, and solver statistics

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
    "Water":   { "id": 1 },
    "Toluene": { "id": 2 },
    "Steel":   { "id": 3 }
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

Every unit also requires a `material` integer field pointing to an `id` in the `materials` section (this identifies the structural material the unit is made of, separate from the fluid composition).

### Recycle streams

Recycle streams require no special configuration. Any stream produced as the `out` of one unit can be used as the `in` of any other unit — including upstream units. The EO solver resolves the full coupled system simultaneously.

```json
"units": {
  "tank_1": { "type": "Tank", "in": ["feed", "recycle"], "out": "after_tank", ... },
  "pipe_1": { "in": "after_tank", "out": "recycle", ... }
}
```

## Quick Start Examples

### New project setup
```bash
pf init
pf plan flowsheets/hydraulic-chain.json   # validate + preview diff
pf apply flowsheets/hydraulic-chain.json  # solve + save snapshot
```

### Iterating on a design
```bash
# Edit flowsheet.json, then:
pf plan flowsheets/hydraulic-chain.json   # see what changed (+/~/-)
pf apply flowsheets/hydraulic-chain.json  # warm-start from last converged state
```

### Run a dynamic simulation (uses pfstate as t=0 if available)
```bash
pf run flowsheets/closed-loop-chain.json
```

### Generate a flowsheet diagram
```bash
pf diagram flowsheets/closed-loop-chain.json
```

## Plan / Apply Workflow Detail

The plan/apply workflow mirrors Terraform's preview-before-execute model:

| Step | Command | What happens |
|------|---------|--------------|
| 1 | `pf init` | Creates `.processforge/config.json` and `outputs/`. Run once. |
| 2 | `pf plan` | Schema + DOF validation, Pint unit checks, structural diff vs. saved state, Mermaid diagram. No solver runs. |
| 3 | `pf apply` | Loads last snapshot as warm-start x₀, runs Newton solve. Falls back to 10-step homotopy continuation if needed. Saves a new snapshot on success. |

### Structural diff output (`pf plan`)

```
+ compressor_2   [Pump]    (added)
~ pump_1         [Pump]    deltaP: 1e5 → 2e5
- old_valve      [Valve]   (removed)
```

### Snapshot versioning (`.pfstate/`)

```
outputs/my-flowsheet.pfstate/
  snapshots/
    0001_2026-04-06T12:00:00Z/   ← first apply
    0002_2026-04-06T13:30:00Z/   ← second apply
  latest                          ← plain text: "0002_2026-04-06T13:30:00Z"
```

Each snapshot stores the converged variable vector (`x`), variable names, and the flowsheet config at that point in time.

### Divergence guardrails

If both the direct Newton solve and the homotopy fallback fail to converge:
1. `latest` is automatically reverted to the previous good snapshot.
2. A `*_divergence.json` report is written with the drifted parameters, final `||F||`, homotopy step history, and the last `x` vector for debugging.


## Running on the cloud 
To run a simulation on any cloud provider (AWS EC2, Google Cloud, etc.) with Docker installed:

Pull the image:
```
docker pull ghcr.io/urjanova/processforge:latest
```
Run a simulation and map an output folder:

```Bash
  docker run -p 8080:8080 \
    -e AWS_ACCESS_KEY_ID=... \
    -e AWS_SECRET_ACCESS_KEY=... \
    -e AWS_DEFAULT_REGION=us-east-1 \
    ghcr.io/urjanova/processforge:latest \
    pf-serve
```
(The -v flag ensures the .h5 or .xdmf files generated by OpenMC/FESTIM are saved to your actual hard drive, not lost inside the container.)
## Logo credit
Google Gemini / Nano Banana

## License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

For licensing inquiries, please [contact the development](https://forms.gle/wUweVnoSqA9VeD7m9) team.
