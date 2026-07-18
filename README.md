# Processforge

![processforge-logo](images/processforge-logo.svg)

A lightweight Python framework for process simulation, coupling hydraulic, thermal, and reactor workflows.

## Install

Install the `pf` command-line tool with [uv](https://docs.astral.sh/uv/):

```bash
uv tool install processforge
```

For CoolProp-backed units:

```bash
uv tool install "processforge[coolprop]"
```

Optional solver backends:

```bash
uv tool install "processforge[eo]"
uv tool install "processforge[eo-casadi]"
uv tool install "processforge[modelica]"
```

## Quick start

1. **Install the tool**

   ```bash
   uv tool install processforge
   ```

2. **Download an example flowsheet**

   ```bash
   curl -O https://raw.githubusercontent.com/urjanova/processforge/master/flowsheets/hydraulic-chain.json
   ```

3. **Initialize, plan, and apply**

   ```bash
   pf init
   pf plan hydraulic-chain.json
   pf apply hydraulic-chain.json
   ```

   `plan` validates the flowsheet (schema, DOF, units) without running the solver; `apply` solves it and stores a snapshot.

4. **Look at the output**

   `pf apply` writes results under `outputs/`:
   - `*_results.zarr` — simulation results store (per-variable arrays, composition flattened)
   - `*_results.zarr.schema.json` — schema file describing streams, variables, dtypes, units, shapes, and run provenance
   - `*_validation.xlsx` — validation report
   - `*.pfstate/` — versioned snapshot store with a `latest` pointer

   For flowsheets with Tank units (dynamic), use `pf run` instead of `pf apply` to solve with the SM solver.

## Python API

```python
from processforge import EOFlowsheet, validate_flowsheet
config = validate_flowsheet("flowsheets/hydraulic-chain.json")
fs = EOFlowsheet(config, backend="scipy")
results = fs.run()
```

## Usage
See the [usage guide](docs/usage.md) for CLI commands and workflows.

Individual providers (e.g. FESTIM, OpenMC) can be run via Docker images using the provider image contract.

## Flowsheets
The core of Processforge is the flowsheet JSON format, which defines materials, units, and recycle rules. See below for more information on the formats and example flowsheets.
- [docs/flowsheets.md](docs/flowsheets.md) : flowsheet JSON format, materials, units, and recycle rules
- [flowsheets/](flowsheets/) : example flowsheets shipped with the repo


## License

BSD 3-Clause License
