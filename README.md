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
   pf init flowsheets/hydraulic-chain.json
   pf plan flowsheets/hydraulic-chain.json
   pf apply flowsheets/hydraulic-chain.json
   ```

   `plan` validates the flowsheet (schema, DOF, units) without running the solver; `apply` solves it and stores a snapshot.

4. **Look at the output**

    `pf apply` (and `pf run`) write a single unified archive under `outputs/`:

    - `outputs/<base>.pfarchive/` — the unified store for a flowsheet's solved state and outputs:
      - `snapshots/` — Zarr store of converged state vectors (one group per successful `pf apply`), each with `x`/`x_delta` arrays and config/var-name/metadata attributes, plus a `latest` pointer (powers warm-start, drift detection, and homotopy).
      - `runs/<run_id>.json` — the run manifest: every stream and unit engine output (values, units, dtypes, shapes) and run provenance (backend, version, flowsheet hash).
      - `outputs/streams/<name>.json` — per-stream timeseries from the solve.
      - `artifacts.json` — content-addressed registry of all output artifacts (local + remote URIs).
      - `index.json` — `field_name → occurrences` index for fast cross-run lookups.
      - `latest_run` — plain-text pointer to the most recent run.
    - `outputs/<base>_divergence.json` — written only when both direct and homotopy solves fail on `pf apply`, capturing drifted params, solver stats, and top residual violators.

    For dynamic flowsheets with Tank units, use `pf run` (not `pf apply`, which is steady-state EO only) to solve with the dynamic engine.

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
