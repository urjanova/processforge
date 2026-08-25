# Usage

## Plan / Apply Workflow

1. `pf init` — creates `.processforge/config.json` and `outputs/`.
2. `pf plan` — validates schema, DOF, units, structural diff, and health-checks any Docker-containerized providers; no solver run.
3. `pf apply` — loads the last snapshot as warm start, runs the solver, and saves a new snapshot on success.

## CLI Quick Start

```bash
pf init
pf plan flowsheets/hydraulic-chain.json
pf apply flowsheets/hydraulic-chain.json
```

Use `pf plan` to validate the flowsheet and preview changes without running the solver. Use `pf apply` to solve, warm-start from the last converged state, and save a new snapshot.

## Commands

```bash
pf init [--path PATH]
pf plan flowsheets/hydraulic-chain.json [--no-diagram] [--output-dir diagrams/]
# For flowsheets using Docker providers (OpenMC, FESTIM), plan probes each
# provider's /health endpoint and fails if a container is unreachable.
# Bring containers up first with: pf init <flowsheet.json>
pf apply flowsheets/hydraulic-chain.json
pf run flowsheets/closed-loop-chain.json [--export-images]
pf diagram flowsheets/hydraulic-chain.json [--format svg] [--output-dir diagrams/]
pf export-modelica flowsheets/hydraulic-chain.json [--output-dir modelica/] [--no-compile]
pf export-fmu flowsheets/hydraulic-chain.json [--output-dir outputs/] [--backend scipy]
```

### Outputs

`pf apply` or `pf run` writes a single unified archive under `outputs/`:

- `outputs/<base>.pfarchive/` — the unified store for a flowsheet's solved state and outputs:
  - `snapshots/` — Zarr store of converged state vectors (one group per successful `pf apply`), each with `x`/`x_delta` arrays and config/var-name/metadata attributes, plus a `latest` pointer.
  - `runs/<run_id>.json` — the run manifest: every stream and unit engine output (values, units, dtypes, shapes) and run provenance (backend, version, flowsheet hash).
  - `outputs/streams/<name>.json` — per-stream timeseries from the solve.
  - `artifacts.json` — content-addressed registry of all output artifacts (local + remote URIs).
  - `index.json` — `field_name → occurrences` index for fast cross-run lookups.
  - `latest_run` — plain-text pointer to the most recent run.
- `outputs/<base>_divergence.json` — written when both direct and homotopy solves fail (on `pf apply`).

## Python API

```python
from processforge import EOFlowsheet, validate_flowsheet

config = validate_flowsheet("flowsheets/closed-loop-chain.json")
fs = EOFlowsheet(config, backend="scipy")
results = fs.run()
```

Optional backends:
- `"pyomo"` — requires `processforge[eo]`
- `"casadi"` — requires `processforge[eo-casadi]`

## Docker (Provider Images)

Individual providers such as FESTIM and OpenMC can be run via Docker images using the provider image contract. Processforge itself is not distributed as a Docker image.

Refer to the provider documentation for available Docker image tags and configuration.
