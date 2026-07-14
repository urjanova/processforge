# Usage

## Plan / Apply Workflow

1. `pf init` — creates `.processforge/config.json` and `outputs/`.
2. `pf plan` — validates schema, DOF, units, and structural diff; no solver run.
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
pf apply flowsheets/hydraulic-chain.json
pf run flowsheets/closed-loop-chain.json [--export-images]
pf diagram flowsheets/hydraulic-chain.json [--format svg] [--output-dir diagrams/]
pf export-modelica flowsheets/hydraulic-chain.json [--output-dir modelica/] [--no-compile]
pf export-fmu flowsheets/hydraulic-chain.json [--output-dir outputs/] [--backend scipy]
```

### Outputs

`pf apply` or `pf run` creates:
- `*_results.zarr` — simulation results store.
- `*_validation.xlsx` — validation report.
- `*.pfstate/` — versioned snapshot store with `latest` pointer.
- `*_divergence.json` — written when both direct and homotopy solves fail.

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

## Cloud and Docker

The repository supports containerized execution and S3-compatible upload.

### Run in Docker

```bash
docker run --rm \
  -v "$(pwd)/flowsheets:/app/flowsheets" \
  ghcr.io/urjanova/processforge:latest \
  pf apply /app/flowsheets/my-flowsheet.json
```

#### Getting outputs onto the host

Run outputs (the results `.zarr` and provider artifacts such as OpenMC statepoint/XML
files) are written under `PROCESSFORGE_OUTPUT_DIR`, which defaults to `/data` in the image.
Bind-mount a host directory at `/data` to collect them:

```bash
docker run --rm \
  -v "$(pwd)/tmp_files:/data" \
  ghcr.io/urjanova/processforge:latest \
  processforge run flowsheets/openmc/msre_eigenvalue.json
# → results appear in ./tmp_files (e.g. msre_eigenvalue_results.zarr, openmc/msre_run/…)
```

Override the location with `-e PROCESSFORGE_OUTPUT_DIR=<path>` (e.g. a subdirectory of the
mount so each run is isolated). A relative `output_dir` in a flowsheet's provider config is
resolved against this root.

### OpenMC in containers

OpenMC workflows require nuclear cross-section data. The container startup script `scripts/fetch_openmc_data.sh` downloads and caches it automatically into `/data` (the same volume also collects run outputs).

### API server

`pf-serve` launches an HTTP API for programmatic submission.

```bash
docker run -d --name pf-api -p 9000:9000 ghcr.io/urjanova/processforge:latest pf-serve
```
