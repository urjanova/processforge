# Provider Docker Images

Guide for building Docker-based Processforge providers.

## Why Docker?

Some simulation codes are not available on PyPI. OpenMC and FESTIM, for
example, are distributed through conda-forge and require system-level
dependencies (MPI compilers, FEniCSx, HDF5) that are difficult to install
alongside a pip-based Python environment. Docker solves this by packaging
each provider and its dependencies into an isolated, reproducible image.

Docker is also the deployment model for cloud and web service scenarios.
A provider container runs independently, exposes a REST API, and can be
scaled, updated, or replaced without touching the `pf` CLI or other
providers.

Even for packages that are available on pip, Docker can be useful when:

- The provider has heavy native dependencies (MPI, FEniCSx, PETSc).
- You need a specific version of a system library pinned.
- You want to isolate the provider's environment from the host.
- You are deploying to a platform that expects containerised services
  (Kubernetes, Railway, AWS ECS).

## How it works

The `pf` CLI and any Processforge-compatible client communicate with
provider containers over HTTP. Each container runs a thin FastAPI server
that wraps the provider's Python library and exposes a standard API.

```
┌────────────┐        HTTP         ┌──────────────────────┐
│  pf CLI    │ ──────────────────► │  processforge-openmc  │
│  (client)  │                     │  :9001                │
└────────────┘                     └──────────────────────┘
```

The flowsheet JSON declares where each provider lives:

```json
{
  "providers": {
    "openmc": {
      "type": "openmc",
      "url": "http://localhost:9001"
    }
  }
}
```

The same contract works whether the client is the `pf` CLI on your
laptop, a FastAPI app on Railway, or any other HTTP client.

## Deployment scenarios

### Local development

`pf init` generates a `docker-compose.yml` and pulls images. Containers
run locally. The flowsheet points to `localhost`:

```json
{
  "providers": {
    "openmc": {
      "type": "openmc",
      "url": "http://localhost:9001"
    }
  }
}
```

```
$ pf init flowsheet.json
$ pf validate flowsheet.json
$ pf run flowsheet.json
```

### Cloud (Railway, ECS, etc.)

Provider containers are deployed separately. The flowsheet points to
their public or internal URLs. No local Docker needed:

```json
{
  "providers": {
    "openmc": {
      "type": "openmc",
      "url": "https://openmc-production.up.railway.app"
    }
  }
}
```

```
$ pf init flowsheet.json    # writes lock file, no compose generated
$ pf validate flowsheet.json  # checks reachability
$ pf run flowsheet.json     # sends POST /run to Railway URL
```

### FastAPI app (same platform as providers)

A FastAPI app deployed on the same platform as the provider containers
can call them via internal networking. The API contract is identical —
the app sends the same `POST /run` request that `pf` sends:

```python
import requests

resp = requests.post(
    "http://openmc-internal:9001/run",
    json={"unit_config": {...}, "materials": {...}},
)
result = resp.json()
```

The `pf` CLI and the FastAPI app are both HTTP clients to the same
provider services. The contract does not change.

## API contract

A provider container must expose two endpoints.

### `GET /health`

Health check. Returns 200 when the provider is ready to accept requests.

```
GET /health

200 OK
{
  "status": "ready",
  "provider_type": "openmc",
  "version": "0.14.0"
}
```

### `POST /run`

Run a simulation. The request contains everything the container needs —
unit configuration, materials, and any provider-specific settings. The
response contains scalar results and metadata.

```
POST /run
Content-Type: application/json

{
  "unit_config": {
    "type": "SolverUnit",
    "sim_type": "eigenvalue_csg",
    "material": 1,
    "solver_config": {
      "batches": 100,
      "inactive": 10,
      "particles": 10000
    }
  },
  "materials": {
    "tungsten": {
      "id": 1,
      "density": 19.3,
      "density_units": "g/cm3",
      "nuclides": [{"name": "W184", "percent": 100.0, "percent_type": "ao"}]
    }
  }
}

200 OK
{
  "status": "completed",
  "sim_type": "eigenvalue_csg",
  "scalars": {
    "k_eff": 1.0032,
    "k_eff_std_dev": 0.0008
  },
  "metadata": {
    "run_dir": "/data/openmc",
    "statepoint_path": "/data/openmc/statepoint.100.h5"
  }
}
```

The container is stateless. Every request is self-contained. The `pf`
CLI handles all validation, initialization, and teardown. The container
just receives a request and returns a result.

## Response format

All endpoints return JSON. Successful responses use HTTP 200. Errors
use HTTP 4xx/5xx with a consistent body:

```json
{
  "detail": "OpenMC cross sections not found at /data/xs/cross_sections.xml"
}
```

## Command support

| Command | Containerized providers | Pip providers |
|---------|------------------------|---------------|
| `pf plan` | Yes (validation only, local) | Yes |
| `pf validate` | Yes (checks reachability) | Yes (checks importability) |
| `pf run` | Yes (POST /run to container) | Yes (runs in-process) |
| `pf apply` | No — not supported | Yes (coolprop, cantera, modelica) |

`pf apply` uses the EO solver which requires in-process access to
thermodynamic property calls. It only works with pip-installable
providers. Use `pf run` for containerized providers.

## Creating a provider Dockerfile

### Minimal template

```dockerfile
FROM python:3.12-slim

WORKDIR /app

# Install provider dependencies.
# For conda-forge packages, use micromamba or conda.
# For pip-only packages, use pip directly.
RUN pip install --no-cache-dir \
    processforge \
    <provider-dependencies>

# Copy the provider server (see "Provider server" below).
COPY provider_server.py .

EXPOSE 9001

CMD ["python", "provider_server.py"]
```

### Example: OpenMC

```dockerfile
FROM mambaorg/micromamba:1.5.8

USER root
RUN apt-get update && apt-get install -y --no-install-recommends build-essential && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# OpenMC is only available on conda-forge
RUN micromamba install -y -n base -c conda-forge openmc && \
    micromamba clean --all --yes

WORKDIR /app
RUN micromamba run -n base pip install processforge fastapi uvicorn

COPY provider_server.py .

ENV OPENMC_DATA_ROOT=/data
VOLUME /data

EXPOSE 9001
CMD ["micromamba", "run", "-n", "base", "python", "provider_server.py"]
```

### Example: FESTIM

```dockerfile
FROM mambaorg/micromamba:1.5.8

USER root
RUN apt-get update && apt-get install -y --no-install-recommends build-essential && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# FESTIM + FEniCSx are only available on conda-forge
RUN micromamba install -y -n base -c conda-forge festim fenics-dolfinx && \
    micromamba clean --all --yes

WORKDIR /app
RUN micromamba run -n base pip install processforge fastapi uvicorn

COPY provider_server.py .

VOLUME /data

EXPOSE 9002
CMD ["micromamba", "run", "-n", "base", "python", "provider_server.py"]
```

### Example: custom provider from pip

Even if your provider is available on pip, you may want a Docker image
for isolation or deployment:

```dockerfile
FROM python:3.12-slim

WORKDIR /app
RUN pip install --no-cache-dir processforge my-custom-provider fastapi uvicorn

COPY provider_server.py .

EXPOSE 9003
CMD ["python", "provider_server.py"]
```

## Provider server

Each container runs a thin FastAPI server with two endpoints.

```python
"""Thin FastAPI server wrapping a Processforge provider."""
from fastapi import FastAPI, HTTPException

app = FastAPI()


@app.get("/health")
def health():
    return {
        "status": "ready",
        "provider_type": "openmc",
    }


@app.post("/run")
def run(body: dict):
    from processforge.types import UnitConfig, MaterialDef

    unit_cfg = UnitConfig.from_dict(body["unit_config"])
    materials = {
        name: MaterialDef.from_dict(mat)
        for name, mat in body.get("materials", {}).items()
    }

    # Initialize provider
    from processforge.providers.registry import get_provider_class
    provider = get_provider_class("openmc")()
    # ... call provider.initialize() with config from body ...

    # Run simulation
    try:
        result = provider.run_simulation(unit_cfg, body.get("inlet", {}))
        return result.as_dict() | {"metadata": result.metadata}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        provider.teardown()


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=9001)
```

## Volume mounts

Provider containers typically need access to:

| Mount | Container path | Purpose |
|-------|---------------|---------|
| Output directory | `/data` | Simulation results, statepoint files |
| Cross-section data | `/data/cross_sections` | OpenMC nuclear data |

The `pf init` command generates compose files with these mounts
configured. The host path defaults to `outputs/` and can be overridden
with `PROCESSFORGE_OUTPUT_DIR`.

For cloud deployments (Railway, ECS), use platform-native volume
attachments instead of Docker Compose mounts.

## Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `PROCESSFORGE_OUTPUT_DIR` | `outputs` | Host directory mounted as `/data` |
| `OPENMC_DATA_ROOT` | `/data` | Root directory for OpenMC data |
| `OPENMC_DATA_URL` | (none) | URL to download cross-section archive on first start |

## Port conventions

| Provider | Default port |
|----------|-------------|
| openmc | 9001 |
| festim | 9002 |
| (future) | 9003+ |

## Registering your provider

After building and pushing your Docker image, register it in the
Processforge provider catalog so `pf init` knows about it.

**File:** `src/processforge/providers/registry.py`

Add an entry to `_PROVIDER_CATALOG`:

```python
"my_provider": {
    "module": "processforge.providers.my_provider",
    "class": "MyProvider",
    "optional_dep": "my_provider",
    "description": "Description of what this provider does",
    "docker_image": "ghcr.io/yourorg/processforge-my-provider:latest",
    "default_port": 9003,
},
```

The `docker_image` field tells `pf init` to generate a Docker Compose
service. The `default_port` is used when the flowsheet omits the `url`
field.

## Flowsheet declaration

Users reference your provider in their flowsheet:

```json
{
  "providers": {
    "my_provider": {
      "type": "my_provider",
      "url": "http://localhost:9003"
    }
  }
}
```

If `url` is omitted, `pf init` defaults to `http://localhost:{default_port}`.

## Testing your container

1. Build the image:
   ```
   docker build -t processforge-my-provider .
   ```

2. Run the container:
   ```
   docker run -p 9003:9003 processforge-my-provider
   ```

3. Test the health endpoint:
   ```
   curl http://localhost:9003/health
   ```

4. Test simulation:
   ```
   curl -X POST http://localhost:9003/run \
     -H "Content-Type: application/json" \
     -d '{"unit_config": {"type": "SolverUnit", "sim_type": "..."}, "materials": {}, "inlet": {}}'
   ```

5. Validate from the `pf` CLI:
   ```
   pf validate my_flowsheet.json
   ```
