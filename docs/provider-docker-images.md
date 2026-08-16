# Running Processforge Providers

Three common setups:

---

## Scenarios at a Glance

| Scenario | You run locally | Provider runs on | What's happening (plain English) |
|----------|-----------------|------------------|----------------------------------|
| **1. Local + OpenMC Docker** | `pf` CLI + `docker compose` | Your machine (Docker) | `pf init` writes a `docker-compose.yml` and pulls the image. You start it with `docker compose up -d`. `pf run` talks to the already-running container, saves results to `outputs/`. |
| **2. Local CLI + Railway Docker** | `pf` CLI | Railway cloud (Docker) | You run `pf run` on your laptop. It sends the job to your Railway URL over HTTPS. Railway runs Docker, computes, sends results back to your `outputs/`. |
| **3. Railway + FastAPI** | `pf` CLI | Railway cloud (FastAPI, no Docker for you) | Same as #2, but Railway runs your Python/FastAPI code directly (no Dockerfile needed). You just push code to GitHub → Railway. |

---

## Scenario 1: Local + OpenMC Docker (Recommended for Development)

**Prerequisites:** Docker (with `docker compose`) installed, `pf` CLI installed (`pip install processforge`)

```bash
# 1. Create flowsheet
pf init my_flowsheet.json

# 2. Edit my_flowsheet.json (pf init creates a template with OpenMC)
#    It already has:
#    "providers": {
#      "openmc": {
#        "type": "openmc",
#        "docker_image": "ghcr.io/urjanova/processforge-openmc:latest",
#        "url": "http://localhost:9001"
#      }
#    }

# 3. Initialise the provider environment (generates docker-compose.yml + pulls image)
pf init my_flowsheet.json

# 4. Start the container(s) — this is you, not pf! Containers keep running
#    until you stop them.
docker compose -f .processforge/<env>/docker-compose.yml up -d

# 5. Run (talks to the already-running container on http://localhost:9001)
pf run my_flowsheet.json
```

**What happens:**

- `pf init <flowsheet.json>` reads the `providers` block and, for every
  *containerized, locally-addressed* provider (type `openmc`/`festim` with no
  URL or a `localhost`/`127.0.0.1` URL), writes a `docker-compose.yml` into a
  per-flowsheet environment directory:

  ```
  .processforge/<flowsheet-hash>/docker-compose.yml
  ```

  The hash is derived from the flowsheet's absolute path, so each flowsheet keeps
  its own provider environment. It then runs `docker compose pull` to fetch the
  images. `pf init` does **not** start the containers.
- You start (and stop) the containers yourself with `docker compose`:

  ```bash
  # List the generated services
  docker compose -f .processforge/<env>/docker-compose.yml config --services

  # Start in the background
  docker compose -f .processforge/<env>/docker-compose.yml up -d

  # Follow logs
  docker compose -f .processforge/<env>/docker-compose.yml logs -f

  # Stop (keeps images/volumes)
  docker compose -f .processforge/<env>/docker-compose.yml down
  ```

  The compose file maps the provider's API port (e.g. `9001:9001` for OpenMC,
  `9002:9002` for FESTIM) and mounts a host directory at `/data` inside the
  container:

  ```yaml
  services:
    openmc:
      image: ghcr.io/urjanova/processforge-openmc:latest
      ports:
        - "9001:9001"
      volumes:
        - "${OPENMC_DATA_ROOT:-outputs}:/data"
  ```

  Override the mount with `OPENMC_DATA_ROOT` to point the container at your
  OpenMC data/workspace instead of `outputs`:

  ```bash
  OPENMC_DATA_ROOT=/path/to/data docker compose -f .processforge/<env>/docker-compose.yml up -d
  ```

- `pf run` (and `pf apply`) does **not** start or stop containers. It assumes
  the container is already running and pings `GET {url}/health`. If a
  containerized provider is unreachable, `pf run` fails fast with:

  ```
  Provider 'openmc' unreachable at http://localhost:9001. Run: pf init <flowsheet>
  ```

  Start the containers (step 4 above) and re-run. To verify availability without
  running a full sim, use `pf validate <flowsheet>` or `pf plan <flowsheet>`.

### How `docker-compose.yml` is generated for local providers

`pf init` is what creates the compose file used to run providers on your
machine. The generation rules (in `src/processforge/compose.py` and
`src/processforge/cli/init.py`):

- **Which providers get a service:** `pf init` walks the flowsheet's
  `providers` block. A provider gets a compose entry only when it is both:

  1. *Containerized* — its `type` has a `docker_image` in the provider
     registry (`openmc`, `festim`, …), and
  2. *Locally addressed* — it has no `url`, or a `localhost` / `127.0.0.1` /
     `0.0.0.0` URL (see `is_local_provider_url`).

  For each such provider, `pf init` records its `docker_image`,
  `default_port`, and `url` and asks the compose generator to write a service.

- **Which providers are skipped:** any containerized provider with a *remote*
  URL (e.g. a Railway `https://…` endpoint) is treated as externally managed.
  `pf init` logs `remote — skipping compose` and writes **no** service for it —
  this is why Scenarios 2 and 3 need no compose file. Non-containerized
  (pip-installable) providers like `coolprop`/`cantera` are also skipped.

- **Where it's written:** a per-flowsheet environment directory keyed by a
  stable hash of the flowsheet's absolute path:

  ```
  .processforge/<flowsheet-hash>/docker-compose.yml
  ```

  (`flowsheet_hash` comes from `lock.flowsheet_env_dir`, so two different
  flowsheets never clobber each other's provider environments.)

- **What each service entry contains** (from `compose.generate_compose`):

  ```yaml
  services:
    <provider_name>:
      image: <docker_image>            # e.g. ghcr.io/urjanova/processforge-openmc:latest
      ports:
        - "<port>:<port>"              # host:container, e.g. 9001:9001
      volumes:
        - "${OPENMC_DATA_ROOT:-outputs}:/data"
  ```

  The file is stamped with an "Auto-generated by pf init — do not edit
  manually" header; re-run `pf init <flowsheet.json>` to regenerate it.

- **Pull + (no) start:** after writing the file, `pf init` runs
  `docker compose pull` to fetch the images. It never starts or stops
  containers — that's your `docker compose up -d` / `down` step above. If
  Docker isn't installed it warns and you pull/run manually.

- **Re-generation / migration:** re-running `pf init` rewrites the compose
  file (e.g. after you change `docker_image`), warning if the environment
  already exists. A legacy single-environment layout
  (`.processforge/docker-compose.yml` + `lock.json` at the root) is migrated
  into the per-flowsheet directory automatically on first `pf init`.

**Use your own OpenMC image:**
```bash
docker build -t my-openmc -f docker/Dockerfile.openmc .
```
Then change `docker_image` in the flowsheet to `my-openmc:latest` and re-run
`pf init my_flowsheet.json` to regenerate the compose file (it will reference
your image instead of the `ghcr.io` one).

### What the provider images do

The prebuilt `ghcr.io/urjanova/processforge-openmc:latest` and
`ghcr.io/urjanova/processforge-festim:latest` images — and the reference
[`docker/Dockerfile.openmc`](../docker/Dockerfile.openmc) and
[`docker/Dockerfile.festim`](../docker/Dockerfile.festim) you can build from —
share the same contract that the generated compose file relies on:

- Both `FROM mambaorg/micromamba:1.5.8`, `pip install "processforge[api]"`, and
  serve the provider API via `python -m processforge.api.serve`.
- **OpenMC** ([`docker/Dockerfile.openmc`](../docker/Dockerfile.openmc)): pins
  `openmc=0.15.3`, sets `PORT=9001` / `EXPOSE 9001`, `OPENMC_DATA_ROOT=/data`,
  `PROVIDER_TYPE=openmc`, `VOLUME /data`, and runs `fetch_openmc_data.sh` as its
  `ENTRYPOINT` to populate the cross-section data mounted at `/data`.
- **FESTIM** ([`docker/Dockerfile.festim`](../docker/Dockerfile.festim)):
  installs `festim fenics-dolfinx` from conda-forge, sets `PORT=9002` /
  `EXPOSE 9002`, `PROVIDER_TYPE=festim`, `VOLUME /data`.

These `PORT`/`EXPOSE` values are exactly the host:container mappings emitted in
the compose `ports:` block, so `http://localhost:9001` (OpenMC) and
`http://localhost:9002` (FESTIM) reach the running containers.

---

## Scenario 2: Local CLI + Your Docker on Railway

**Prerequisites:** Railway account, `pf` CLI locally, Docker image pushed to a registry (GHCR, Docker Hub)

### One-time setup: Push your provider to Railway
```bash
# 1. Build and push your image
docker build -t ghcr.io/youruser/processforge-my-provider:latest .
docker push ghcr.io/youruser/processforge-my-provider:latest

# 2. Create Railway service
#    - New Project → Deploy from Image
#    - Image: ghcr.io/youruser/processforge-my-provider:latest
#    - Port: 9003 (or your provider's port)
#    - Railway gives you: https://my-provider.up.railway.app
```

### Run from your laptop
```bash
# 1. Create flowsheet
pf init my_flowsheet.json

# 2. Edit to point at Railway (NO docker_image field)
{
  "providers": {
    "my_provider": {
      "type": "my_provider",
      "url": "https://my-provider.up.railway.app"
    }
  },
  "units": { ... }
}

# 3. Run
pf run my_flowsheet.json
```

**What happens:** Your laptop sends HTTP to Railway. Railway's Docker container runs the sim. Results download to your `outputs/`.

---

## Scenario 3: Railway + FastAPI (No Docker on Your Side)

**Prerequisites:** Railway account, `pf` CLI locally, provider code on GitHub

### One-time setup: Deploy code to Railway
```bash
# 1. Push provider code to GitHub (with pyproject.toml, main.py, etc.)
git push origin main

# 2. Railway: New Project → GitHub Repo
#    - Railway auto-detects Python, runs `pip install -e .`
#    - Start command: uvicorn my_provider:app --host 0.0.0.0 --port $PORT
#    - Railway gives you: https://my-provider.up.railway.app
```

### Run from your laptop
```bash
pf init my_flowsheet.json
# Edit flowsheet - same as Scenario 2, just URL, no docker_image
pf run my_flowsheet.json
```

**What happens:** Identical to Scenario 2, but Railway runs your Python process directly (no Docker layer). Simpler for pure-Python providers.

---

## Provider API Contract (All Scenarios)
Your provider must implement:

| Endpoint | Method | Request | Response |
|----------|--------|---------|----------|
| `/health` | GET | — | `{"status": "ready", "provider_type": "openmc"}` |
| `/run` | POST | UnitConfig + Materials | Scalars + metadata |

Full spec: [provider-api.openapi.json](../provider-api.openapi.json)

---

## Quick Test Any Provider
```bash
# Health
curl https://your-provider-url/health

# Test run (replace with valid request for your provider)
curl -X POST https://your-provider-url/run \
  -H "Content-Type: application/json" \
  -d '{"unit_config": {"type": "SolverUnit", "sim_type": "..."}, "materials": {}, "inlet": {}}'
```

---

## Registering a New Provider Type
Add to `src/processforge/providers/registry.py`:
```python
"my_provider": {
    "module": "processforge.providers.my_provider",
    "class": "MyProvider",
    "optional_dep": "my_provider",
    "description": "My custom provider",
    "docker_image": "ghcr.io/you/processforge-my-provider:latest",  # used in Scenario 1
    "default_port": 9003,
},
```
- `docker_image` = used by `pf init` for **Scenario 1** (local Docker)
- `default_port` = fallback if `url` omitted in flowsheet