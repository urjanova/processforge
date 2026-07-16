# Running Processforge Providers

Three common setups:

---

## Scenarios at a Glance

| Scenario | You run locally | Provider runs on | What's happening (plain English) |
|----------|-----------------|------------------|----------------------------------|
| **1. Local + OpenMC Docker** | `pf` CLI | Your machine (Docker) | You run `pf run`, it starts an OpenMC container on your laptop, runs the sim, stops it. Results land in `outputs/`. |
| **2. Local CLI + Railway Docker** | `pf` CLI | Railway cloud (Docker) | You run `pf run` on your laptop. It sends the job to your Railway URL over HTTPS. Railway runs Docker, computes, sends results back to your `outputs/`. |
| **3. Railway + FastAPI** | `pf` CLI | Railway cloud (FastAPI, no Docker for you) | Same as #2, but Railway runs your Python/FastAPI code directly (no Dockerfile needed). You just push code to GitHub → Railway. |

---

## Scenario 1: Local + OpenMC Docker (Recommended for Development)

**Prerequisites:** Docker installed, `pf` CLI installed (`pip install processforge`)

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

# 3. Run
pf run my_flowsheet.json
```

**What happens:** `pf run` generates `docker-compose.yml`, starts the OpenMC container on port 9001, POSTs your simulation, saves results to `outputs/`, stops the container.

**Use your own OpenMC image:**
```bash
docker build -t my-openmc -f docker/Dockerfile.openmc .
```
Then change `docker_image` in flowsheet to `my-openmc:latest`.

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