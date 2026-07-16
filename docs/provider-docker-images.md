# Provider Docker Images

Guide for building Docker-based Processforge providers.

## Provider Checklist
- [ ] Dockerfile builds & runs
- [ ] `GET /health` returns `{"status": "ready", "provider_type": "..."}`
- [ ] `POST /run` accepts UnitConfig + Materials, returns scalars + metadata
- [ ] Errors return 4xx/5xx with `{"detail": "..."}`
- [ ] Registered in `registry.py` with `docker_image` and `default_port`
- [ ] Cross-section/data volumes documented (if needed)

## Quick Test
```bash
docker build -t my-provider .
docker run -p 9003:9003 my-provider
curl localhost:9003/health
curl -X POST localhost:9003/run -d @test-request.json
pf validate flowsheet.json  # with your image in providers.my_provider.docker_image
```

## API Contract
See [provider-api.openapi.json](../provider-api.openapi.json) for the full OpenAPI 3.1 specification.

## Deployment Scenarios

| Scenario | `docker_image` | `url` | `pf init` behavior |
|----------|----------------|-------|-------------------|
| Local (catalog default) | omitted | omitted | Generates compose, uses catalog image & `localhost:{default_port}` |
| Local (custom image) | your image | omitted | Generates compose, uses your image & `localhost:{default_port}` |
| Cloud (Railway/ECS) | your image | your URL | Writes lock file only, no compose |
| Same-platform FastAPI | omitted | internal URL | No Docker needed |

## Minimal Dockerfile Template
```dockerfile
FROM python:3.12-slim
WORKDIR /app
COPY pyproject.toml ./
RUN pip install --no-cache-dir "processforge[fastapi]"
COPY my_provider.py .
EXPOSE 9003
CMD ["uvicorn", "my_provider:app", "--host", "0.0.0.0", "--port", "9003"]
```

Full examples: [`docker/Dockerfile.openmc`](../docker/Dockerfile.openmc), [`docker/Dockerfile.festim`](../docker/Dockerfile.festim)

## Volumes

| Container path | Purpose | Host default |
|----------------|---------|--------------|
| `/data` | Simulation outputs | `outputs/` (`PROCESSFORGE_OUTPUT_DIR`) |
| `/data/cross_sections` | OpenMC nuclear data | `outputs/cross_sections` |

Generated compose files configure these mounts. Override with `PROCESSFORGE_OUTPUT_DIR`.

## Environment Variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `PROCESSFORGE_OUTPUT_DIR` | `outputs` | Host directory mounted as `/data` |
| `OPENMC_DATA_ROOT` | `/data` | Root for OpenMC data |
| `OPENMC_DATA_URL` | (none) | URL to download cross-section archive on first start |

## Port Conventions

| Provider | Default port |
|----------|-------------|
| openmc | 9001 |
| festim | 9002 |
| (future) | 9003+ |

## Registering Your Provider

Edit `src/processforge/providers/registry.py`:
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

## Flowsheet Declaration
```json
{
  "providers": {
    "my_provider": {
      "type": "my_provider",
      "docker_image": "ghcr.io/yourorg/processforge-my-provider:latest",
      "url": "http://localhost:9003"
    }
  }
}
```
- `docker_image` omitted → catalog default
- `url` omitted → `http://localhost:{default_port}`