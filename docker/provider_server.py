"""Thin FastAPI server wrapping a Processforge provider.

Each provider container runs this server.  The provider type is resolved
from the ``PROVIDER_TYPE`` environment variable (set in the Dockerfile).

Endpoints
---------
* ``GET /health``  — readiness probe, returns provider metadata.
* ``POST /run``    — run a simulation, return scalar results.
"""
from __future__ import annotations

import os
import traceback

from fastapi import FastAPI, HTTPException

app = FastAPI(title="Processforge Provider")

PROVIDER_TYPE = os.environ.get("PROVIDER_TYPE", "")
DEFAULT_PORT = int(os.environ.get("PORT", "9000"))


@app.get("/health")
def health():
    return {
        "status": "ready",
        "provider_type": PROVIDER_TYPE,
    }


@app.post("/run")
def run(body: dict):
    from processforge.types import UnitConfig, MaterialDef

    try:
        unit_cfg = UnitConfig.from_dict(body["unit_config"])
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Invalid unit_config: {exc}") from exc

    materials_raw = body.get("materials", {})
    materials = {
        name: MaterialDef.from_dict(mat)
        for name, mat in materials_raw.items()
    }

    from processforge.providers.registry import get_provider_class

    try:
        provider_cls = get_provider_class(PROVIDER_TYPE)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    provider = provider_cls()

    # Build the provider config for initialize(). Prefer the real config the
    # CLI serialised into the request body (carries cross_sections, output_dir,
    # url, …); fall back to a minimal container-default config otherwise.
    provider_config_raw = body.get("provider_config")
    if provider_config_raw:
        from processforge.types import provider_config_from_dict

        provider_config = provider_config_from_dict(provider_config_raw)
    else:
        provider_config = type("Cfg", (), {
            "url": None,
            "output_dir": os.environ.get("PROCESSFORGE_OUTPUT_DIR", "/data"),
            "cross_sections": None,
            "type": PROVIDER_TYPE,
        })()

    # Build a minimal flowsheet config holding materials.
    flowsheet_config = type("FS", (), {"materials": materials})()

    try:
        provider.initialize(provider_config, flowsheet_config)
        result = provider.run_simulation(unit_cfg, body.get("inlet", {}))
        return result.as_dict() | {"metadata": result.metadata}
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Simulation failed: {exc}\n{traceback.format_exc()}",
        ) from exc
    finally:
        provider.teardown()


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=DEFAULT_PORT)
