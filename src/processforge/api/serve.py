"""Processforge provider HTTP API.

Run as a module so deployment platforms that execute
``python -m processforge.api.serve`` can boot the server directly::

    python -m processforge.api.serve

The listen port is taken from the ``PORT`` environment variable
(default ``9000``), and the provider type from ``PROVIDER_TYPE``.
"""

from __future__ import annotations

import os
import traceback

from fastapi import FastAPI, HTTPException
from loguru import logger

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
        _maybe_upload_run_outputs(result)
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


def _maybe_upload_run_outputs(result) -> None:
    """Upload the provider's run_dir to S3 when configured.

    No-op unless ``S3_BUCKET`` is set in the container environment and the
    provider recorded a ``run_dir`` in its result metadata. Upload failures are
    logged and never fail the run.
    """
    run_dir = result.metadata.get("run_dir")
    if not run_dir:
        return
    try:
        from processforge.utils.s3_upload import upload_directory_to_s3

        s3_uri = upload_directory_to_s3(run_dir)
        if s3_uri:
            result.metadata["s3_uri"] = s3_uri
    except Exception as exc:  # noqa: BLE001
        logger.warning(f"S3 upload of {run_dir} failed: {exc}")


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=DEFAULT_PORT)
