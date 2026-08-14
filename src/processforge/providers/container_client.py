"""CLI-side HTTP client for containerized providers (OpenMC, FESTIM, …).

A provider that runs in a Docker container does its real work inside the
container (in ``processforge.api.serve``).  The CLI side must NOT import the
backend library (OpenMC, FESTIM, …) locally — instead it speaks HTTP to the
container's ``provider_server.py`` instance.  :class:`ContainerProviderClient`
implements the :class:`~processforge.providers.base.AbstractProvider` contract
as a thin HTTP wrapper around that server.

Request/response contract (matches ``processforge.api.serve``)
----------------------------------------------------------------
* ``POST /run`` with a JSON body
  ``{"unit_config", "materials", "inlet", "output_dir", "provider_config"}``.
* JSON response ``{"status", "sim_type", <scalars…>, "metadata"}`` deserialised
  into a :class:`~processforge.types.SimulationResult`.
"""
from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import TYPE_CHECKING, Optional

from loguru import logger

from .base import AbstractProvider
from .registry import get_provider_default_port

if TYPE_CHECKING:
    from processforge.types import (
        FlowsheetConfig,
        MaterialDef,
        ProviderConfig,
        SimulationResult,
        UnitConfig,
    )


class ContainerProviderClient(AbstractProvider):
    """CLI-side HTTP client for a provider running in a Docker container.

    Responsibilities
    ----------------
    * **Health check** — ``initialize()`` verifies the container is reachable
      via ``GET /health``.
    * **Simulation dispatch** — ``run_simulation()`` serialises the unit config
      and material registry and ``POST``s them to the container's ``/run``
      endpoint, then deserialises the JSON response into a ``SimulationResult``.
    * **Material validation** — delegated to the provider class via
      ``validate_material()`` (it is a classmethod and needs no local backend).
    """

    def __init__(self, provider_type: str):
        self._ptype = provider_type
        self._url: Optional[str] = None
        self._provider_output_dir: str = "outputs"
        self._provider_config: Optional["ProviderConfig"] = None
        self._materials: dict = {}
        self._initialized: bool = False

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_url(self, provider_config: "ProviderConfig") -> str:
        """Return the service URL, deriving a default from the catalog if needed."""
        url = getattr(provider_config, "url", None)
        if not url:
            port = get_provider_default_port(self._ptype) or 9000
            url = f"http://localhost:{port}"
        return url.rstrip("/")

    @staticmethod
    def _serialize_unit_config(unit_config: "UnitConfig") -> dict:
        """Serialise a ``UnitConfig`` into the ``unit_config`` block of the body."""
        return {
            "type": unit_config.type,
            "in": unit_config.inputs,
            "material": unit_config.material,
            "provider": unit_config.provider,
            "out": unit_config.out,
            "sim_type": unit_config.sim_type,
            "solver_config": unit_config.solver_config,
            **unit_config.extra,
        }

    @staticmethod
    def _serialize_materials(materials: dict[str, "MaterialDef"]) -> dict:
        """Serialise the material registry into the ``materials`` block of the body."""
        out: dict = {}
        for name, mat in materials.items():
            out[name] = {
                "id": mat.id,
                "density": mat.density,
                "density_units": mat.density_units,
                "temperature": mat.temperature,
                "depletable": mat.depletable,
                "nuclides": mat.nuclides,
                **mat.extra,
            }
        return out

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def initialize(
        self,
        provider_config: "ProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        """Store config, resolve URL, verify the Docker service is reachable."""
        self._url = self._resolve_url(provider_config)

        out_dir = os.path.expandvars(
            getattr(provider_config, "output_dir", "outputs") or "outputs"
        )
        if not os.path.isabs(out_dir):
            root = os.environ.get("PROCESSFORGE_OUTPUT_DIR", "outputs")
            out_dir = os.path.join(root, out_dir)
        self._provider_output_dir = out_dir

        # Material registry — serialised into /run and used by validate_material().
        self._materials = dict(flowsheet_config.materials)
        self._provider_config = provider_config

        try:
            req = urllib.request.Request(f"{self._url}/health", method="GET")
            with urllib.request.urlopen(req, timeout=10) as resp:
                health = json.loads(resp.read().decode())
            logger.info(
                f"ContainerProviderClient connected to {self._url} — "
                f"status={health.get('status')}"
            )
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            raise RuntimeError(
                f"Provider '{self._ptype}' service unreachable at {self._url}. "
                f"Ensure the container is running: pf init <flowsheet>"
            ) from exc

        self._initialized = True

    def get_thermo_properties(self, stream: dict) -> dict:
        raise NotImplementedError(
            "ContainerProviderClient does not support stream thermodynamics."
        )

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        """Return ``None`` — the container provider uses ``run_simulation``."""
        return None

    def teardown(self) -> None:
        """Release provider state."""
        self._initialized = False

    def run_simulation(self, unit_config: "UnitConfig", inlet: dict) -> "SimulationResult":
        """Serialise the request, POST it to the container, return the result."""
        from processforge.types import SimulationResult

        if not self._url:
            raise RuntimeError(
                "ContainerProviderClient has not been initialized — call initialize() first."
            )

        provider_config_dump = (
            self._provider_config.model_dump()
            if self._provider_config is not None
            else {"type": unit_config.provider or self._ptype, "url": self._url}
        )

        body = {
            "unit_config": self._serialize_unit_config(unit_config),
            "materials": self._serialize_materials(self._materials),
            "inlet": inlet,
            "output_dir": self._provider_output_dir,
            "provider_config": provider_config_dump,
        }

        logger.info(
            f"ContainerProviderClient: POST {self._url}/run "
            f"sim_type='{unit_config.sim_type}'"
        )

        payload = json.dumps(body).encode()
        req = urllib.request.Request(
            f"{self._url}/run",
            data=payload,
            headers={"Content-Type": "application/json"},
            method="POST",
        )

        try:
            with urllib.request.urlopen(req, timeout=900) as resp:
                data = json.loads(resp.read().decode())
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode() if exc.fp else str(exc)
            raise RuntimeError(
                f"Provider '{self._ptype}' service returned HTTP {exc.code}: {detail}"
            ) from exc
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            raise RuntimeError(
                f"Provider '{self._ptype}' service unreachable at {self._url}: {exc}"
            ) from exc

        scalars = {
            k: v
            for k, v in data.items()
            if k not in ("status", "sim_type", "metadata")
        }
        return SimulationResult(
            status=data.get("status", "completed"),
            sim_type=data.get("sim_type", unit_config.sim_type or ""),
            scalars=scalars,
            metadata=data.get("metadata", {}),
        )


