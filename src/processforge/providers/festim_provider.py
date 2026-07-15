"""FESTIM provider — HTTP client for the FESTIM Docker service.

FESTIM runs as an external Docker container with a FastAPI server
(``processforge-festim``).  This module is the client side: it serialises
the simulation request, sends it over HTTP, and deserialises the response.

The container's API contract (``GET /health``, ``POST /run``) is documented
in ``docs/provider-docker-images.md``.
"""
from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import TYPE_CHECKING, Optional

from loguru import logger

from .base import AbstractProvider
from .registry import (
    get_provider_default_port,
    register_provider,
)

if TYPE_CHECKING:
    from processforge.types import (
        FestimProviderConfig,
        FlowsheetConfig,
        SimulationResult,
        UnitConfig,
    )

_DEFAULT_PORT = 9002


class FestimProvider(AbstractProvider):
    """HTTP client for the FESTIM Docker service.

    Responsibilities
    ----------------
    * **Health check** — ``initialize()`` verifies the container is reachable
      via ``GET /health``.
    * **Simulation dispatch** — ``run_simulation()`` serialises the unit config
      and materials, sends ``POST /run``, and deserialises the JSON response
      into a ``SimulationResult``.
    * **Material validation** — ``validate_material()`` enforces FESTIM-specific
      required fields (``D_0``, ``E_D``) locally so flowsheet validation works
      without contacting the service.
    """

    def __init__(self):
        self._url: Optional[str] = None
        self._provider_output_dir: str = "outputs/festim"
        self._materials: dict = {}
        self._initialized: bool = False

    def _resolve_url(self, provider_config: "FestimProviderConfig") -> str:
        """Return the service URL, deriving a default from the catalog if needed."""
        url = provider_config.url
        if not url:
            port = get_provider_default_port("festim") or _DEFAULT_PORT
            url = f"http://localhost:{port}"
        return url.rstrip("/")

    def _health_check(self, url: str) -> dict:
        """GET /health — returns the JSON body on success."""
        req = urllib.request.Request(f"{url}/health", method="GET")
        with urllib.request.urlopen(req, timeout=10) as resp:
            return json.loads(resp.read().decode())

    def initialize(
        self,
        provider_config: "FestimProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        """Store config, resolve URL, verify the Docker service is reachable."""
        self._url = self._resolve_url(provider_config)

        # Resolve output dir (kept for metadata; the container manages its own output).
        out_dir = os.path.expandvars(provider_config.output_dir)
        if not os.path.isabs(out_dir):
            root = os.environ.get("PROCESSFORGE_OUTPUT_DIR", "outputs")
            out_dir = os.path.join(root, out_dir)
        self._provider_output_dir = out_dir

        # Material registry — used by validate_material() and serialised in /run.
        self._materials = {}
        for mat_name, mat_def in flowsheet_config.materials.items():
            self._materials[mat_name] = mat_def

        # Verify the service is up.
        try:
            health = self._health_check(self._url)
            logger.info(
                f"FestimProvider connected to {self._url} — "
                f"status={health.get('status')}"
            )
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            raise RuntimeError(
                f"FESTIM service unreachable at {self._url}. "
                f"Ensure the container is running: docker compose up {exc}"
            ) from exc

        self._initialized = True

    def get_thermo_properties(self, stream: dict) -> dict:
        raise NotImplementedError("FestimProvider does not support stream thermodynamics.")

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        """Return ``None`` — FESTIM uses ``run_simulation`` via ``SolverUnit``."""
        return None

    def teardown(self) -> None:
        """Release provider state."""
        self._initialized = False

    @classmethod
    def validate_material(cls, mat_name: str, mat_def, unit_cfg) -> list:
        """Validate FESTIM-specific material properties.

        Rules
        -----
        * ``extra.D_0`` is required (diffusion pre-exponential).
        * ``extra.E_D`` is required (diffusion activation energy).

        Returns:
            List of error strings (empty = valid).
        """
        errors = []
        extra = mat_def.extra or {}

        if "D_0" not in extra:
            errors.append(
                f"Material '{mat_name}' is missing 'extra.D_0' "
                "(diffusion pre-exponential, required for FESTIM)."
            )
        if "E_D" not in extra:
            errors.append(
                f"Material '{mat_name}' is missing 'extra.E_D' "
                "(diffusion activation energy, required for FESTIM)."
            )

        return errors

    def run_simulation(self, unit_config: "UnitConfig", inlet: dict) -> "SimulationResult":
        """POST /run to the FESTIM Docker service and return the result.

        Serialises ``unit_config`` and ``materials`` into the JSON body
        expected by the container's ``POST /run`` endpoint, sends the request,
        and deserialises the response into a ``SimulationResult``.
        """
        from processforge.types import SimulationResult

        if not self._url:
            raise RuntimeError("FestimProvider has not been initialized — call initialize() first.")

        body = {
            "unit_config": {
                "type": unit_config.type,
                "sim_type": unit_config.sim_type,
                "material": unit_config.material,
                "solver_config": unit_config.solver_config or {},
                **unit_config.extra,
            },
            "materials": {
                name: {
                    "id": mat.id,
                    "density": mat.density,
                    "density_units": mat.density_units,
                    "temperature": mat.temperature,
                    "depletable": mat.depletable,
                    "nuclides": mat.nuclides,
                    **mat.extra,
                }
                for name, mat in self._materials.items()
            },
            "inlet": inlet,
            "output_dir": self._provider_output_dir,
        }

        logger.info(
            f"FestimProvider: POST {self._url}/run "
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
            with urllib.request.urlopen(req, timeout=300) as resp:
                data = json.loads(resp.read().decode())
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode() if exc.fp else str(exc)
            raise RuntimeError(
                f"FESTIM service returned HTTP {exc.code}: {detail}"
            ) from exc
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            raise RuntimeError(
                f"FESTIM service unreachable at {self._url}: {exc}"
            ) from exc

        result = SimulationResult(
            status=data.get("status", "completed"),
            sim_type=data.get("sim_type", unit_config.sim_type or ""),
            scalars=data.get("scalars", {}),
            metadata=data.get("metadata", {}),
        )

        logger.info(
            f"FestimProvider: '{result.sim_type}' completed — "
            f"scalars={list(result.scalars.keys())}"
        )

        return result


# Register the provider
register_provider("festim", FestimProvider)
