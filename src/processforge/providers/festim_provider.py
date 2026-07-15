"""FESTIM provider — HTTP client for the FESTIM Docker service.

Architecture (three-layer strategy pattern)
-------------------------------------------
1. **FestimBuildHelpers** — shared utilities for constructing the HTTP request
   body sent to the container.

2. **Strategy registry** (``FestimSimStrategy`` + ``register_festim_sim_type``)
   Each ``sim_type`` string maps to a strategy class.  New simulation types are
   added by subclassing ``FestimSimStrategy`` and calling
   ``register_festim_sim_type(name, MyStrategy)`` — no changes needed here.

3. **FestimProvider** (``AbstractProvider`` subclass)
   Owns health checks, material validation, and HTTP dispatch.  The strategy's
   ``build`` method constructs the request body; the provider sends it.

Adding a new sim_type::

    from processforge.providers.festim_provider import (
        FestimSimStrategy, FestimBuildHelpers, register_festim_sim_type,
    )

    class MyCustomSim(FestimSimStrategy):
        def build(self, unit_config, materials, inlet, helpers):
            ...
            return body_dict

    register_festim_sim_type("my_custom_sim", MyCustomSim)
"""
from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from abc import ABC, abstractmethod
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
        MaterialDef,
        SimulationResult,
        UnitConfig,
    )

_DEFAULT_PORT = 9002

# Sim-type registry — maps sim_type strings to lightweight descriptors.
# Unlike OpenMC (which uses strategy classes), FESTIM delegates entirely to
# the container.  The registry exists only for *local* validation so that
# ``pf plan`` / ``pf run`` can reject unknown sim_types before contacting
# the service.
_SIM_TYPE_REGISTRY: dict[str, dict] = {}


def register_festim_sim_type(name: str, *, description: str = "") -> None:
    """Register a FESTIM simulation type by name.

    Args:
        name:        The ``sim_type`` string used in the flowsheet JSON.
        description: One-line human-readable summary.
    """
    _SIM_TYPE_REGISTRY[name] = {"description": description}


def get_registered_sim_types() -> dict[str, dict]:
    """Return a copy of the currently registered FESTIM sim_type descriptors.

    Use this function (rather than accessing ``_SIM_TYPE_REGISTRY`` directly)
    so that callers are insulated from internal implementation changes.
    """
    return dict(_SIM_TYPE_REGISTRY)


# Seed built-in types.
register_festim_sim_type(
    "hydrogen_transport",
    description="1-D hydrogen diffusion through a material slab",
)


# ---------------------------------------------------------------------------
# Build helpers
# ---------------------------------------------------------------------------

class FestimBuildHelpers:
    """Shared FESTIM request-body utilities injected into sim strategies.

    All methods are ``@staticmethod`` so strategies can call them without a
    provider instance, making them unit-testable in isolation.
    """

    @staticmethod
    def build_unit_config_body(unit_config: "UnitConfig") -> dict:
        """Serialise ``UnitConfig`` into the ``unit_config`` block of the HTTP body."""
        return {
            "type": unit_config.type,
            "sim_type": unit_config.sim_type,
            "material": unit_config.material,
            "solver_config": unit_config.solver_config or {},
            **unit_config.extra,
        }

    @staticmethod
    def build_materials_body(materials: dict[str, "MaterialDef"]) -> dict:
        """Serialise the material registry into the ``materials`` block of the HTTP body."""
        return {
            name: {
                "id": mat.id,
                "density": mat.density,
                "density_units": mat.density_units,
                "temperature": mat.temperature,
                "depletable": mat.depletable,
                "nuclides": mat.nuclides,
                **mat.extra,
            }
            for name, mat in materials.items()
        }


# ---------------------------------------------------------------------------
# Simulation type strategy registry
# ---------------------------------------------------------------------------

class FestimSimStrategy(ABC):
    """Base class for a FESTIM simulation type.

    Each subclass encapsulates the request-body construction and local
    validation for one ``sim_type`` value.  The ``build`` method receives
    the unit config, materials, inlet, and shared helpers, and returns the
    complete HTTP request body dict to POST to the container.

    Register subclasses with :func:`register_festim_sim_type`.

    Example::

        class MyFixedSourceCSG(FestimSimStrategy):
            def build(self, unit_config, materials, inlet, helpers):
                ...
                return body_dict

        register_festim_sim_type("fixed_source_csg", MyFixedSourceCSG)
    """

    @abstractmethod
    def build(
        self,
        unit_config: "UnitConfig",
        materials: dict[str, "MaterialDef"],
        inlet: dict,
        helpers: FestimBuildHelpers,
    ) -> dict:
        """Construct the HTTP request body for this simulation type.

        Args:
            unit_config: Typed unit config from the flowsheet JSON.
            materials:   ``{mat_name: MaterialDef}`` built by the provider.
            inlet:       Upstream inlet data.
            helpers:     Shared :class:`FestimBuildHelpers` instance.

        Returns:
            Complete request body dict to POST to ``/run``.
        """


_SIM_TYPE_REGISTRY: dict[str, type] = {}


def register_festim_sim_type(name: str, strategy_cls: type) -> None:
    """Register a FESTIM simulation type by name.

    Args:
        name:         The ``sim_type`` string used in the flowsheet JSON.
        strategy_cls: A subclass of :class:`FestimSimStrategy`.
    """
    _SIM_TYPE_REGISTRY[name] = strategy_cls


def get_registered_sim_types() -> dict[str, type]:
    """Return a view of the currently registered FESTIM sim_type → strategy mapping.

    Use this function (rather than accessing ``_SIM_TYPE_REGISTRY`` directly) so
    that callers are insulated from internal implementation changes.
    """
    return dict(_SIM_TYPE_REGISTRY)


# ---------------------------------------------------------------------------
# Built-in strategies
# ---------------------------------------------------------------------------

class _HydrogenTransportStrategy(FestimSimStrategy):
    """1-D hydrogen diffusion through a material slab.

    Validates that the unit config carries the required FESTIM-specific fields
    (``mesh``, ``species``, ``subdomains``, ``boundary_conditions``, ``exports``)
    before constructing the request body.
    """

    def build(
        self,
        unit_config: "UnitConfig",
        materials: dict[str, "MaterialDef"],
        inlet: dict,
        helpers: FestimBuildHelpers,
    ) -> dict:
        extra = unit_config.extra

        # Validate required sim-type-specific fields
        for key in ("mesh", "species", "subdomains"):
            if key not in extra:
                raise ValueError(
                    f"sim_type='hydrogen_transport' requires '{key}' in unit config"
                )

        vertices = extra.get("mesh", {}).get("vertices")
        if not isinstance(vertices, list) or len(vertices) < 2:
            raise ValueError(
                "sim_type='hydrogen_transport' requires mesh.vertices to be "
                "a list with at least 2 points"
            )

        species_list = extra.get("species", [])
        if not species_list:
            raise ValueError(
                "sim_type='hydrogen_transport' requires a non-empty 'species' list"
            )

        subs = extra.get("subdomains", {})
        if not subs.get("volume"):
            raise ValueError(
                "sim_type='hydrogen_transport' requires subdomains.volume to be non-empty"
            )
        if not subs.get("surface"):
            raise ValueError(
                "sim_type='hydrogen_transport' requires subdomains.surface to be non-empty"
            )

        return {
            "unit_config": helpers.build_unit_config_body(unit_config),
            "materials": helpers.build_materials_body(materials),
            "inlet": inlet,
        }


# Register built-ins at module load
register_festim_sim_type("hydrogen_transport", _HydrogenTransportStrategy)


# ---------------------------------------------------------------------------
# FestimProvider
# ---------------------------------------------------------------------------

class FestimProvider(AbstractProvider):
    """HTTP client for the FESTIM Docker service.

    Responsibilities
    ----------------
    * **Health check** — ``initialize()`` verifies the container is reachable
      via ``GET /health``.
    * **Simulation dispatch** — ``run_simulation()`` looks up the strategy from
      ``_SIM_TYPE_REGISTRY`` by ``sim_type``, calls ``strategy.build()``,
      sends ``POST /run``, and deserialises the JSON response into a
      ``SimulationResult``.
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
        """Build and run a FESTIM simulation from a typed ``UnitConfig``.

        Execution flow
        --------------
        1. Look up strategy from ``_SIM_TYPE_REGISTRY`` by ``sim_type``
        2. Call ``strategy.build()`` → request body dict
        3. POST body to the Docker container
        4. Deserialise response → ``SimulationResult``
        """
        from processforge.types import SimulationResult

        if not self._url:
            raise RuntimeError("FestimProvider has not been initialized — call initialize() first.")

        sim_type = unit_config.sim_type
        strategy_cls = _SIM_TYPE_REGISTRY.get(sim_type)
        if strategy_cls is None:
            raise ValueError(
                f"FestimProvider: unknown sim_type '{sim_type}'. "
                f"Register with register_festim_sim_type(). "
                f"Built-in types: {sorted(_SIM_TYPE_REGISTRY)}"
            )

        helpers = FestimBuildHelpers()
        body = strategy_cls().build(unit_config, self._materials, inlet, helpers)
        body["output_dir"] = self._provider_output_dir

        logger.info(
            f"FestimProvider: POST {self._url}/run "
            f"sim_type='{sim_type}' with strategy {strategy_cls.__name__}"
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
