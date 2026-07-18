"""Framework-level typed models for the provider interface.

These types form the contract between the processforge framework and providers.
They contain only fields the framework itself understands; provider-specific
properties (e.g. nuclides for OpenMC) travel in ``extra``.

Adding a new provider never requires changes to these models.
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
from pydantic import BaseModel, ConfigDict, Field


class MaterialDef(BaseModel):
    """Framework-level representation of a material loaded from the flowsheet JSON.

    Fields mirror the schema's ``materials`` section. Provider-specific
    properties (e.g. nuclides for OpenMC) are captured in ``extra`` so
    providers can access them without the framework needing to know about them.

    Usage::

        mat = MaterialDef.from_dict(config["materials"]["tungsten"])
        mat.get("D_0")                               # unified accessor
    """

    id: int
    density: Optional[float] = None
    density_units: Optional[str] = None
    temperature: Optional[float] = None
    depletable: bool = False
    nuclides: list = []
    extra: dict = {}

    @classmethod
    def from_dict(cls, d: dict) -> "MaterialDef":
        """Construct from a raw dict (as loaded from JSON).

        Known framework fields populate typed attributes; everything else
        goes into ``extra`` for provider-specific access.
        """
        known_keys = cls.model_fields.keys() - {"extra"}
        known = {k: v for k, v in d.items() if k in known_keys}
        extra = {k: v for k, v in d.items() if k not in known_keys}
        # An explicit "extra" object in the JSON holds provider-specific
        # properties; flatten it so providers read e.g. extra["D_0"] directly
        # rather than extra["extra"]["D_0"].
        nested = extra.pop("extra", None)
        if isinstance(nested, dict):
            extra = {**nested, **extra}
        return cls(**known, extra=extra)

    def get(self, key: str, default=None):
        """Unified accessor — checks typed fields first, then ``extra``.

        Lets provider code use ``mat_def.get("D_0")`` without caring
        whether the field is a typed attribute or a provider-specific extra.
        """
        known_keys = type(self).model_fields.keys() - {"extra"}
        if key in known_keys:
            return getattr(self, key, default)
        return self.extra.get(key, default)


class UnitConfig(BaseModel):
    """Framework-level representation of a unit definition loaded from JSON.

    ``in`` is a Python keyword so inlet stream names use ``inputs``.
    Provider-specific unit fields beyond the framework schema go into ``extra``.

    Usage::

        ucfg = UnitConfig.from_dict(config["units"]["membrane"])
        ucfg.type          # "SolverUnit"
        ucfg.sim_type      # "heat_2d"
        ucfg.solver_config # opaque dict passed to provider
    """

    type: str
    material: Optional[int] = None
    provider: Optional[str] = None
    inputs: list = []                                # JSON "in" key
    out: Optional[str] = None
    retentate_out: Optional[str] = None
    permeate_out: Optional[str] = None
    sim_type: Optional[str] = None
    solver_config: dict = {}
    extra: dict = {}

    @classmethod
    def from_dict(cls, d: dict) -> "UnitConfig":
        """Construct from a raw unit config dict."""
        raw_in = d.get("in", [])
        inputs = [raw_in] if isinstance(raw_in, str) else list(raw_in or [])
        known_keys = cls.model_fields.keys() - {"extra", "inputs"}
        known = {k: v for k, v in d.items() if k in known_keys}
        extra = {k: v for k, v in d.items() if k not in known_keys and k != "in"}
        return cls(inputs=inputs, **known, extra=extra)


# ---------------------------------------------------------------------------
# Per-provider configuration models
# ---------------------------------------------------------------------------


class CoolPropProviderConfig(BaseModel):
    """Configuration for the built-in CoolProp provider (no extra fields)."""

    type: str = "coolprop"

    @classmethod
    def from_dict(cls, d: dict) -> "CoolPropProviderConfig":
        return cls()


class CanteraProviderConfig(BaseModel):
    """Configuration for the Cantera thermochemistry provider.

    Flowsheet JSON example::

        "providers": {
            "cantera": {"type": "cantera", "mechanism": "gri30.yaml", "phase": null}
        }
    """

    type: str = "cantera"
    mechanism: str = "gri30.yaml"
    phase: Optional[str] = None

    @classmethod
    def from_dict(cls, d: dict) -> "CanteraProviderConfig":
        return cls(
            mechanism=d.get("mechanism", "gri30.yaml"),
            phase=d.get("phase") or None,
        )


class ModelicaProviderConfig(BaseModel):
    """Configuration for the OpenModelica FMU provider.

    Flowsheet JSON example::

        "providers": {
            "modelica": {"type": "modelica", "output_dir": "outputs"}
        }
    """

    type: str = "modelica"
    output_dir: str = "outputs"

    @classmethod
    def from_dict(cls, d: dict) -> "ModelicaProviderConfig":
        return cls(output_dir=d.get("output_dir", "outputs"))


class OpenMCProviderConfig(BaseModel):
    """Configuration for the OpenMC Monte Carlo neutronics provider.

    Paths support environment variable expansion (``${VAR}`` / ``$VAR``), which
    is the recommended way to manage data files in Docker and Railway deployments.

    Flowsheet JSON example::

        "providers": {
            "openmc": {
                "type": "openmc",
                "output_dir": "outputs/openmc",
                "cross_sections": "${OPENMC_DATA_ROOT}/cross_sections/cross_sections.xml"
            }
        }

    Custom Docker image::

        "providers": {
            "openmc": {
                "type": "openmc",
                "docker_image": "my-org/custom-openmc:v2"
            }
        }

    Container deployment pattern:

    * Set ``OPENMC_DATA_ROOT`` to the mount point of your data volume (default ``/data``).
    * Set ``OPENMC_DATA_URL`` to download a cross-section ``.tar.gz`` archive on first start.
    * ``OPENMC_CROSS_SECTIONS`` is exported automatically by the startup script — do not
      set it manually when using ``OPENMC_DATA_URL``.

    Local usage::

        OPENMC_DATA_ROOT=/path/to/openmc_data processforge run flowsheet.json
    """

    type: str = "openmc"
    output_dir: str = "outputs/openmc"
    cross_sections: Optional[str] = None
    docker_image: Optional[str] = None

    @classmethod
    def from_dict(cls, d: dict) -> "OpenMCProviderConfig":
        return cls(
            output_dir=d.get("output_dir", "outputs/openmc"),
            cross_sections=d.get("cross_sections"),
            docker_image=d.get("docker_image"),
        )


class FestimProviderConfig(BaseModel):
    """Configuration for the FESTIM hydrogen transport provider.

    FESTIM runs as an external Docker service. The ``url`` field points to
    the service's HTTP endpoint (default ``http://localhost:9002``).

    Flowsheet JSON examples::

        "providers": {
            "festim": {
                "type": "festim",
                "url": "http://localhost:9002"
            }
        }

    Custom Docker image::

        "providers": {
            "festim": {
                "type": "festim",
                "docker_image": "my-org/custom-festim:v2"
            }
        }

    Cloud deployment::

        "providers": {
            "festim": {
                "type": "festim",
                "url": "https://festim-production.up.railway.app"
            }
        }
    """

    type: str = "festim"
    url: Optional[str] = None
    output_dir: str = "outputs/festim"
    docker_image: Optional[str] = None

    @classmethod
    def from_dict(cls, d: dict) -> "FestimProviderConfig":
        return cls(
            url=d.get("url"),
            output_dir=d.get("output_dir", "outputs/festim"),
            docker_image=d.get("docker_image"),
        )


#: Union of all known provider config types — extend when adding a new provider.
ProviderConfig = Union[
    CoolPropProviderConfig,
    CanteraProviderConfig,
    ModelicaProviderConfig,
    OpenMCProviderConfig,
    FestimProviderConfig,
]

_PROVIDER_CONFIG_REGISTRY: dict[str, type] = {
    "coolprop": CoolPropProviderConfig,
    "cantera": CanteraProviderConfig,
    "modelica": ModelicaProviderConfig,
    "openmc": OpenMCProviderConfig,
    "festim": FestimProviderConfig,
}


def provider_config_from_dict(d: dict) -> ProviderConfig:
    """Parse a raw provider config dict into the appropriate typed model.

    Args:
        d: Raw provider config block from the flowsheet JSON
           (must contain a ``"type"`` key).

    Raises:
        ValueError: If ``"type"`` is missing or names an unknown provider.
    """
    ptype = d.get("type")
    if not ptype:
        raise ValueError("Provider config is missing the required 'type' field.")
    cls = _PROVIDER_CONFIG_REGISTRY.get(ptype)
    if cls is None:
        raise ValueError(
            f"Unknown provider type '{ptype}'. "
            f"Known types: {sorted(_PROVIDER_CONFIG_REGISTRY)}"
        )
    return cls.from_dict(d)


# ---------------------------------------------------------------------------
# Full flowsheet configuration
# ---------------------------------------------------------------------------


class FlowsheetConfig(BaseModel):
    """Typed representation of a complete flowsheet configuration.

    ``from_dict()`` converts the raw JSON dict into typed nested models:
    providers → per-provider config classes, units → ``UnitConfig``,
    materials → ``MaterialDef``.  Unknown top-level keys land in ``extra``
    (e.g. ``_config_path`` injected at runtime by the Modelica runner).

    Usage::

        cfg = FlowsheetConfig.from_dict(raw)
        cfg.default_provider          # Optional[str]
        cfg.materials["tungsten"]     # MaterialDef
        cfg.providers["openmc"]       # OpenMCProviderConfig
    """

    providers: dict = {}             # dict[str, ProviderConfig]
    default_provider: Optional[str] = None
    streams: dict = {}
    units: dict = {}                 # dict[str, UnitConfig]
    materials: dict = {}             # dict[str, MaterialDef]
    material_mixes: dict = {}
    simulation: dict = {}
    metadata: Optional[dict] = None
    extra: dict = {}                 # runtime / unknown fields

    @classmethod
    def from_dict(cls, d: dict) -> "FlowsheetConfig":
        """Construct from a raw flowsheet config dict (as loaded from JSON)."""
        known_keys = cls.model_fields.keys() - {"extra"}
        providers = {
            name: provider_config_from_dict(cfg)
            for name, cfg in d.get("providers", {}).items()
        }
        units = {
            name: UnitConfig.from_dict(cfg)
            for name, cfg in d.get("units", {}).items()
        }
        materials = {
            name: MaterialDef.from_dict(mat)
            for name, mat in d.get("materials", {}).items()
        }
        extra = {k: v for k, v in d.items() if k not in known_keys}
        return cls(
            providers=providers,
            default_provider=d.get("default_provider"),
            streams=d.get("streams", {}),
            units=units,
            materials=materials,
            material_mixes=d.get("material_mixes", {}),
            simulation=d.get("simulation", {}),
            metadata=d.get("metadata"),
            extra=extra,
        )

    def get(self, key: str, default=None):
        """Dict-compatible accessor for backwards-compatible code paths.

        Checks typed fields first, then ``extra``.  Returns ``default`` when
        the key is absent from both.  Exists primarily so helper functions that
        pre-date this model (e.g. ``_derive_model_name`` in transpiler.py)
        keep working unchanged.
        """
        known_keys = type(self).model_fields.keys() - {"extra"}
        if key in known_keys:
            val = getattr(self, key)
            return val if val is not None else default
        return self.extra.get(key, default)


class SimulationResult(BaseModel):
    """Result returned by ``provider.run_simulation()``.

    ``scalars`` holds named scalar outputs (hydrogen_flux, k_effective, …)
    that downstream units or result writers can consume.
    ``metadata`` holds non-scalar information such as output file paths.
    """

    status: str                                  # "completed" | "failed"
    sim_type: str
    scalars: dict = {}                           # {name: float}
    metadata: dict = {}                          # e.g. {"xdmf_files": [...]}

    def as_dict(self) -> dict:
        """Flat dict for flowsheet result storage."""
        return {"status": self.status, "sim_type": self.sim_type, **self.scalars}


class RunInfo(BaseModel):
    """Provenance metadata for a simulation run.

    This structure is written to the Zarr ``run_info`` group by
    ``save_results_zarr`` and is produced by ``build_run_info``.
    """

    git_hash: str
    timestamp: str
    python_version: str
    platform: str
    processforge_version: str
    mode: str
    backend: str
    pkg_versions: dict[str, str] = {}
    x0: Optional[list[float]] = None
    var_names: Optional[list[str]] = None

    def as_dict(self) -> dict:
        """Dict view for compatibility with older call sites."""
        return self.model_dump()


class SnapshotState(BaseModel):
    """Typed snapshot payload loaded from a .pfstate store."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    config: dict
    x: list[float]
    var_names: list[str] = []
    snapshot_id: str = ""
    timestamp: str = ""
    metadata: dict = {}
    x_delta: np.ndarray = Field(default_factory=lambda: np.array([]))
    parent_snapshot_id: str | None = None

    def as_dict(self) -> dict:
        """Dict view for compatibility with older call sites."""
        return {
            "config": self.config,
            "x": list(self.x),
            "var_names": list(self.var_names),
            "snapshot_id": self.snapshot_id,
            "timestamp": self.timestamp,
            "metadata": self.metadata,
            "x_delta": self.x_delta.tolist() if len(self.x_delta) else [],
            "parent_snapshot_id": self.parent_snapshot_id,
        }


class StreamTimeseries(BaseModel):
    """Typed representation of a single stream's simulation timeseries.

    Each stream carries scalar property lists (T, P, flowrate) and a
    component-mole-fraction dict ``z``.  ``extra='allow'`` permits units to
    attach additional properties (e.g. ``phase``, ``beta``) without breaking
    the model.
    """

    model_config = ConfigDict(extra="allow")

    time: list[float] = []
    T: list[float] = []
    P: list[float] = []
    flowrate: list[float] = []
    z: dict[str, list[float]] = {}


class MergedInletTimeseries(BaseModel):
    """Result of flow-weighted merging of multiple inlet stream timeseries.

    Carries the same shape as :class:`StreamTimeseries` but is produced
    exclusively by ``Flowsheet._merge_inlet_timeseries``.  Dict-like
    access (``[]``, ``.get()``, ``.items()``, ``in``) is supported so
    downstream code that previously consumed plain dicts works unchanged.
    """

    time: list[float]
    T: list[float]
    P: list[float]
    flowrate: list[float]
    z: dict[str, list[float]]

    # -- dict-compatible helpers for downstream consumers ------------------

    def get(self, key: str, default=None):
        return getattr(self, key) if key in MergedInletTimeseries.model_fields else default

    def items(self):
        return ((name, getattr(self, name)) for name in MergedInletTimeseries.model_fields)

    def keys(self):
        return MergedInletTimeseries.model_fields.keys()

    def __contains__(self, key: str) -> bool:
        return key in MergedInletTimeseries.model_fields

    def __iter__(self):
        return iter(MergedInletTimeseries.model_fields)
