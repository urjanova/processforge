"""Framework-level typed dataclasses for the provider interface.

These types form the contract between the processforge framework and providers.
They contain only fields the framework itself understands; provider-specific
properties (e.g. D_0/E_D for Festim, nuclides for OpenMC) travel in ``extra``.

Adding a new provider never requires changes to these dataclasses.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Union


@dataclass
class MaterialDef:
    """Framework-level representation of a material loaded from the flowsheet JSON.

    Fields mirror the schema's ``materials`` section. Provider-specific
    properties (D_0/E_D for Festim, nuclides for OpenMC, …) are captured
    in ``extra`` so providers can access them without the framework needing
    to know about them.

    Usage::

        mat = MaterialDef.from_dict(config["materials"]["tungsten"])
        diffusivity = mat.extra.get("D_0")          # Festim field
        mat.get("D_0")                               # unified accessor
    """

    id: int
    density: Optional[float] = None
    density_units: Optional[str] = None
    temperature: Optional[float] = None
    depletable: bool = False
    nuclides: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    _KNOWN_FIELDS: frozenset = field(
        default=frozenset({
            "id", "density", "density_units",
            "temperature", "depletable", "nuclides",
        }),
        init=False,
        repr=False,
        compare=False,
    )

    @classmethod
    def from_dict(cls, d: dict) -> "MaterialDef":
        """Construct from a raw dict (as loaded from JSON).

        Known framework fields populate typed attributes; everything else
        goes into ``extra`` for provider-specific access.
        """
        known_keys = frozenset({
            "id", "density", "density_units",
            "temperature", "depletable", "nuclides",
        })
        known = {k: v for k, v in d.items() if k in known_keys}
        extra = {k: v for k, v in d.items() if k not in known_keys}
        return cls(**known, extra=extra)

    def get(self, key: str, default=None):
        """Unified accessor — checks typed fields first, then ``extra``.

        Lets provider code use ``mat_def.get("D_0")`` without caring
        whether the field is a typed attribute or a provider-specific extra.
        """
        known_keys = frozenset({
            "id", "density", "density_units",
            "temperature", "depletable", "nuclides",
        })
        if key in known_keys:
            return getattr(self, key, default)
        return self.extra.get(key, default)


@dataclass
class UnitConfig:
    """Framework-level representation of a unit definition loaded from JSON.

    ``in`` is a Python keyword so inlet stream names use ``inputs``.
    Provider-specific unit fields beyond the framework schema go into ``extra``.

    Usage::

        ucfg = UnitConfig.from_dict(config["units"]["membrane"])
        ucfg.type          # "FestimMembrane"
        ucfg.sim_type      # "heat_2d"
        ucfg.solver_config # opaque dict passed to provider
    """

    type: str
    material: int
    provider: Optional[str] = None
    inputs: list = field(default_factory=list)      # JSON "in" key
    out: Optional[str] = None
    retentate_out: Optional[str] = None
    permeate_out: Optional[str] = None
    sim_type: Optional[str] = None
    solver_config: dict = field(default_factory=dict)
    extra: dict = field(default_factory=dict)

    _KNOWN_FIELDS: frozenset = field(
        default=frozenset({
            "type", "material", "provider", "out", "retentate_out",
            "permeate_out", "sim_type", "solver_config",
        }),
        init=False,
        repr=False,
        compare=False,
    )

    @classmethod
    def from_dict(cls, d: dict) -> "UnitConfig":
        """Construct from a raw unit config dict."""
        known_keys = frozenset({
            "type", "material", "provider", "out", "retentate_out",
            "permeate_out", "sim_type", "solver_config",
        })
        raw_in = d.get("in", [])
        inputs = [raw_in] if isinstance(raw_in, str) else list(raw_in or [])
        known = {k: v for k, v in d.items() if k in known_keys}
        extra = {k: v for k, v in d.items() if k not in known_keys and k != "in"}
        return cls(inputs=inputs, **known, extra=extra)


# ---------------------------------------------------------------------------
# Per-provider configuration dataclasses
# ---------------------------------------------------------------------------

@dataclass
class CoolPropProviderConfig:
    """Configuration for the built-in CoolProp provider (no extra fields)."""

    type: str = "coolprop"

    @classmethod
    def from_dict(cls, d: dict) -> "CoolPropProviderConfig":
        return cls()


@dataclass
class CanteraProviderConfig:
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


@dataclass
class FestimProviderConfig:
    """Configuration for the FESTIM hydrogen-transport provider.

    ``materials`` holds inline material overrides that take priority over the
    flowsheet's global ``materials`` section.

    Flowsheet JSON example::

        "providers": {
            "festim": {"type": "festim", "output_dir": "outputs/festim_2d"}
        }
    """

    type: str = "festim"
    output_dir: str = "output"
    materials: dict = field(default_factory=dict)

    @classmethod
    def from_dict(cls, d: dict) -> "FestimProviderConfig":
        return cls(
            output_dir=d.get("output_dir", "output"),
            materials=d.get("materials", {}),
        )


@dataclass
class ModelicaProviderConfig:
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


@dataclass
class OpenMCProviderConfig:
    """Configuration for the OpenMC Monte Carlo neutronics provider.

    Flowsheet JSON example::

        "providers": {
            "openmc": {
                "type": "openmc",
                "output_dir": "outputs/openmc",
                "cross_sections": "/path/to/cross_sections.xml"
            }
        }
    """

    type: str = "openmc"
    output_dir: str = "outputs/openmc"
    cross_sections: Optional[str] = None

    @classmethod
    def from_dict(cls, d: dict) -> "OpenMCProviderConfig":
        return cls(
            output_dir=d.get("output_dir", "outputs/openmc"),
            cross_sections=d.get("cross_sections"),
        )


#: Union of all known provider config types — extend when adding a new provider.
ProviderConfig = Union[
    CoolPropProviderConfig,
    CanteraProviderConfig,
    FestimProviderConfig,
    ModelicaProviderConfig,
    OpenMCProviderConfig,
]

_PROVIDER_CONFIG_REGISTRY: dict[str, type] = {
    "coolprop": CoolPropProviderConfig,
    "cantera": CanteraProviderConfig,
    "festim": FestimProviderConfig,
    "modelica": ModelicaProviderConfig,
    "openmc": OpenMCProviderConfig,
}


def provider_config_from_dict(d: dict) -> ProviderConfig:
    """Parse a raw provider config dict into the appropriate typed dataclass.

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

@dataclass
class FlowsheetConfig:
    """Typed representation of a complete flowsheet configuration.

    ``from_dict()`` converts the raw JSON dict into typed nested dataclasses:
    providers → per-provider config classes, units → ``UnitConfig``,
    materials → ``MaterialDef``.  Unknown top-level keys land in ``extra``
    (e.g. ``_config_path`` injected at runtime by the Modelica runner).

    Usage::

        cfg = FlowsheetConfig.from_dict(raw)
        cfg.default_provider          # Optional[str]
        cfg.materials["tungsten"]     # MaterialDef
        cfg.providers["festim"]       # FestimProviderConfig
    """

    providers: dict = field(default_factory=dict)   # dict[str, ProviderConfig]
    default_provider: Optional[str] = None
    streams: dict = field(default_factory=dict)
    units: dict = field(default_factory=dict)        # dict[str, UnitConfig]
    materials: dict = field(default_factory=dict)    # dict[str, MaterialDef]
    material_mixes: dict = field(default_factory=dict)
    simulation: dict = field(default_factory=dict)
    metadata: Optional[dict] = None
    extra: dict = field(default_factory=dict)        # runtime / unknown fields

    _KNOWN_FIELDS: frozenset = field(
        default=frozenset({
            "providers", "default_provider", "streams", "units",
            "materials", "material_mixes", "simulation", "metadata",
        }),
        init=False,
        repr=False,
        compare=False,
    )

    @classmethod
    def from_dict(cls, d: dict) -> "FlowsheetConfig":
        """Construct from a raw flowsheet config dict (as loaded from JSON)."""
        known_keys = frozenset({
            "providers", "default_provider", "streams", "units",
            "materials", "material_mixes", "simulation", "metadata",
        })
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
        pre-date this dataclass (e.g. ``_derive_model_name`` in transpiler.py)
        keep working unchanged.
        """
        known_keys = frozenset({
            "providers", "default_provider", "streams", "units",
            "materials", "material_mixes", "simulation", "metadata",
        })
        if key in known_keys:
            val = getattr(self, key)
            return val if val is not None else default
        return self.extra.get(key, default)


@dataclass
class SimulationResult:
    """Result returned by ``provider.run_simulation()``.

    ``scalars`` holds named scalar outputs (hydrogen_flux, k_effective, …)
    that downstream units or result writers can consume.
    ``metadata`` holds non-scalar information such as output file paths.
    """

    status: str                                  # "completed" | "failed"
    sim_type: str
    scalars: dict = field(default_factory=dict)  # {name: float}
    metadata: dict = field(default_factory=dict) # e.g. {"xdmf_files": [...]}

    def as_dict(self) -> dict:
        """Flat dict for flowsheet result storage."""
        return {"status": self.status, "sim_type": self.sim_type, **self.scalars}


@dataclass
class RunInfo:
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
    pkg_versions: dict[str, str] = field(default_factory=dict)
    x0: Optional[list[float]] = None
    var_names: Optional[list[str]] = None

    def as_dict(self) -> dict:
        """Dict view for compatibility with older call sites."""
        return {
            "git_hash": self.git_hash,
            "timestamp": self.timestamp,
            "python_version": self.python_version,
            "platform": self.platform,
            "processforge_version": self.processforge_version,
            "mode": self.mode,
            "backend": self.backend,
            "pkg_versions": dict(self.pkg_versions),
            "x0": list(self.x0) if self.x0 is not None else None,
            "var_names": list(self.var_names) if self.var_names is not None else None,
        }


@dataclass
class SnapshotState:
    """Typed snapshot payload loaded from a .pfstate store."""

    config: dict
    x: list[float]
    var_names: list[str] = field(default_factory=list)
    snapshot_id: str = ""
    timestamp: str = ""

    def as_dict(self) -> dict:
        """Dict view for compatibility with older call sites."""
        return {
            "config": self.config,
            "x": list(self.x),
            "var_names": list(self.var_names),
            "snapshot_id": self.snapshot_id,
            "timestamp": self.timestamp,
        }
