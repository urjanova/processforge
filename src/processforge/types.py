"""Framework-level typed dataclasses for the provider interface.

These types form the contract between the processforge framework and providers.
They contain only fields the framework itself understands; provider-specific
properties (e.g. D_0/E_D for Festim, nuclides for OpenMC) travel in ``extra``.

Adding a new provider never requires changes to these dataclasses.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


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

    friendly_material_id: int
    density: Optional[float] = None
    density_units: Optional[str] = None
    temperature: Optional[float] = None
    depletable: bool = False
    nuclides: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    _KNOWN_FIELDS: frozenset = field(
        default=frozenset({
            "friendly_material_id", "density", "density_units",
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
            "friendly_material_id", "density", "density_units",
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
            "friendly_material_id", "density", "density_units",
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
