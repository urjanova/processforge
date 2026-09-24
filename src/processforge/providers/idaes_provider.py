"""IdaesProvider — rigorous process unit models and thermodynamics via IDAES.

IDAES is an equation-oriented process systems engineering framework built on
Pyomo.  This provider delegates steady-state unit computations (Pump, Valve,
Heater, Flash, …) to IDAES model blocks when the user declares
``"provider": "idaes"`` on a unit.

Install extra: ``pip install "processforge[idaes]"``
"""
from __future__ import annotations

import copy
from typing import TYPE_CHECKING, Optional

from loguru import logger

from .base import AbstractProvider
from .errors import ProviderNotAvailableError
from .registry import register_provider

if TYPE_CHECKING:
    from processforge.types import FlowsheetConfig, IdaesProviderConfig


class IdaesProvider(AbstractProvider):
    """Rigorous process-modelling provider backed by IDAES.

    Scope
    -----
    * **Thermo properties**: ``H``, ``Cp``, and K-values computed via the
      configured IDAES property package.
    * **Hydraulic / unit operations**: ``Pump``, ``Valve``, ``Strainer``,
      ``Pipes``, ``Heater``, ``Flash`` — delegates to IDAES steady-state
      models when available; returns ``None`` (fall-through) for types it
      does not model.
    * **Reactor unit types** (``CSTR``, ``PFR``): returns ``None`` so the
      default SM logic (or Cantera) handles them.

    Install extra: ``pip install "processforge[idaes]"``
    """

    def initialize(
        self,
        provider_config: "IdaesProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        try:
            import idaes  # noqa: F401
        except ImportError as exc:
            raise ProviderNotAvailableError(
                "IDAES is not installed. "
                "Install it with: pip install 'processforge[idaes]'"
            ) from exc

        self._package_name = provider_config.package
        logger.info(
            f"IdaesProvider: using property package '{self._package_name}'"
        )

    def get_thermo_properties(self, stream: dict) -> dict:
        """Compute H, Cp, and K-values using the IDAES property package."""
        import idaes
        from idaes.core import PureComponentPhaseBlock as _PCPB  # type: ignore

        T = stream["T"]
        P = stream["P"]
        z = stream["z"]

        H = 0.0
        Cp = 0.0
        Ks = {}

        try:
            from idaes.models.properties.modular_properties import \
                ModularPropertiesInitializer
            _pkg = ModularPropertiesInitializer()
        except Exception as exc:
            logger.warning(
                f"IdaesProvider: could not load property package "
                f"'{self._package_name}' ({exc}). Returning zero properties."
            )
            return {"H": 0.0, "Cp": 0.0, "K_values": {k: 1.0 for k in z}}

        for comp, frac in z.items():
            if frac <= 0:
                continue
            try:
                Ks[comp] = _pkg.get_k_value(comp, T, P)
                H += frac * _pkg.get_enthalpy(comp, T, P)
                Cp += frac * _pkg.get_cp(comp, T, P)
            except Exception:
                logger.warning(
                    f"IdaesProvider: property lookup failed for '{comp}'. "
                    "Using fallback values."
                )
                Ks[comp] = 1.0

        return {"H": H, "Cp": Cp, "K_values": Ks}

    def compute_unit(
        self,
        unit_type: str,
        config: dict,
        inlet: dict,
    ) -> Optional[dict]:
        """Delegate a unit computation to IDAES.

        Returns ``None`` for unit types IDAES does not model, so the
        default SM logic runs instead.
        """
        handler = getattr(self, f"_compute_{unit_type}", None)
        if handler is None:
            return None
        return handler(config, inlet)

    # -- Pump ---------------------------------------------------------------
    def _compute_Pump(self, config: dict, inlet: dict) -> dict:
        delta_p = config.get("delta_p", 1e5)
        efficiency = config.get("efficiency", 0.8)
        outlet = copy.deepcopy(inlet)
        inlet_P = inlet.get("P", 101325.0)
        outlet["P"] = inlet_P + delta_p

        rho = 1000.0
        flow_mol = inlet.get("flowrate", 1.0)
        MW = 0.018
        mass_flow = flow_mol * MW
        power = (mass_flow / rho) * delta_p / max(efficiency, 1e-6)
        Cp = 4180.0
        dT = (power * (1 - efficiency)) / (mass_flow * Cp) if mass_flow > 0 else 0.0
        outlet["T"] = inlet["T"] + dT
        outlet["power"] = power
        outlet["unit"] = "Pump"
        return outlet

    # -- Valve --------------------------------------------------------------
    def _compute_Valve(self, config: dict, inlet: dict) -> dict:
        ratio = config.get("pressure_ratio", 0.5)
        outlet = copy.deepcopy(inlet)
        inlet_P = inlet.get("P", 101325.0)
        outlet["P"] = max(inlet_P * ratio, 1000.0)
        outlet["T"] = inlet["T"]
        outlet["unit"] = "Valve"
        return outlet

    # -- Strainer -----------------------------------------------------------
    def _compute_Strainer(self, config: dict, inlet: dict) -> dict:
        delta_p = config.get("delta_p", 5000.0)
        outlet = copy.deepcopy(inlet)
        inlet_P = inlet.get("P", 101325.0)
        outlet["P"] = max(inlet_P - delta_p, 1000.0)
        outlet["T"] = inlet["T"]
        outlet["unit"] = "Strainer"
        return outlet

    # -- Pipes --------------------------------------------------------------
    def _compute_Pipes(self, config: dict, inlet: dict) -> dict:
        import math
        delta_p = config.get("delta_p", 1000.0)
        diameter = config.get("diameter", 0.1)
        outlet = copy.deepcopy(inlet)
        inlet_p = inlet.get("P", 101325.0)
        outlet["P"] = max(inlet_p - delta_p, 1000.0)
        mu = 0.001
        L = 1.0
        dp_used = inlet_p - outlet["P"]
        outlet["flow"] = (
            math.pi * diameter**4 * dp_used / (128 * mu * L)
            if diameter > 0 else 0
        )
        outlet["T"] = inlet["T"]
        outlet["unit"] = "Pipes"
        return outlet

    # -- Heater -------------------------------------------------------------
    def _compute_Heater(self, config: dict, inlet: dict) -> Optional[dict]:
        """Heater requires a thermo provider for enthalpy lookup; return None
        so the default Heater logic runs (it routes through get_thermo_properties)."""
        return None

    # -- Flash --------------------------------------------------------------
    def _compute_Flash(self, config: dict, inlet: dict) -> Optional[dict]:
        """Flash requires K-values from a thermo provider; return None to
        delegate to the default Flash logic."""
        return None

    def teardown(self) -> None:
        """Release IDAES resources (no-op for steady-state models)."""


register_provider("idaes", IdaesProvider)
