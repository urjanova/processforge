"""CanteraProvider — thermochemistry and reactor kinetics via Cantera."""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from loguru import logger

from .base import AbstractProvider
from .cantera_jacobian import CanteraJacobianMixin
from .registry import register_provider

if TYPE_CHECKING:
    from processforge.types import CanteraProviderConfig, FlowsheetConfig

# Unit types this provider can compute natively.
_REACTOR_UNIT_TYPES = {"CSTR", "PFR", "IdealGasReactor"}


class CanteraProvider(AbstractProvider, CanteraJacobianMixin):
    """Provides thermochemical properties and reactor simulations via Cantera.

    Scope
    -----
    * **Thermo properties**: replaces CoolProp for ``H``, ``Cp``, and K-values
      using a Cantera ``Solution`` loaded from a mechanism file.
    * **Reactor unit types**: ``CSTR``, ``PFR``, ``IdealGasReactor`` — runs a
      ``ct.ReactorNet`` for the given residence time and returns the outlet
      state.
    * **Standard unit types** (``Pump``, ``Valve``, etc.): returns ``None``
      from ``compute_unit`` so the default SM logic runs unchanged.

    Cantera is not symbolically differentiable, so this provider is
    **incompatible with the Pyomo and CasADi EO backends**.  Use
    ``backend="scipy"`` (the default) when combining Cantera with the EO solver.

    Flowsheet JSON declaration::

        "providers": {
            "cantera": {
                "type": "cantera",
                "mechanism": "gri30.yaml",
                "phase": null
            }
        }

    Install extra: ``pip install "processforge[cantera]"``
    """

    # Exposed so the SM factory can discover and register reactor unit classes.
    unit_types: dict[str, type] = {}

    def initialize(
        self,
        provider_config: "CanteraProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        try:
            import cantera as ct
        except ImportError as exc:
            raise RuntimeError(
                "Cantera is not installed. "
                "Install it with: pip install 'processforge[cantera]'"
            ) from exc

        self._ct = ct
        mechanism = provider_config.mechanism
        phase = provider_config.phase
        if phase:
            self._gas = ct.Solution(mechanism, phase)
        else:
            self._gas = ct.Solution(mechanism)
        logger.info(
            f"CanteraProvider: loaded mechanism '{mechanism}' "
            f"({len(self._gas.species_names)} species)"
        )

        # Populate unit_types lazily (avoids importing unit modules at provider
        # import time, which could cause circular deps).
        if not CanteraProvider.unit_types:
            from processforge.units.cstr import CSTR
            from processforge.units.pfr import PFR
            CanteraProvider.unit_types = {"CSTR": CSTR, "PFR": PFR, "IdealGasReactor": CSTR}

    def get_thermo_properties(self, stream: dict) -> dict:
        """Compute H, Cp, and K-values using the loaded Cantera mechanism."""
        gas = self._gas
        T = stream["T"]
        P = stream["P"]
        z = stream["z"]

        # Set gas state using mole fractions.
        composition = ", ".join(f"{k}:{v}" for k, v in z.items() if v > 0)
        try:
            gas.TPX = T, P, composition
        except Exception as exc:
            logger.warning(
                f"CanteraProvider: failed to set gas state ({exc}). "
                "Returning zero properties."
            )
            return {"H": 0.0, "Cp": 1.0, "K_values": {k: 1.0 for k in z}}

        H = gas.enthalpy_mole  # J/kmol → convert to J/mol
        Cp = gas.cp_mole       # J/(kmol·K) → J/(mol·K)

        # K-values via chemical equilibrium at constant T, P.
        K_values = self._compute_K_values(z, T, P, composition)

        return {"H": H / 1000.0, "Cp": Cp / 1000.0, "K_values": K_values}

    def compute_unit(
        self,
        unit_type: str,
        config: dict,
        inlet: dict,
    ) -> Optional[dict]:
        """Run reactor units via Cantera; return None for standard unit types."""
        if unit_type not in _REACTOR_UNIT_TYPES:
            return None
        return self._run_reactor(unit_type, config, inlet)

    def teardown(self) -> None:
        self._gas = None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _compute_K_values(
        self, z: dict, T: float, P: float, composition: str
    ) -> dict:
        """Estimate K-values from equilibrium vs. feed composition.

        Uses a simple ratio of equilibrium mole fractions to feed fractions
        as a proxy for vapour–liquid K-values.  For gas-phase kinetics the
        K-values are typically all 1.0 (single-phase), so this is adequate.
        """
        ct = self._ct
        try:
            gas_eq = ct.Solution(self._gas.source)
            gas_eq.TPX = T, P, composition
            gas_eq.equilibrate("TP")
            K_values: dict = {}
            for name, frac in z.items():
                try:
                    idx = gas_eq.species_index(name)
                    x_eq = gas_eq.X[idx]
                    x_feed = max(frac, 1e-30)
                    K_values[name] = x_eq / x_feed
                except Exception:
                    K_values[name] = 1.0
        except Exception as exc:
            logger.debug(f"CanteraProvider: equilibrium K-value calc failed ({exc}); using K=1.0")
            K_values = {k: 1.0 for k in z}
        return K_values

    def _run_reactor(self, unit_type: str, config: dict, inlet: dict) -> dict:
        """Run a Cantera IdealGasReactor for CSTR/PFR/IdealGasReactor units."""
        ct = self._ct
        gas = self._gas

        T = inlet.get("T", 300.0)
        P = inlet.get("P", 101325.0)
        z = inlet.get("z", {})
        flowrate = inlet.get("flowrate", 1.0)
        residence_time = config.get("residence_time", 1.0)

        composition = ", ".join(f"{k}:{v}" for k, v in z.items() if v > 0)
        try:
            gas.TPX = T, P, composition
        except Exception as exc:
            logger.warning(
                f"CanteraProvider: reactor inlet state error ({exc}). "
                "Returning inlet unchanged."
            )
            return dict(inlet)

        reactor = ct.IdealGasReactor(gas)
        net = ct.ReactorNet([reactor])

        try:
            net.advance(residence_time)
        except Exception as exc:
            logger.warning(
                f"CanteraProvider: reactor integration failed ({exc}). "
                "Returning inlet unchanged."
            )
            return dict(inlet)

        # Build outlet from post-reaction gas state.
        thermo = reactor.thermo
        out_z = {
            name: float(thermo.X[thermo.species_index(name)])
            for name in z
            if name in thermo.species_names
        }
        # Re-normalise if any species were missing.
        total = sum(out_z.values()) or 1.0
        out_z = {k: v / total for k, v in out_z.items()}

        logger.debug(
            f"CanteraProvider: {unit_type} T_out={thermo.T:.1f} K, "
            f"P_out={thermo.P:.0f} Pa"
        )
        return {
            "T": float(thermo.T),
            "P": float(thermo.P),
            "flowrate": flowrate,
            "z": out_z,
        }


# Self-register so the registry knows about "cantera" without a hard import.
register_provider("cantera", CanteraProvider)
