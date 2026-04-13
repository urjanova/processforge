"""CoolPropProvider — default provider wrapping the existing thermo.py functions."""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from .base import AbstractProvider

if TYPE_CHECKING:
    from processforge.types import CoolPropProviderConfig, FlowsheetConfig


class CoolPropProvider(AbstractProvider):
    """Default thermodynamic provider backed by CoolProp.

    Wraps the existing ``processforge.thermo`` functions with no behaviour
    change.  All existing flowsheets that omit a ``providers`` block continue
    to work exactly as before.
    """

    def initialize(
        self,
        provider_config: "CoolPropProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        try:
            import CoolProp.CoolProp as CP
            self.CP = CP
        except ImportError:
            raise ImportError(
                "CoolProp is not installed. To use the CoolProp thermodynamics provider, "
                "please install it by running `pip install \"processforge[coolprop]\"`"
            )

    def get_thermo_properties(self, stream: dict) -> dict:
        """Calculate thermodynamic properties using CoolProp."""
        from loguru import logger

        z = stream["z"]
        T = stream["T"]
        P = stream["P"]

        H = 0.0
        Cp = 0.0
        for comp, frac in z.items():
            if frac <= 0.0:
                continue
            
            try:
                H += frac * self.CP.PropsSI("HMOLAR", "T", T, "P", P, comp)
            except ValueError:
                logger.warning(f"Component '{comp}' not found in CoolProp. Enthalpy contribution is 0.")

            try:
                Cp += frac * self.CP.PropsSI("Cpmolar", "T", T, "P", P, comp)
            except ValueError:
                logger.warning(f"Component '{comp}' not found in CoolProp. Cp contribution is 0.")

        Ks = {}
        for comp in z.keys():
            try:
                fugL = self.CP.PropsSI("fugL", "T", T, "P", P, comp)
                fugV = self.CP.PropsSI("fugV", "T", T, "P", P, comp)
                Ks[comp] = fugL / fugV if fugV != 0 else 1.0
            except Exception:
                logger.warning(f"Could not calculate K-value for '{comp}'. Using fallback K=1.0.")
                Ks[comp] = 1.0

        return {
            "H": H,
            "Cp": Cp,
            "K_values": Ks,
        }

    def compute_unit(
        self,
        unit_type: str,
        config: dict,
        inlet: dict,
    ) -> Optional[dict]:
        # CoolProp provider never intercepts unit computation —
        # always fall through to the default SM unit logic.
        return None

    def teardown(self) -> None:
        pass
