"""CoolPropProvider — default provider wrapping the existing thermo.py functions."""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from .base import AbstractProvider
from .errors import ProviderNotAvailableError

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
        except ImportError as exc:
            raise ProviderNotAvailableError(
                "CoolProp is not installed. To use the CoolProp thermodynamics provider, "
                "please install it by running `pip install \"processforge[coolprop]\"`"
            ) from exc

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
        from CoolProp.CoolProp import AbstractState

        for comp in z.keys():
            try:
                st = AbstractState("PR", comp)
                st.specify_phase(self.CP.iphase_liquid)
                st.update(self.CP.PT_INPUTS, P, T)
                phi_L = st.fugacity_coefficient(0)
                st.unspecify_phase()
                st.specify_phase(self.CP.iphase_gas)
                st.update(self.CP.PT_INPUTS, P, T)
                phi_V = st.fugacity_coefficient(0)
                st.unspecify_phase()
                Ks[comp] = phi_L / phi_V if phi_V != 0 else 1.0
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
