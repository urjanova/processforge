"""CoolPropProvider — default provider wrapping the existing thermo.py functions."""
from __future__ import annotations

import math
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

    ``get_thermo_properties`` mixes pure-component CoolProp values with
    mole-fraction weighting (log-linear for viscosity) and additionally
    returns ``molar_volume`` [m³/mol], ``density`` [kg/m³], ``molar_mass``
    [kg/mol] and ``viscosity`` [Pa·s] when CoolProp can supply them for the
    stream's components.
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
        molar_mass = 0.0
        molar_volume = 0.0
        ln_visc = 0.0
        mass_ok = vol_ok = visc_ok = True
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

            # Optional derived/transport properties (best-effort per component)
            try:
                molar_mass += frac * self.CP.PropsSI("MOLARMASS", "T", T, "P", P, comp)
            except Exception:  # noqa: BLE001
                mass_ok = False
            try:
                dmolar = self.CP.PropsSI("DMOLAR", "T", T, "P", P, comp)
                if dmolar > 0:
                    molar_volume += frac / dmolar
                else:
                    vol_ok = False
            except Exception:  # noqa: BLE001
                vol_ok = False
            try:
                visc = self.CP.PropsSI("VISCOSITY", "T", T, "P", P, comp)
                if visc > 0:
                    ln_visc += frac * math.log(visc)
                else:
                    visc_ok = False
            except Exception:  # noqa: BLE001
                visc_ok = False

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

        props = {
            "H": H,
            "Cp": Cp,
            "K_values": Ks,
        }
        if mass_ok and molar_mass > 0:
            props["molar_mass"] = molar_mass
            if vol_ok and molar_volume > 0:
                props["molar_volume"] = molar_volume
                props["density"] = molar_mass / molar_volume
        if visc_ok and ln_visc != 0.0:
            props["viscosity"] = math.exp(ln_visc)
        return props

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
