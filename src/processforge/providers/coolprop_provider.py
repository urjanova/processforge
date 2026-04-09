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
        # CoolProp needs no warm-up; nothing to do.
        pass

    def get_thermo_properties(self, stream: dict) -> dict:
        """Delegate to the existing CoolProp-backed thermo functions."""
        from processforge.thermo import get_enthalpy_molar, get_Cp_molar, get_K_values

        z = stream["z"]
        T = stream["T"]
        P = stream["P"]
        return {
            "H": get_enthalpy_molar(z, T, P),
            "Cp": get_Cp_molar(z, T, P),
            "K_values": get_K_values(list(z.keys()), T, P),
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
