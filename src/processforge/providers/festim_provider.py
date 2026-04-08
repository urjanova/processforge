from __future__ import annotations

from typing import Optional

from loguru import logger

from .base import AbstractProvider
from .registry import register_provider


class FestimProvider(AbstractProvider):
    """Hydrogen transport provider backed by FESTIM.
    
    This provider acts as a bridge for material property lookup and provider registration
    for Festim membrane units. The dynamic finite element calculation is primarily 
    handled within the unit itself.
    """

    def __init__(self):
        self._materials = {}
        self._initialized = False

    def initialize(self, provider_config: dict, flowsheet_config: dict) -> None:
        """Initialize the FESTIM provider with material properties."""
        try:
            import festim  # noqa: F401
        except ImportError as exc:
            raise RuntimeError(
                "FESTIM is not installed. "
                "Install with: pip install 'festim'"
            ) from exc

        # Load material configurations if provided in the provider block
        self._materials = provider_config.get("materials", {})
        self._initialized = True
        logger.info(f"FestimProvider initialized with {len(self._materials)} materials")

    def get_material_props(self, material_id: str) -> dict:
        """Get the material properties required by FESTIM.
        
        If not explicitly in the provider config, attempts to provide some sane defaults or
        throws an error if strict checking is desired. Here we just return what's configured.
        """
        props = self._materials.get(material_id)
        if not props:
            # Provide some default properties (e.g. for steel/iron) if not specified
            # or raise an error. We'll return a default set for demonstration.
            logger.warning(f"Material '{material_id}' not explicitly configured in FestimProvider. Using defaults.")
            return {
                "D_0": 4.1e-7,
                "E_D": 0.39,
                "S_0": 4.28e-5,
                "E_S": 0.26
            }
        return props

    def get_thermo_properties(self, stream: dict) -> dict:
        """Fall back to CoolProp — Festim focuses on hydrogen transport."""
        from processforge.thermo import get_Cp_molar, get_enthalpy_molar, get_K_values

        z, T, P = stream["z"], stream["T"], stream["P"]
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
        """Compute outlet stream for FESTIM units.
        
        Since FestimMembrane handles its own dynamic evaluation and ODE stepping,
        this steady-state intercept can just return None to let the unit run itself.
        """
        return None

    def teardown(self) -> None:
        """Release any provider-level resources."""
        self._initialized = False

# Register the provider
register_provider("festim", FestimProvider)
