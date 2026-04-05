import copy

from ..eo.units.strainer_eo import StrainerEOMixin
from .provider_mixin import ProviderMixin


class Strainer(ProviderMixin, StrainerEOMixin):
    """
    Simplified pressure-drop element.
    Removes a fixed ΔP from the inlet pressure.
    Adiabatic, isenthalpic.
    """

    def __init__(self, name, deltaP=5000.0):
        self.name = name
        self.deltaP = deltaP

    def _run_impl(self, inlet):
        """
        Run the strainer unit process.

        This method simulates the operation of a strainer unit by adjusting the pressure
        of the inlet stream and returning the modified outlet stream. The pressure is
        reduced by the deltaP value, but not below 1000.0 Pa. The temperature remains
        unchanged, and the unit identifier is set to "Strainer".

        Args:
            inlet (dict): A dictionary representing the inlet stream properties, expected
                to contain keys like "P" (pressure in Pa, defaults to 101325.0 if missing)
                and "T" (temperature).

        Returns:
            dict: A dictionary representing the outlet stream properties, with updated
                "P", unchanged "T", and added "unit" key.
        """

        outlet = copy.deepcopy(inlet)
        inlet_P = inlet.get("P", 101325.0)
        outlet["P"] = max(inlet_P - self.deltaP, 1000.0)
        outlet["T"] = inlet["T"]
        outlet["unit"] = "Strainer"
        return outlet
