import copy

class Strainer:
    """
    Simplified pressure-drop element.
    Removes a fixed ΔP from the inlet pressure.
    Adiabatic, isenthalpic.
    """

    def __init__(self, name, deltaP=5000.0):
        self.name = name
        self.deltaP = deltaP

    def run(self, inlet):
        outlet = copy.deepcopy(inlet)
        inlet_P = inlet.get("P", 101325.0)
        outlet["P"] = max(inlet_P - self.deltaP, 1000.0)
        outlet["T"] = inlet["T"]
        outlet["unit"] = "Strainer"
        return outlet
