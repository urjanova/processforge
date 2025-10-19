import copy

class Valve:
    """
    Simple steady-state valve model.
    Reduces pressure according to a fixed ratio (P_out = ratio * P_in).
    Assumes isenthalpic (no temperature change).
    """

    def __init__(self, name, pressure_ratio=0.5):
        self.name = name
        self.pressure_ratio = pressure_ratio

    def run(self, inlet):
        outlet = copy.deepcopy(inlet)
        inlet_P = inlet.get("P", 101325.0)
        outlet["P"] = max(inlet_P * self.pressure_ratio, 1000.0)
        outlet["T"] = inlet["T"]
        outlet["unit"] = "Valve"
        return outlet
