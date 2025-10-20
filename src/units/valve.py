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
        """
        Simulates the operation of a valve unit by adjusting the pressure of the inlet stream.
        This method takes an inlet stream dictionary, applies a pressure reduction based on the
        valve's pressure ratio, ensures the outlet pressure does not drop below 1000.0 Pa,
        and returns the modified outlet stream with the same temperature and updated unit label.
        Parameters:
        inlet (dict): A dictionary representing the inlet stream, expected to contain keys like
                      'P' (pressure in Pa, defaults to 101325.0 if missing) and 'T' (temperature).
        Returns:
        dict: A deep copy of the inlet dictionary with updated 'P' (pressure), unchanged 'T' (temperature),
              and 'unit' set to "Valve".
        """

        outlet = copy.deepcopy(inlet)
        inlet_P = inlet.get("P", 101325.0)
        outlet["P"] = max(inlet_P * self.pressure_ratio, 1000.0)
        outlet["T"] = inlet["T"]
        outlet["unit"] = "Valve"
        return outlet
