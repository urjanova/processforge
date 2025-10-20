import copy

class Pump:
    """
    Simple steady-state pump model.
    Adds a fixed pressure rise ΔP to the inlet stream.
    Optionally computes shaft power and outlet temperature rise (adiabatic efficiency).
    """

    def __init__(self, name, deltaP=1e5, efficiency=0.8):
        self.name = name
        self.deltaP = deltaP
        self.efficiency = efficiency

    def run(self, inlet):
        """
        Simulates the operation of a pump unit by processing the inlet stream.
        This method takes an inlet stream dictionary, calculates the outlet conditions
        based on the pump's deltaP and efficiency, and returns the modified outlet stream.
        Assumptions include water properties at 298 K for density and specific heat.
        Args:
            inlet (dict): A dictionary representing the inlet stream with keys such as
                'P' (pressure in Pa, default 101325.0), 'flowrate' (molar flow rate in mol/s, default 1.0),
                and 'T' (temperature in K).
        Returns:
            dict: The outlet stream dictionary with updated 'P', 'T', 'power', and 'unit' keys.
                - 'P': Outlet pressure (inlet P + deltaP).
                - 'T': Outlet temperature (adjusted for adiabatic inefficiency).
                - 'power': Shaft power in J/s.
                - 'unit': Set to "Pump".
        """

        outlet = copy.deepcopy(inlet)
        inlet_P = inlet.get("P", 101325.0)
        outlet["P"] = inlet_P + self.deltaP

        # Density assumption (water at 298 K)
        rho = 1000.0  # kg/m³
        R = 8.314
        flow_mol = inlet.get("flowrate", 1.0)
        MW = 0.018  # kg/mol
        mass_flow = flow_mol * MW  # kg/s

        # Shaft power (J/s)
        power = (mass_flow / rho) * self.deltaP / max(self.efficiency, 1e-6)

        # Approximate temperature rise (adiabatic inefficiency)
        Cp = 4180.0  # J/kg-K for water
        dT = (power * (1 - self.efficiency)) / (mass_flow * Cp)
        outlet["T"] = inlet["T"] + dT

        outlet["power"] = power
        outlet["unit"] = "Pump"
        return outlet
