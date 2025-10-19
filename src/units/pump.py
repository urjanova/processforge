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
