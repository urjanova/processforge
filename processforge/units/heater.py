import numpy as np
from scipy.optimize import brentq
from ..thermo import get_enthalpy_molar

class Heater:
    def __init__(self, name, params):
        self.name = name
        self.inlet = params["in"]
        self.outlet = params["out"]
        self.duty = params.get("duty", 0.0)  # W
        self.flowrate = params.get("flowrate", 1.0)  # mol/s (default 1)

    def run(self, streams):
        feed = streams[self.inlet]
        T_in, P, z = feed["T"], feed["P"], feed["z"]

        # Enthalpy in
        H_in = get_enthalpy_molar(z, T_in, P)

        # Target enthalpy out
        H_target = H_in + self.duty / self.flowrate

        # Define residual function
        def residual(T):
            return get_enthalpy_molar(z, T, P) - H_target

        # Solve for outlet T (reasonable bounds: -100 K to 2000 K)
        try:
            T_out = brentq(residual, 100.0, 2000.0)
        except ValueError:
            # fallback if solver fails
            T_out = T_in

        streams[self.outlet] = {"T": T_out, "P": P, "z": z}
        return streams
