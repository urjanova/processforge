import copy
import math


class Pipes:
    """
    Simple steady-state pipes model.
    Represents pipes connecting units, with a fixed pressure drop due to friction.
    Assumes isenthalpic (no temperature change).
    """

    def __init__(self, name, delta_p=1000.0, diameter=0.1):
        self.name = name
        self.delta_p = delta_p
        self.diameter = diameter

    def run(self, inlet):
        """
        Processes the inlet stream through a pipe unit, calculating outlet properties including pressure drop and flow rate.
        This method assumes laminar flow and uses a simplified Hagen-Poiseuille equation for flow calculation.
        The pressure drop is applied, but not below 1000.0 Pa to avoid unrealistic values.
        Temperature and other properties are passed through unchanged except for the unit identifier.
        
        Handles both single-value snapshots and timeseries (lists).
        
        Args:
            inlet (dict): A dictionary containing inlet stream properties, such as 'P' (pressure in Pa),
                          'T' (temperature), and potentially others. Defaults for 'P' is 101325.0 Pa if not provided.
        Returns:
            dict: A dictionary representing the outlet stream with updated 'P' (pressure), 'flow' (volumetric flow rate in m³/s),
                  'T' (temperature, unchanged), and 'unit' set to "Pipes". Other keys from inlet are deep-copied.
        Notes:
            - Assumes mu = 0.001 Pa·s (viscosity), L = 1.0 m (length).
            - Flow is set to 0 if diameter is not positive.
            - Uses copy.deepcopy to ensure outlet is independent of inlet.
        """
        outlet = copy.deepcopy(inlet)
        inlet_p = inlet.get("P", 101325.0)
        
        # Handle both scalar and list inputs
        if isinstance(inlet_p, list):
            outlet["P"] = [max(p - self.delta_p, 1000.0) for p in inlet_p]
            # Calculate flow assuming laminar flow, mu = 0.001 Pa.s, L = 1 m
            mu = 0.001
            L = 1.0
            delta_p_used = [inlet_p[i] - outlet["P"][i] for i in range(len(inlet_p))]
            outlet["flow"] = [
                math.pi * self.diameter**4 * dp / (128 * mu * L) if self.diameter > 0 else 0
                for dp in delta_p_used
            ]
        else:
            outlet["P"] = max(inlet_p - self.delta_p, 1000.0)
            mu = 0.001
            L = 1.0
            delta_p_used = inlet_p - outlet["P"]
            outlet["flow"] = (
                math.pi * self.diameter**4 * delta_p_used / (128 * mu * L)
                if self.diameter > 0 else 0
            )
        
        outlet["T"] = inlet["T"]
        outlet["unit"] = "Pipes"
        return outlet

    def run_dynamic(self, inlet_stream_ts, time_range, t_eval, solver):
        """
        Simulates dynamic operation of the pipes unit over a time step dt.
        For this simple model, assumes no accumulation and behaves like steady-state.

        Parameters:
        inlet (dict): Inlet stream dictionary.
        time_range (tuple): Time range for simulation (not used in this simple model).
        t_eval (array-like): Time points for evaluation (not used in this simple model).
        solver (object): Solver object (not used in this simple model).

        Returns:
        dict: Outlet stream dictionary, same as run method.
        """
        return self.run(inlet_stream_ts)
