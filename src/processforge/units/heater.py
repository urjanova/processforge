import numpy as np
from scipy.optimize import brentq
from ..thermo import get_enthalpy_molar

class Heater:
    """
    Represents a simple heater unit in a process simulation.
    The Heater class models a unit that adds a specified amount of heat (duty)
    to a process stream, raising its temperature. It assumes no pressure drop
    and no change in composition.
    Attributes:
        name (str): The unique name of the heater unit.
        inlet (str): The key for the inlet stream in the streams dictionary.
        outlet (str): The key for the outlet stream in the streams dictionary.
        duty (float): The heat duty to be added to the stream, in Watts (W).
        flowrate (float): The molar flow rate of the stream, in mol/s.
    """
    """
    Initializes a Heater instance.
    Args:
        name (str): The name of the heater unit.
        params (dict): A dictionary of parameters for the heater.
            Expected keys:
            - "in" (str): The name of the inlet stream.
            - "out" (str): The name of the outlet stream.
            - "duty" (float, optional): Heat duty in Watts. Defaults to 0.0.
            - "flowrate" (float, optional): Molar flow rate in mol/s.
                                                Defaults to 1.0.
    """
    """
    Executes the heater simulation for a single time step.
    Calculates the outlet stream temperature by performing an energy balance.
    It determines the required outlet enthalpy based on the inlet enthalpy and
    the specified heat duty, then solves for the temperature that corresponds
    to this outlet enthalpy at constant pressure.
    Args:
        streams (dict): A dictionary containing all process streams, where each
                        stream is a dictionary with "T" (temperature),
                        "P" (pressure), and "z" (composition).
    Returns:
        dict: The updated streams dictionary, including the newly calculated
                outlet stream for this heater.
    """

    def __init__(self, name, params):
        self.name = name
        self.inlet = params["in"]
        self.outlet = params["out"]
        self.duty = params.get("duty", 0.0)  # W
        self.flowrate = params.get("flowrate", 1.0)  # mol/s (default 1)

    def run(self, streams):
        """
        Simulates the heater unit operation by calculating the outlet temperature based on enthalpy balance.
        This method takes the inlet stream properties, computes the enthalpy change due to the specified duty,
        and solves for the outlet temperature using a root-finding algorithm. If the solver fails, it defaults
        to the inlet temperature.
        Args:
            streams (dict): A dictionary containing stream data, where keys are stream names and values are
                    dictionaries with keys 'T' (temperature in K), 'P' (pressure in Pa), and 'z' (mole fractions).
        Returns:
            dict: The updated streams dictionary with the outlet stream added or modified, containing 'T', 'P', and 'z'.
        Raises:
            None explicitly, but relies on external functions like get_enthalpy_molar and brentq.
        Notes:
            - Assumes self.inlet, self.outlet, self.duty, and self.flowrate are set in the class instance.
            - Temperature bounds for solving are fixed between 100 K and 2000 K.
            - If the enthalpy solver fails, outlet temperature is set to inlet temperature as a fallback.
        """

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
