from .. import thermo

class Flash:
    """Represents an isothermal flash drum unit operation in a process flowsheet.

    This class models the separation of a multi-component feed stream into
    a vapor and a liquid phase. The separation is performed at a specified
    pressure and the temperature of the inlet stream.

    Attributes:
        name (str): The unique name of the flash drum unit.
        inlet (str): The name/identifier of the inlet stream.
        out_vap (str): The name/identifier for the outlet vapor stream.
        out_liq (str): The name/identifier for the outlet liquid stream.
        P (float): The operating pressure of the flash drum.
    """

    """Initializes the Flash unit with its configuration.

    Args:
        name (str): The unique name for the flash drum unit.
        params (dict): A dictionary containing the configuration parameters.
            Expected keys are:
            - "in" (str): The name of the inlet stream.
            - "out_vap" (str): The name of the outlet vapor stream.
            - "out_liq" (str): The name of the outlet liquid stream.
            - "P" (float): The operating pressure of the flash drum.
    """

    """Executes the flash calculation for the unit.

    This method performs an isothermal flash calculation using the Rachford-Rice
    equation to determine the vapor fraction (beta). It then calculates the
    compositions of the resulting liquid and vapor phases. The outlet streams
    are added to the main streams dictionary.

    Args:
        streams (dict): A dictionary of all process streams, where keys are
            stream names and values are dictionaries of stream properties
            (e.g., {"T": temp, "P": pressure, "z": composition_dict}).

    Returns:
        dict: The updated streams dictionary containing the newly calculated
                outlet liquid and vapor streams.
    """


    def __init__(self, name, params):
        self.name = name
        self.inlet = params["in"]
        self.out_vap = params["out_vap"]
        self.out_liq = params["out_liq"]
        self.P = params["P"]

    def run(self, streams):
        """Performs an isothermal flash calculation on an inlet stream.
        This method takes a dictionary of all process streams, identifies the
        inlet stream for this flash unit, and performs a vapor-liquid
        equilibrium calculation at the inlet temperature and the unit's
        specified pressure.
        The calculation involves:
        1.  Retrieving K-values (vapor-liquid equilibrium constants) for each
            component at the specified temperature and pressure.
        2.  Solving the Rachford-Rice equation to determine the vapor
            fraction (β) of the outlet mixture.
        3.  Calculating the compositions of the resulting liquid (x) and
            vapor (y) phases.
        4.  Creating two new outlet streams (liquid and vapor) and adding them
            to the streams dictionary.
        Args:
            streams (dict): A dictionary where keys are stream names and values
                are dictionaries containing stream properties like temperature ('T'),
                pressure ('P'), and composition ('z').
        Returns:
            dict: The updated streams dictionary containing the newly created
                liquid and vapor outlet streams from the flash calculation.
        """

        feed = streams[self.inlet]
        T = feed["T"]
        P = self.P
        z = feed["z"]  # dict of mole fractions

        # Get K-values for components
        Ks = thermo.get_K_values(z.keys(), T, P)

        # Solve for vapor fraction β
        beta = thermo.rachford_rice(z, Ks)

        # Compute liquid & vapor compositions
        x = {i: z[i] / (1 + beta*(Ks[i]-1)) for i in z}
        y = {i: Ks[i] * x[i] for i in z}

        # Normalize to sum=1
        x_tot = sum(x.values())
        y_tot = sum(y.values())
        x = {i: x[i]/x_tot for i in x}
        y = {i: y[i]/y_tot for i in y}

        streams[self.out_liq] = {"T": T, "P": P, "z": x, "phase": "liq", "beta": beta}
        streams[self.out_vap] = {"T": T, "P": P, "z": y, "phase": "vap", "beta": beta}

        return streams
