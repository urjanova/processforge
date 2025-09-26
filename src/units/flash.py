from .. import thermo

class Flash:


    def __init__(self, name, params):
        self.name = name
        self.inlet = params["in"]
        self.out_vap = params["out_vap"]
        self.out_liq = params["out_liq"]
        self.P = params["P"]

    def run(self, streams):
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
