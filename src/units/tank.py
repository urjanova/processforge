import numpy as np
from scipy.integrate import solve_ivp
from ..thermo import get_enthalpy_molar, get_Cp_molar
from dataclasses import dataclass


@dataclass
class TankState:
    T: float
    P: float
    z: dict
    n_total: float


class Tank:
    """
    Simple well-mixed molar tank.
    Parameters expected in params:
      - vol: (optional) placeholder, not used explicitly here
      - inlet: name of inlet stream (must exist in flowsheet streams)
      - outlet_flow: molar flowrate out [mol/s] (constant)
      - initial_n: {comp: moles} initial contents (dict)
      - initial_T: initial temperature [K]
      - P: pressure in Pa (assumed constant)
      - duty: heat duty [W] applied to tank (positive heats)
      - flow_in: molar inlet flowrate [mol/s] (used if inlet stream not providing flowrate) 
    """

    def __init__(self, name, **params):
        self.name = name
        self.params = params
        self.inlet = params.get("inlet")
        self.outflow = params.get("outlet_flow", 1.0)
        self.P = params.get("P", 101325)
        self.duty = params.get("duty", 0.0)
        # initial composition: dict of component -> moles
        self.initial_n = params.get("initial_n", {})
        self.initial_T = params.get("initial_T", 300.0)
        self.flow_in = params.get("flow_in", 1.0)  # mol/s if inlet doesn't specify

    def run(self, inlet):
        """
        Steady-state tank operation.
        For steady-state, the tank acts as a well-mixed vessel at equilibrium.
        The outlet composition equals the tank composition (assumed same as inlet for steady-state),
        outlet flow rate is fixed by the tank parameter, and temperature/pressure are from the tank.
        
        Args:
            inlet (dict): Inlet stream with T, P, z, flowrate
            
        Returns:
            dict: Outlet stream with tank's outlet flow rate and mixed conditions
        """
        import copy
        outlet = copy.deepcopy(inlet)
        
        # Tank outlet conditions
        outlet["flowrate"] = self.outflow
        outlet["P"] = self.P
        
        # For steady-state with no duty, temperature stays approximately the same
        # With duty, we'd need energy balance - simplified here
        if self.duty != 0:
            # Simplified: assume Cp ~ 75 J/mol/K for water-like mixture
            Cp_approx = 75.0
            F_in = inlet.get("flowrate", self.flow_in)
            if F_in > 1e-10:
                dT = self.duty / (F_in * Cp_approx)
                outlet["T"] = inlet.get("T", self.initial_T) + dT
            else:
                outlet["T"] = inlet.get("T", self.initial_T)
        
        outlet["unit"] = "Tank"
        return outlet

    def _ode_rhs(self, t, y, inlet_stream):
        """
        Computes the right-hand side of the ordinary differential equations (ODEs) for the tank unit.
        This method calculates the time derivatives of the molar amounts of species and the temperature
        based on inlet and outlet flows, enthalpy balances, and heat duty. It is used in numerical
        integration of the tank's dynamic behavior.
        Parameters
        ----------
        t : float
            Current time (in seconds). Not used in calculations but required for ODE solver interface.
        y : array-like
            State vector containing molar amounts of species (n) followed by temperature (T).
            Shape: (len(comps) + 1,), where comps are the components from inlet_stream.
        inlet_stream : dict
            Dictionary containing inlet stream properties:
            - "z": dict of mole fractions for each component.
            - "flowrate": inlet molar flowrate (mol/s), defaults to self.flow_in if not provided.
            - "T": inlet temperature (K).
            - "P": inlet pressure (Pa), defaults to self.P if not provided.
        Returns
        -------
        np.ndarray
            Array of time derivatives: [dn_dt[0], dn_dt[1], ..., dn_dt[N], dT_dt],
            where N is the number of components. Shape: (len(comps) + 1,).
        Notes
        -----
        - Molar balances are computed for each species based on inlet and outlet flows.
        - Energy balance accounts for enthalpy differences and heat duty (self.duty).
        - Cp_mix is the molar heat capacity of the mixture; a small positive value is enforced if zero or negative.
        - Assumes constant outlet flowrate (self.outflow) and tank pressure (self.P).
        """
        
        n = y[:-1]
        T = y[-1]
        comps = list(inlet_stream["z"].keys())
        n_tot = max(1e-12, sum(n))
        x = {comps[i]: max(0.0, n[i] / n_tot) for i in range(len(comps))}

        # inlet properties
        Fin = inlet_stream.get("flowrate", self.flow_in)  # mol/s
        z_in = inlet_stream["z"]
        T_in = inlet_stream["T"]
        P_in = inlet_stream.get("P", self.P)

        # outlet flow (constant)
        Fout = self.outflow

        # species molar balances
        dn_dt = np.zeros(len(comps))
        for i, comp in enumerate(comps):
            dn_dt[i] = Fin * z_in.get(comp, 0.0) - Fout * x.get(comp, 0.0)

        # enthalpies
        H_in = get_enthalpy_molar(z_in, T_in, P_in)  # J/mol
        H_out = get_enthalpy_molar(x, T, self.P)  # J/mol

        # Cp mix
        Cp_mix = get_Cp_molar(x, T, self.P)  # J/mol/K
        if Cp_mix <= 0:
            Cp_mix = 1e-6

        # energy balance (molar basis)
        # dT/dt = (Fin*H_in - Fout*H_out + Q) / (n_tot * Cp_mix)
        Q = self.duty
        dT_dt = (Fin * H_in - Fout * H_out + Q) / (n_tot * Cp_mix)

        return np.concatenate([dn_dt, [dT_dt]])

    def run_dynamic(self, inlet_stream_ts, t_span, t_eval, solver):
        """
        inlet_stream_ts: dict of timeseries data for the inlet stream
        t_span: (t0, tf)
        t_eval: list or array of times
        solver: solver object (currently unused but good for future)
        Returns: dict of timeseries data for the outlet stream
        """
        # For a dynamic tank, the inlet conditions are assumed constant for now.
        # We take the initial value of the inlet timeseries.
        # A more complex model could use a time-varying inlet.
        inlet_snapshot = {}
        for prop, values in inlet_stream_ts.items():
            if prop == "time":
                continue
            if prop == "z":
                # 'z' is a dict of lists, so we build a dict of the first values
                inlet_snapshot['z'] = {comp: ts[0] for comp, ts in values.items()}
            else:
                # Other properties are lists of values
                inlet_snapshot[prop] = values[0]

        comps = list(inlet_snapshot["z"].keys())

        # initial state vector
        n0 = [self.initial_n.get(c, 0.0) for c in comps]
        y0 = np.array(n0 + [self.initial_T])

        # integrate
        sol = solve_ivp(
            self._ode_rhs,
            t_span,
            y0,
            t_eval=t_eval,
            args=(inlet_snapshot,),
            rtol=1e-6,
            atol=1e-8,
        )

        # Reformat results into a timeseries dictionary for the outlet stream
        outlet_ts = {
            "time": sol.t.tolist(),
            "T": sol.y[-1, :].tolist(),
            "P": [self.P] * len(sol.t),
            "flowrate": [self.outflow] * len(sol.t),
        }
        # Add composition timeseries
        z_ts = {c: [] for c in comps}
        for i in range(len(sol.t)):
            n_vec = sol.y[:-1, i]
            n_tot = max(1e-12, sum(n_vec))
            for j, comp in enumerate(comps):
                z_ts[comp].append(n_vec[j] / n_tot)
        outlet_ts["z"] = z_ts

        return outlet_ts
