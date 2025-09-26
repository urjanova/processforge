import numpy as np
from scipy.integrate import solve_ivp
from ..thermo import get_enthalpy_molar, get_Cp_molar

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
    def __init__(self, name, params):
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

    def _ode_rhs(self, t, y, inlet_stream):
        """
        y: [n1, n2, ..., T]
        inlet_stream: dict with keys T, P, z (mole fractions) and flowrate optional
        """
        n = y[:-1]
        T = y[-1]
        comps = list(inlet_stream["z"].keys())
        n_tot = max(1e-12, sum(n))
        x = {comps[i]: max(0.0, n[i]/n_tot) for i in range(len(comps))}

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

    def run_dynamic(self, streams, t_span, t_eval):
        """
        streams: dict of current streams (must include inlet stream)
        t_span: (t0, tf)
        t_eval: list or array of times
        Returns: dict time -> streams snapshot (each snapshot is dict of stream->data)
        """
        if self.inlet not in streams:
            raise ValueError(f"Inlet stream {self.inlet} not found in flowsheet streams")

        inlet = streams[self.inlet]
        comps = list(inlet["z"].keys())
        # initial state vector
        n0 = [self.initial_n.get(c, 0.0) for c in comps]
        y0 = np.array(n0 + [self.initial_T])

        # integrate
        sol = solve_ivp(self._ode_rhs, t_span, y0, t_eval=t_eval, args=(inlet,), rtol=1e-6, atol=1e-8)

        results = {}
        for idx, t in enumerate(sol.t):
            y = sol.y[:, idx]
            n = y[:-1]
            T = y[-1]
            n_tot = max(1e-12, sum(n))
            x = {comps[i]: max(0.0, n[i]/n_tot) for i in range(len(comps))}
            # assemble tank stream state
            results[t] = {
                self.name: {
                    "T": float(T),
                    "P": self.P,
                    "z": x,
                    "n_total": float(n_tot)
                }
            }
        return results
