from __future__ import annotations

from typing import Any
from loguru import logger
import numpy as np

from .provider_mixin import ProviderMixin


class FestimMembrane(ProviderMixin):
    """FESTIM-backed membrane for hydrogen permeation processing.
    
    This unit uses the FESTIM solver under the hood to calculate dynamic hydrogen 
    concentration profiles across a material boundary. It then updates stream quantities
    by subtracting the "lost" permeated amount.

    Flowsheet JSON definition expects:
    - thickness [m]
    - mesh_cells (int)
    - area [m^2]
    - material_id
    """

    def __init__(self, name: str, **params):
        self.name = name
        self.params = params
        
        # Configuration
        self.thickness = params.get("thickness", 1e-3)
        self.mesh_cells = params.get("mesh_cells", 100)
        self.area = params.get("area", 1.0)
        self.material_id = params.get("material", "default")
        
        # State tracking
        self.current_time = 0.0
        self._festim_model: Any = None
        self._inner_pressure: Any = None
        self._temp_obj: Any = None

    def _setup_festim_model(self):
        """Sets up the mesh and materials once to avoid recreation overhead."""
        if self._festim_model is not None:
            return
            
        import festim as F
        
        self._festim_model = F.Simulation()
        self._festim_model.mesh = F.MeshFromVertices(
            np.linspace(0, self.thickness, self.mesh_cells)
        )
        
        # Retrieve material properties from provider
        # self.material is set by flowsheet._build_units() after construction;
        # fall back to self.material_id which is parsed from params in __init__.
        mat_id = getattr(self, "material", self.material_id)
        mat_props = {}
        if hasattr(self, "_provider") and hasattr(self._provider, "get_material_props"):
            mat_props = self._provider.get_material_props(mat_id)
        else:
            # Sane default values if provider is missing or doesn't support get_material_props
            mat_props = {"D_0": 4.1e-7, "E_D": 0.39, "S_0": 4.28e-5, "E_S": 0.26}
            
        self._festim_model.materials = F.Material(
            id=1,
            D_0=mat_props["D_0"],
            E_D=mat_props["E_D"]
        )
        
        # We will dynamically update these boundary values during evaluate()
        self._inner_pressure = F.as_constant(0.0)

        # Inner Wall: Sieverts' law driven by the process fluid
        # Outer Wall: Recombination or Dirichlet (assuming zero concentration to environment)
        self._festim_model.boundary_conditions = [
            F.SievertsBC(
                surfaces=1,
                S_0=mat_props["S_0"],
                E_S=mat_props["E_S"],
                pressure=self._inner_pressure
            ),
            F.DirichletBC(surfaces=2, value=0, field=0)
        ]

        self._temp_obj = F.Temperature(300.0)
        self._festim_model.T = self._temp_obj
        
        # Transient settings but manually stepped
        self._festim_model.settings = F.Settings(
            absolute_tolerance=1e-10, relative_tolerance=1e-10,
            transient=True, final_time=1e9 # controlled by dt later
        )
        self._festim_model.dt = F.Stepsize(initial_value=1.0)
        
        self._festim_model.initialise()
        logger.info(f"FestimMembrane '{self.name}': initialized underlying FESTIM model.")

    def _run_impl(self, inlet: dict) -> dict:
        """Run single steady-state timestep natively."""
        import festim as F
        
        # 1. Extract current conditions
        T = inlet.get("T", 300.0)
        P_total = inlet.get("P", 101325.0)
        flowrate = inlet.get("flowrate", 0.0)
        z = inlet.get("z", {})
        
        y_H2 = z.get("H2", 0.0)
        if y_H2 <= 0.0:
            # No H2 to permeate
            return {
                "retentate_out": dict(inlet),
                "permeate_out": {"T": T, "P": P_total, "flowrate": 0.0, "z": {}}
            }

        # 2. Update model
        self._setup_festim_model()
        self._temp_obj.value = T  # applied when initialise() calls create_functions below
        self._inner_pressure.assign(P_total * y_H2)

        # For steady state, we can simulate a large transient or switch transient=False
        # For simplicity, we just force a steady simulation configuration transiently
        self._festim_model.settings.transient = False
        self._festim_model.initialise() # Re-init to pickup the false transient setting
        self._festim_model.run()
        
        # 3. Get molar loss rate
        flux_1d = F.SurfaceFlux(field="solute", surface=2).compute(self._festim_model.h_transport_problem.mobile.solution)
        # flux is in [molecules/m2/s] or [mol/m2/s] depending on Festim config. We assume mol/m2/s per our constants.
        loss_rate = flux_1d * self.area
        
        # Clamp to avoid going negative flow
        h2_flow_in = flowrate * y_H2
        actual_loss = min(loss_rate, h2_flow_in)
        
        # Calculate retentate state
        retentate = dict(inlet)
        new_flowrate = flowrate - actual_loss
        new_z = dict(z)
        if new_flowrate > 0:
            new_z["H2"] = (h2_flow_in - actual_loss) / new_flowrate
            # renormalize other species
            for comp in new_z:
                if comp != "H2":
                    new_z[comp] = (flowrate * z[comp]) / new_flowrate
        retentate["flowrate"] = new_flowrate
        retentate["z"] = new_z
        
        # Permeate state
        permeate = {"T": T, "P": P_total, "flowrate": actual_loss, "z": {"H2": 1.0}}
        
        return {"retentate_out": retentate, "permeate_out": permeate}

    def run_dynamic(self, inlet_stream_ts, t_span, t_eval, solver):
        """Perform a dynamic run over the full time span."""
        import festim as F
        self._setup_festim_model()
        
        # We need to manually step FESTIM for each timestep in t_eval
        num_steps = len(t_eval)
        
        # Pre-allocate output time-series arrays
        res_retentate = {
            "time": t_eval.tolist(),
            "T": [0.0] * num_steps,
            "P": [0.0] * num_steps,
            "flowrate": [0.0] * num_steps,
            "z": {c: [0.0] * num_steps for c in inlet_stream_ts["z"].keys()}
        }
        res_permeate = {
            "time": t_eval.tolist(),
            "T": [0.0] * num_steps,
            "P": [0.0] * num_steps,
            "flowrate": [0.0] * num_steps,
            "z": {"H2": [1.0] * num_steps}
        }
        
        # Initial conditions from the first element of inlet_stream_ts
        self.current_time = t_span[0]
        
        # Enable transient simulation
        self._festim_model.settings.transient = True
        self._festim_model.initialise()
        
        logger.info(f"FestimMembrane '{self.name}': integrating ODE dynamically from {t_span[0]} to {t_span[1]}")

        for i in range(num_steps):
            # Extract current timestep condition
            T = inlet_stream_ts["T"][i]
            P = inlet_stream_ts["P"][i]
            flowrate = inlet_stream_ts["flowrate"][i]
            y_H2 = inlet_stream_ts["z"].get("H2", [0.0]*num_steps)[i]
            
            # Step size
            if i < num_steps - 1:
                dt = t_eval[i+1] - t_eval[i]
            else:
                dt = 1.0 # arbitrary small step for last point
                
            # Update temperature directly on the FEniCS Function (initialise already ran)
            from fenics import Constant, interpolate
            V = self._temp_obj.T.function_space()
            self._temp_obj.T.assign(interpolate(Constant(T), V))
            self._temp_obj.T_n.assign(self._temp_obj.T)
            self._inner_pressure.assign(P * y_H2)
            
            self._festim_model.dt.value.assign(dt)
            self._festim_model.settings.final_time = self.current_time + dt
            
            # Run FESTIM one dt step
            if y_H2 > 0:
                self._festim_model.run()
                
            self.current_time += dt
            
            # Molar loss
            try:
                flux_1d = F.SurfaceFlux(field="solute", surface=2).compute(self._festim_model.h_transport_problem.mobile.solution)
                loss_rate = abs(flux_1d) * self.area  # Ensure positive magnitude
            except Exception:
                loss_rate = 0.0
                
            h2_flow_in = flowrate * y_H2
            actual_loss = min(loss_rate, h2_flow_in)
            
            new_flowrate = flowrate - actual_loss
            
            res_retentate["T"][i] = T
            res_retentate["P"][i] = P
            res_retentate["flowrate"][i] = new_flowrate
            
            if new_flowrate > 0:
                for comp, vals in inlet_stream_ts["z"].items():
                    if comp == "H2":
                        res_retentate["z"][comp][i] = (h2_flow_in - actual_loss) / new_flowrate
                    else:
                        res_retentate["z"][comp][i] = (flowrate * vals[i]) / new_flowrate
            
            res_permeate["T"][i] = T
            res_permeate["P"][i] = P
            res_permeate["flowrate"][i] = actual_loss
            res_permeate["z"]["H2"][i] = 1.0

        return {"retentate_out": res_retentate, "permeate_out": res_permeate}
