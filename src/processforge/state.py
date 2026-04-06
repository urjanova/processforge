import os
import json
import zarr
import numpy as np
from loguru import logger
from .eo.flowsheet import EOFlowsheet

class StateManager:
    """Manages the .pfstate Zarr store and handles drift detection."""

    def __init__(self, state_path: str):
        self.state_path = state_path

    def save_state(self, config: dict, x_converged: np.ndarray, var_names: list[str]):
        """Save the configured flowsheet and converged variables to .pfstate."""
        store = zarr.storage.LocalStore(self.state_path)
        root = zarr.group(store=store, overwrite=True)
        
        # Save config as JSON string attribute
        root.attrs["config"] = json.dumps(config)
        
        # Save variables
        root.create_dataset("x", data=x_converged)
        root.attrs["var_names"] = var_names
        logger.info(f"Saved .pfstate to {self.state_path}")

    def load_state(self) -> dict | None:
        """Load the state from .pfstate if it exists."""
        if not os.path.exists(self.state_path):
            return None
            
        try:
            store = zarr.storage.LocalStore(self.state_path)
            root = zarr.group(store=store, overwrite=False)
            config = json.loads(root.attrs["config"])
            x_converged = root["x"][:]
            var_names = root.attrs["var_names"]
            return {"config": config, "x": x_converged, "var_names": var_names}
        except Exception as e:
            logger.warning(f"Failed to load state from {self.state_path}: {e}")
            return None

    def detect_drift(self, current_config: dict, state: dict) -> list[str]:
        """Compare current flowsheet config (desired) with state config (actual).
        Returns a list of drifted parameter paths.
        """
        old_config = state["config"]
        drifted = []

        # Check stream parameters
        for s_name, s_data in current_config.get("streams", {}).items():
            old_s = old_config.get("streams", {}).get(s_name, {})
            for param in ["T", "P", "flowrate"]:
                if s_data.get(param) != old_s.get(param):
                    drifted.append(f"streams.{s_name}.{param}")
            
            # Check compositions
            cur_z = s_data.get("z", {})
            old_z = old_s.get("z", {})
            for comp in cur_z:
                if cur_z.get(comp) != old_z.get(comp):
                    drifted.append(f"streams.{s_name}.z.{comp}")

        # Check unit parameters
        for u_name, u_cfg in current_config.get("units", {}).items():
            old_u = old_config.get("units", {}).get(u_name, {})
            for param, val in u_cfg.get("parameters", {}).items():
                if val != old_u.get("parameters", {}).get(param):
                    drifted.append(f"units.{u_name}.parameters.{param}")

        return drifted

