"""EOSolver: selects and runs the appropriate EO backend."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from loguru import logger

from .backends import ScipyBackend, PyomoBackend, CasADiBackend
from .backends.base import AbstractEOBackend

if TYPE_CHECKING:
    from .jacobian import GlobalJacobianManager

_BACKENDS: dict[str, type[AbstractEOBackend]] = {
    "scipy": ScipyBackend,
    "pyomo": PyomoBackend,
    "casadi": CasADiBackend,
}


class EOSolver:
    """
    Selects and runs the EO backend specified by ``backend`` name.

    Args:
        backend: One of ``"scipy"``, ``"pyomo"``, ``"casadi"``.
        tol: Convergence tolerance on ``||F||_inf``.
        max_iter: Maximum Newton iterations (Scipy) or IPOPT iterations (Pyomo).
    """

    def __init__(
        self,
        backend: str = "pyomo",
        tol: float = 1e-6,
        max_iter: int = 50,
    ) -> None:
        if backend not in _BACKENDS:
            raise ValueError(
                f"Unknown backend '{backend}'. Choose from: {list(_BACKENDS)}"
            )
        self.backend_name = backend
        self.tol = tol
        self.max_iter = max_iter
        self._backend: AbstractEOBackend = _BACKENDS[backend]()

    def solve(
        self,
        manager: "GlobalJacobianManager",
        x0: np.ndarray,
    ) -> tuple[np.ndarray, bool, dict]:
        """
        Solve the global system F(x) = 0.

        Returns:
            (x_solution, converged, stats)
        """
        logger.info(
            f"EOSolver: solving {manager.n_vars}-variable system "
            f"via '{self.backend_name}' backend"
        )
        x_sol, converged, stats = self._backend.solve(
            manager, x0, tol=self.tol, max_iter=self.max_iter
        )
        if not converged:
            logger.warning(
                f"EOSolver: system did not converge "
                f"(final ||F||_inf = {stats.get('final_norm', '?'):.3e})"
            )
        return x_sol, converged, stats

def solve_with_homotopy(
    fs: "processforge.eo.flowsheet.EOFlowsheet",
    manager: "GlobalJacobianManager",
    solver: "EOSolver",
    state: dict,
    drifted: list[str],
) -> tuple["np.ndarray", bool, dict]:
    """Try standard solve; if fails, fall back to step-wise homotopy."""
    old_config = state["config"]
    current_config = fs.config

    # First, let's just attempt a regular solve using state's x as x0 
    x0 = np.array(state["x"])
    x_sol, converged, stats = solver.solve(manager, x0)
    
    if converged:
        logger.info("Standard solver converged with previous state's warm guess.")
        return x_sol, converged, stats
        
    logger.warning("Standard solve failed. Invoking Homotopy 'step-wise apply' solver...")
    
    # We incrementally step the parameters
    steps = 10
    from copy import deepcopy
    
    # Create interpolations for drifted numerical parameters
    def update_config_value(cfg, path, value):
        parts = path.split('.')
        d = cfg
        for p in parts[:-1]:
            d = d[p]
        d[parts[-1]] = value
        
    def get_config_value(cfg, path):
        parts = path.split('.')
        d = cfg
        for p in parts[:-1]:
            if p not in d:
                return None
            d = d[p]
        return d.get(parts[-1], None)
        
    drifts = []
    for d in drifted:
        old_val = get_config_value(old_config, d)
        new_val = get_config_value(current_config, d)
        if isinstance(old_val, (int, float)) and isinstance(new_val, (int, float)):
            drifts.append((d, old_val, new_val))
            
    if not drifts:
        logger.error("No continuous numerical parameters to interpolate. Homotopy fails.")
        return x_sol, False, stats

    x_current = np.array(state["x"])
    interpolated_config = deepcopy(old_config)
    
    for step in range(1, steps + 1):
        alpha = step / steps
        logger.info(f"--- Homotopy Step {step}/{steps} (alpha={alpha:.1f}) ---")
        
        # Interpolate config
        for d, old_val, new_val in drifts:
            val = old_val + alpha * (new_val - old_val)
            update_config_value(interpolated_config, d, val)
        
        # We need the solver to use the new parameter values, but how?
        # The variables manager encapsulates the flowsheet config during _build().
        # So we have to rebuild the manager and unit objects with the interpolated config.
        # This requires recreating the flowsheet model
        
        from .flowsheet import EOFlowsheet
        
        step_fs = EOFlowsheet(interpolated_config, backend=fs.backend)
        from processforge.providers.manager import teardown_providers
        step_manager = step_fs._build()
        try:
            step_x0 = x_current.copy()
            step_fs.x0 = step_x0
            step_fs.var_names = step_fs._build_var_names(step_manager)
            step_solver = EOSolver(backend=fs.backend, tol=solver.tol, max_iter=solver.max_iter)
            
            x_current, conv, s_stats = step_solver.solve(step_manager, step_x0)
            
            if not conv:
                logger.error(f"Homotopy failed at step {step} (alpha={alpha:.1f})")
                return x_current, False, s_stats
        finally:
            teardown_providers(step_fs._provider_map)
            
    logger.info("Homotopy step-wise solve successfully reached target config.")
    return x_current, True, {"method": "homotopy", "steps": steps}

def solve_with_homotopy(
    fs, manager, solver, state, drifted
) -> "tuple[np.ndarray, bool, dict]":
    """Try standard solve; if fails, fall back to step-wise homotopy."""
    old_config = state["config"]
    current_config = fs.config

    # Standard solve first
    from loguru import logger
    import numpy as np
    
    x0 = np.array(state["x"])
    x_sol, converged, stats = solver.solve(manager, x0)
    if converged:
        logger.info("Standard solver converged with previous state's warm guess.")
        return x_sol, converged, stats
        
    logger.warning("Standard solve failed. Invoking Homotopy 'step-wise apply' solver...")
    
    def update_config_value(cfg, path, value):
        from copy import deepcopy
        parts = path.split('.')
        d = cfg
        for p in parts[:-1]:
            d = d[p]
        d[parts[-1]] = float(value)
        
    def get_config_value(cfg, path):
        parts = path.split('.')
        d = cfg
        for p in parts[:-1]:
            if p not in d:
                return None
            d = d[p]
        return d.get(parts[-1], None)

    drifts = []
    for d in drifted:
        old_val = get_config_value(old_config, d)
        new_val = get_config_value(current_config, d)
        if isinstance(old_val, (int, float)) and isinstance(new_val, (int, float)):
            drifts.append((d, old_val, new_val))
            
    if not drifts:
        logger.error("No continuous numerical parameters to interpolate. Homotopy fails.")
        return x_sol, False, stats

    x_current = np.array(state["x"])
    
    from copy import deepcopy
    interpolated_config = deepcopy(old_config)
    steps = 10
    
    from .flowsheet import EOFlowsheet
    from processforge.providers.manager import teardown_providers
    
    for step in range(1, steps + 1):
        alpha = step / steps
        logger.info(f"--- Homotopy Step {step}/{steps} (alpha={alpha:.1f}) ---")
        
        for d, old_val, new_val in drifts:
            val = old_val + alpha * (new_val - old_val)
            update_config_value(interpolated_config, d, val)
        
        step_fs = EOFlowsheet(interpolated_config, backend=fs.backend)
        try:
            step_manager = step_fs._build()
            step_x0 = x_current.copy()
            step_fs.x0 = step_x0
            step_fs.var_names = step_fs._build_var_names(step_manager)
            
            step_solver = EOSolver(backend=fs.backend, tol=solver.tol, max_iter=solver.max_iter)
            x_current, step_converged, s_stats = step_solver.solve(step_manager, step_x0)
            
            if not step_converged:
                logger.error(f"Homotopy failed to converge at step {step} (alpha={alpha:.1f})")
                return x_current, False, s_stats
        finally:
            teardown_providers(step_fs._provider_map)
            
    logger.info("Homotopy step-wise solve successfully reached target config.")
    # Re-run the final step exactly on current manager to populate the outer variables identically.
    x_sol, conv, stats = solver.solve(manager, x_current)
    return x_sol, conv, stats
