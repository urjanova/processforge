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
        if manager.n_vars == 0:
            logger.info("EOSolver: empty system (no stream variables) — trivially converged.")
            return x0, True, {"iterations": 0, "final_norm": 0.0}
        x_sol, converged, stats = self._backend.solve(
            manager, x0, tol=self.tol, max_iter=self.max_iter
        )
        if not converged:
            logger.warning(
                f"EOSolver: system did not converge "
                f"(final ||F||_inf = {stats.get('final_norm', '?'):.3e})"
            )
        return x_sol, converged, stats

def compute_residual_breakdown(
    manager: "GlobalJacobianManager",
    x: "np.ndarray",
    var_names: list,
    top_n: int = 10,
) -> list:
    """Return top-N largest residuals with per-variable labels.

    Args:
        manager: The assembled Jacobian manager (must still be live).
        x: Current solution vector.
        var_names: Human-readable labels (``stream/field`` format).
        top_n: How many violators to return.

    Returns:
        List of dicts sorted descending by residual magnitude.
    """
    F = manager.evaluate_F(x)
    entries = []
    for i, label in enumerate(var_names):
        if i >= len(F):
            break
        parts = label.split("/", 1)
        entries.append({
            "index": i,
            "var_name": label,
            "stream": parts[0] if len(parts) == 2 else "",
            "field": parts[1] if len(parts) == 2 else label,
            "residual": float(abs(F[i])),
        })
    entries.sort(key=lambda e: e["residual"], reverse=True)
    return entries[:top_n]


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
