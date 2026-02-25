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
