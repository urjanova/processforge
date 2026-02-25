"""AbstractEOBackend: base class for all EO solver backends."""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..jacobian import GlobalJacobianManager


class AbstractEOBackend(ABC):
    """Abstract base for EO solver backends (Pyomo, CasADi, Scipy)."""

    @abstractmethod
    def solve(
        self,
        manager: "GlobalJacobianManager",
        x0: np.ndarray,
        tol: float = 1e-6,
        max_iter: int = 50,
    ) -> tuple[np.ndarray, bool, dict]:
        """
        Solve the global EO system F(x) = 0.

        Args:
            manager: Assembled ``GlobalJacobianManager``.
            x0: Initial guess vector.
            tol: Convergence tolerance on ||F||_inf.
            max_iter: Maximum Newton iterations.

        Returns:
            (x_solution, converged, stats) where ``stats`` is a dict with
            iteration count and final residual norm.
        """
