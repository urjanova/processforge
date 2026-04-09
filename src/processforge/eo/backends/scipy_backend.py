"""ScipyBackend: Newton-Raphson with finite-difference Jacobian (no extra deps)."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.sparse.linalg as spla
from loguru import logger

from .base import AbstractEOBackend

if TYPE_CHECKING:
    from ..jacobian import GlobalJacobianManager


class ScipyBackend(AbstractEOBackend):
    """
    Pure-scipy EO backend.

    Uses ``GlobalJacobianManager.evaluate_F`` and ``evaluate_J`` (finite
    differences) to drive a Newton-Raphson loop with Armijo backtracking.

    No additional dependencies beyond numpy and scipy.
    """

    def solve(
        self,
        manager: "GlobalJacobianManager",
        x0: np.ndarray,
        tol: float = 1e-6,
        max_iter: int = 50,
        log_every: int = 5,
    ) -> tuple[np.ndarray, bool, dict]:
        x = x0.copy()
        converged = False
        final_norm = float("inf")

        for iteration in range(max_iter):
            F = manager.evaluate_F(x)
            norm_F = float(np.max(np.abs(F)))
            if iteration % log_every == 0:
                logger.info(f"  NR iter {iteration:3d}: ||F||_inf = {norm_F:.3e}")
            else:
                logger.debug(f"NR iter {iteration}: ||F||_inf = {norm_F:.3e}")

            if norm_F < tol:
                converged = True
                final_norm = norm_F
                logger.info(
                    f"EO converged in {iteration} iterations "
                    f"(||F||_inf = {norm_F:.3e})"
                )
                break

            J = manager.evaluate_J(x)

            try:
                dx = spla.spsolve(J, -F)
            except Exception as exc:  # singular / ill-conditioned
                logger.warning(f"Linear solve failed at iter {iteration}: {exc}")
                break

            # Armijo backtracking line search
            alpha = 1.0
            x_new = x + alpha * dx
            F_new = manager.evaluate_F(x_new)
            norm_new = float(np.max(np.abs(F_new)))

            backtrack = 0
            while norm_new >= norm_F and backtrack < 20:
                alpha *= 0.5
                x_new = x + alpha * dx
                F_new = manager.evaluate_F(x_new)
                norm_new = float(np.max(np.abs(F_new)))
                backtrack += 1

            if alpha < 1e-12:
                logger.warning(
                    f"Line search stalled at iter {iteration} "
                    f"(alpha={alpha:.2e})"
                )
                break

            x = x_new
            final_norm = norm_new
        else:
            logger.warning(
                f"EO did not converge after {max_iter} iterations "
                f"(||F||_inf = {final_norm:.3e})"
            )

        stats = {"iterations": iteration + 1, "final_norm": final_norm}
        return x, converged, stats
