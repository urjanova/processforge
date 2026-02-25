"""CasADiBackend: solves the EO system using CasADi symbolic AD + nlpsol."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from loguru import logger

from .base import AbstractEOBackend

if TYPE_CHECKING:
    from ..jacobian import GlobalJacobianManager


class CasADiBackend(AbstractEOBackend):
    """
    Accelerator EO backend.  Builds a CasADi symbolic system by calling
    ``unit.get_casadi_residuals()`` on each unit, then assembles a root-finding
    problem solved with ``casadi.rootfinder`` (Newton / KINSOL).

    Requires: ``casadi>=3.6``.

    Units that raise ``NotImplementedError`` in ``get_casadi_residuals()``
    (e.g. Heater, Flash) are not supported — fall back to ScipyBackend.
    """

    def solve(
        self,
        manager: "GlobalJacobianManager",
        x0: np.ndarray,
        tol: float = 1e-6,
        max_iter: int = 50,
    ) -> tuple[np.ndarray, bool, dict]:
        try:
            import casadi as ca
        except ImportError as exc:
            raise ImportError(
                "CasADiBackend requires 'casadi'. "
                "Install with: pip install 'processforge[eo-casadi]'"
            ) from exc

        logger.info("Building CasADi symbolic system")
        comps = manager.components
        n = manager.n_vars

        # Declare symbolic variable vector
        x_sym = ca.SX.sym("x", n)

        # Populate CasADi-mode pointers in each StreamVar
        for name, sv in manager._streams.items():
            o = sv.global_offset
            sv.ca_T = x_sym[o]
            sv.ca_P = x_sym[o + 1]
            sv.ca_F = x_sym[o + 2]
            sv.ca_z = {c: x_sym[o + 3 + i] for i, c in enumerate(comps)}

        # Assemble residual vector
        F_sym_parts: list = []

        # Feed boundary conditions
        for fname, feed in manager._feed_values.items():
            sv = manager._streams[fname]
            o = sv.global_offset
            F_sym_parts.append(x_sym[o] - feed.get("T", 298.15))
            F_sym_parts.append(x_sym[o + 1] - feed.get("P", 101325.0))
            F_sym_parts.append(x_sym[o + 2] - feed.get("flowrate", 1.0))
            z_feed = feed.get("z", {})
            for i, c in enumerate(comps):
                F_sym_parts.append(x_sym[o + 3 + i] - z_feed.get(c, 0.0))

        # Unit residuals
        for entry in manager._unit_entries:
            unit = entry["unit"]
            config = entry["config"]
            inlet_names = entry["inlet_names"]
            outlet_names = entry["outlet_names"]

            inlet_syms = self._make_inlet_syms(
                manager, x_sym, inlet_names, comps
            )

            for outlet_name in outlet_names:
                sv_out = manager._streams[outlet_name]
                outlet_syms = {
                    "T": sv_out.ca_T,
                    "P": sv_out.ca_P,
                    "F": sv_out.ca_F,
                    "z": dict(sv_out.ca_z),
                }
                try:
                    residuals = unit.get_casadi_residuals(
                        inlet_syms, outlet_syms, config
                    )
                except NotImplementedError:
                    raise NotImplementedError(
                        f"Unit '{getattr(unit, 'name', unit)}' does not support "
                        "CasADiBackend. Use backend='scipy' for flowsheets "
                        "containing Heater or Flash units."
                    )
                F_sym_parts.extend(residuals)

        F_sym = ca.vertcat(*F_sym_parts)

        # Build root-finding function and solve
        rfn = ca.rootfinder("rfn", "newton", {"x": x_sym, "g": F_sym},
                            {"abstol": tol, "max_iter": max_iter})
        try:
            sol = rfn(x0)
            x_sol = np.array(sol).flatten()
            F_check = manager.evaluate_F(x_sol)
            converged = float(np.max(np.abs(F_check))) < tol * 10
            stats = {"final_norm": float(np.max(np.abs(F_check)))}
            if converged:
                logger.info("CasADi rootfinder converged")
            else:
                logger.warning(
                    f"CasADi rootfinder may not have converged "
                    f"(||F||_inf = {stats['final_norm']:.3e})"
                )
            return x_sol, converged, stats
        except Exception as exc:
            logger.error(f"CasADi rootfinder failed: {exc}")
            return x0.copy(), False, {"error": str(exc)}

    @staticmethod
    def _make_inlet_syms(
        manager: "GlobalJacobianManager",
        x_sym,
        inlet_names: list[str],
        comps: list[str],
    ) -> dict:
        """Build inlet symbol dict, merging multiple inlets symbolically."""
        import casadi as ca

        if len(inlet_names) == 1:
            sv = manager._streams[inlet_names[0]]
            return {
                "T": sv.ca_T, "P": sv.ca_P, "F": sv.ca_F,
                "z": dict(sv.ca_z),
            }

        # Symbolic flow-weighted merge
        F_total = sum(manager._streams[n].ca_F for n in inlet_names)
        T_merged = sum(
            manager._streams[n].ca_T * manager._streams[n].ca_F
            for n in inlet_names
        ) / (F_total + 1e-20)
        P_merged = sum(
            manager._streams[n].ca_P * manager._streams[n].ca_F
            for n in inlet_names
        ) / (F_total + 1e-20)
        z_merged = {
            c: sum(
                manager._streams[n].ca_z.get(c, ca.SX(0)) * manager._streams[n].ca_F
                for n in inlet_names
            ) / (F_total + 1e-20)
            for c in comps
        }
        return {"T": T_merged, "P": P_merged, "F": F_total, "z": z_merged}
