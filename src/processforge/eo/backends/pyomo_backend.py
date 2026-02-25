"""PyomoBackend: solves the EO system via Pyomo ConcreteModel + IPOPT."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from loguru import logger

from .base import AbstractEOBackend

if TYPE_CHECKING:
    from ..jacobian import GlobalJacobianManager


class PyomoBackend(AbstractEOBackend):
    """
    Primary EO backend.  Builds a Pyomo ``ConcreteModel``, populates it by
    calling ``unit.build_pyomo_block()`` on each registered unit, then solves
    with IPOPT via ``SolverFactory("ipopt")``.

    Requires: ``pyomo`` (and IPOPT binary on PATH).

    Units that raise ``NotImplementedError`` in ``build_pyomo_block()``
    (e.g. Heater, Flash) are not supported — fall back to ScipyBackend.
    """

    def solve(
        self,
        manager: "GlobalJacobianManager",
        x0: np.ndarray,
        tol: float = 1e-6,
        max_iter: int = 300,
    ) -> tuple[np.ndarray, bool, dict]:
        try:
            import pyomo.environ as pyo
        except ImportError as exc:
            raise ImportError(
                "PyomoBackend requires 'pyomo'. "
                "Install with: pip install 'processforge[eo]'"
            ) from exc

        logger.info("Building Pyomo model")
        m = pyo.ConcreteModel()

        # --- declare stream variables ---
        comps = manager.components
        stream_names = list(manager._streams.keys())

        m.streams = pyo.Set(initialize=stream_names)
        m.components = pyo.Set(initialize=comps)

        m.T = pyo.Var(m.streams, domain=pyo.PositiveReals, bounds=(200.0, 1000.0))
        m.P = pyo.Var(m.streams, domain=pyo.PositiveReals, bounds=(100.0, 1e8))
        m.F = pyo.Var(m.streams, domain=pyo.NonNegativeReals, bounds=(0.0, 1e6))
        m.z = pyo.Var(m.streams, m.components, bounds=(0.0, 1.0))

        # Populate Pyomo-mode pointers in each StreamVar
        for name, sv in manager._streams.items():
            sv.pyo_T = m.T[name]
            sv.pyo_P = m.P[name]
            sv.pyo_F = m.F[name]
            sv.pyo_z = {c: m.z[name, c] for c in comps}

        # Initialise variables from x0
        for name, sv in manager._streams.items():
            vals = sv.to_vector_slice(x0)
            m.T[name].set_value(vals["T"])
            m.P[name].set_value(vals["P"])
            m.F[name].set_value(vals["flowrate"])
            for c in comps:
                m.z[name, c].set_value(vals["z"].get(c, 0.0))

        # --- feed boundary conditions ---
        feed_constraints = {}
        for fname, feed in manager._feed_values.items():
            feed_constraints[f"bc_T_{fname}"] = pyo.Constraint(
                expr=m.T[fname] == feed.get("T", 298.15)
            )
            feed_constraints[f"bc_P_{fname}"] = pyo.Constraint(
                expr=m.P[fname] == feed.get("P", 101325.0)
            )
            feed_constraints[f"bc_F_{fname}"] = pyo.Constraint(
                expr=m.F[fname] == feed.get("flowrate", 1.0)
            )
            z_feed = feed.get("z", {})
            for c in comps:
                feed_constraints[f"bc_z_{fname}_{c}"] = pyo.Constraint(
                    expr=m.z[fname, c] == z_feed.get(c, 0.0)
                )
        for cname, con in feed_constraints.items():
            m.add_component(cname, con)

        # --- unit model constraints ---
        for entry in manager._unit_entries:
            unit = entry["unit"]
            config = entry["config"]
            inlet_names = entry["inlet_names"]
            outlet_names = entry["outlet_names"]

            if len(inlet_names) == 1:
                inlet_sv = manager._streams[inlet_names[0]]
            else:
                # Multi-inlet: build a merged inlet StreamVar (virtual)
                # For Pyomo we just use the first outlet as the merge point;
                # full multi-inlet Pyomo support is a future enhancement
                inlet_sv = manager._streams[inlet_names[0]]
                logger.warning(
                    "PyomoBackend: multi-inlet merge not yet fully symbolic; "
                    "using first inlet only."
                )

            for outlet_name in outlet_names:
                outlet_sv = manager._streams[outlet_name]
                try:
                    unit.build_pyomo_block(m, inlet_sv, outlet_sv, config)
                except NotImplementedError:
                    raise NotImplementedError(
                        f"Unit '{getattr(unit, 'name', unit)}' does not support "
                        "PyomoBackend. Use backend='scipy' for flowsheets "
                        "containing Heater or Flash units."
                    )

        # --- solve ---
        solver = pyo.SolverFactory("ipopt")
        solver.options["tol"] = tol
        solver.options["max_iter"] = max_iter
        logger.info("Calling IPOPT")
        result = solver.solve(m, tee=False)

        ok = str(result.solver.termination_condition) == "optimal"
        if not ok:
            logger.warning(
                f"IPOPT termination: {result.solver.termination_condition}"
            )

        # Extract solution into x vector
        x_sol = x0.copy()
        for name, sv in manager._streams.items():
            vals = {
                "T": pyo.value(m.T[name]),
                "P": pyo.value(m.P[name]),
                "flowrate": pyo.value(m.F[name]),
                "z": {c: pyo.value(m.z[name, c]) for c in comps},
            }
            sv.write_to_vector(x_sol, vals)

        stats = {
            "termination": str(result.solver.termination_condition),
            "iterations": getattr(
                result.solver, "iterations", -1
            ),
        }
        return x_sol, ok, stats
