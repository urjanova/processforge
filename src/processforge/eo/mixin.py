"""EOUnitModelMixin: abstract interface for unit models to expose EO residuals."""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .stream_var import StreamVar


class EOUnitModelMixin(ABC):
    """
    Mixin that steady-state unit classes inherit to expose equation-oriented
    residual equations alongside their existing sequential-modular ``run()`` method.

    A unit with one inlet and one outlet exposes ``n_vars_out = N_c + 3`` residuals
    (one per outlet variable: T, P, F, z[0..N_c-1]).

    For **explicit units** (Pump, Valve, Strainer, Pipes) whose SM ``run()`` already
    computes the correct outlet analytically, the default ``get_scipy_residuals()``
    delegates to ``run()`` to avoid duplicating equation logic::

        residuals[j] = outlet_vals[j] - run(inlet_vals)[j]   for j in 0..n_vars-1

    For **implicit/thermodynamic units** (Heater, Flash) that call CoolProp
    internally, subclasses must override ``get_scipy_residuals()`` directly with the
    raw equation residuals.  These units raise ``NotImplementedError`` for
    ``build_pyomo_block()`` and ``get_casadi_residuals()``.
    """

    # Subclasses may override to declare multi-outlet units (e.g. Flash = 2)
    eo_n_outlet_streams: int = 1

    @abstractmethod
    def build_pyomo_block(
        self,
        m: Any,
        inlet_var: "StreamVar",
        outlet_var: "StreamVar",
        config: dict,
    ) -> None:
        """
        Add this unit's equations as Pyomo Constraint objects to model ``m``.

        Constraint naming convention::

            m.add_component(f"{self.name}_eq_T",  pyo.Constraint(...))
            m.add_component(f"{self.name}_eq_P",  pyo.Constraint(...))
            m.add_component(f"{self.name}_eq_F",  pyo.Constraint(...))
            m.add_component(f"{self.name}_eq_z_{comp}", pyo.Constraint(...))

        Args:
            m: Pyomo ``ConcreteModel`` where stream variables already exist.
            inlet_var: ``StreamVar`` for the inlet stream.
            outlet_var: ``StreamVar`` for the outlet stream.
            config: Unit config dict from the flowsheet JSON.
        """

    @abstractmethod
    def get_casadi_residuals(
        self,
        inlet_syms: dict,
        outlet_syms: dict,
        config: dict,
    ) -> list:
        """
        Return a list of CasADi SX residual expressions (one per outlet variable).

        Args:
            inlet_syms:  ``{"T": ca.SX, "P": ca.SX, "F": ca.SX, "z": {comp: ca.SX}}``
            outlet_syms: same structure for outlet variables.
            config: Unit config dict.

        Returns:
            list of ``ca.SX`` expressions of length ``outlet_var.n_vars``.
        """

    def get_scipy_residuals(
        self,
        inlet_vals: dict,
        outlet_vals: dict,
        config: dict,
    ) -> list[float]:
        """
        Evaluate residuals as plain Python floats (for ScipyBackend / FD Jacobian).

        Default implementation delegates to the unit's SM ``run()`` method, making
        explicit units zero-copy: the EO residual is ``outlet - run(inlet)``.

        Subclasses with implicit thermodynamic equations (Heater, Flash) must
        override this to return raw equation residuals without calling ``run()``.

        Args:
            inlet_vals:  stream dict snapshot  ``{"T": float, "P": float,
                         "flowrate": float, "z": {comp: float}}``.
            outlet_vals: same structure for the current outlet variable values.
            config: Unit config dict.

        Returns:
            list of floats of length ``N_c + 3``.
        """
        expected = self.run(inlet_vals)  # type: ignore[attr-defined]
        comps = getattr(self, "components", list(inlet_vals.get("z", {}).keys()))
        residuals = [
            outlet_vals["T"] - expected["T"],
            outlet_vals["P"] - expected["P"],
            outlet_vals["flowrate"] - expected["flowrate"],
        ]
        for c in comps:
            residuals.append(
                outlet_vals["z"].get(c, 0.0) - expected["z"].get(c, 0.0)
            )
        return residuals
