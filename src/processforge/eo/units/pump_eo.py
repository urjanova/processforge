"""PumpEOMixin: EO residuals and Pyomo/CasADi constraints for the Pump unit."""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ..mixin import EOUnitModelMixin

if TYPE_CHECKING:
    from ..stream_var import StreamVar

# Water properties used in pump model (matches SM pump.py assumptions)
_RHO = 1000.0    # kg/m³
_CP = 4180.0     # J/kg·K
_MW = 0.018      # kg/mol


class PumpEOMixin(EOUnitModelMixin):
    """
    EO interface for Pump.

    Explicit unit — ``get_scipy_residuals`` delegates to SM ``run()``.
    Pyomo and CasADi backends implement the same adiabatic pump equations
    symbolically.
    """

    def build_pyomo_block(
        self,
        m: Any,
        inlet_var: "StreamVar",
        outlet_var: "StreamVar",
        config: dict,
    ) -> None:
        try:
            import pyomo.environ as pyo
        except ImportError as exc:
            raise ImportError("PyomoBackend requires 'pyomo'.") from exc

        delta_p = self.deltaP  # type: ignore[attr-defined]
        eta = max(self.efficiency, 1e-6)  # type: ignore[attr-defined]
        dt_pump = delta_p / (_RHO * _CP * eta)  # K per (Pa / (rho*Cp*eta))

        m.add_component(
            f"{self.name}_eq_P",  # type: ignore[attr-defined]
            pyo.Constraint(expr=outlet_var.pyo_P == inlet_var.pyo_P + delta_p),
        )
        m.add_component(
            f"{self.name}_eq_T",  # type: ignore[attr-defined]
            pyo.Constraint(expr=outlet_var.pyo_T == inlet_var.pyo_T + dt_pump),
        )
        m.add_component(
            f"{self.name}_eq_F",  # type: ignore[attr-defined]
            pyo.Constraint(expr=outlet_var.pyo_F == inlet_var.pyo_F),
        )
        comps = inlet_var.components
        for c in comps:
            m.add_component(
                f"{self.name}_eq_z_{c}",  # type: ignore[attr-defined]
                pyo.Constraint(
                    expr=outlet_var.pyo_z[c] == inlet_var.pyo_z[c]
                ),
            )

    def get_casadi_residuals(
        self,
        inlet_syms: dict,
        outlet_syms: dict,
        config: dict,
    ) -> list:
        delta_p = self.deltaP  # type: ignore[attr-defined]
        eta = max(self.efficiency, 1e-6)  # type: ignore[attr-defined]
        dt_pump = delta_p / (_RHO * _CP * eta)

        residuals = [
            outlet_syms["P"] - (inlet_syms["P"] + delta_p),
            outlet_syms["T"] - (inlet_syms["T"] + dt_pump),
            outlet_syms["F"] - inlet_syms["F"],
        ]
        for c in outlet_syms["z"]:
            residuals.append(outlet_syms["z"][c] - inlet_syms["z"][c])
        return residuals
