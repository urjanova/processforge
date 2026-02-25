"""ValveEOMixin: EO residuals and Pyomo/CasADi constraints for the Valve unit."""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ..mixin import EOUnitModelMixin

if TYPE_CHECKING:
    from ..stream_var import StreamVar


class ValveEOMixin(EOUnitModelMixin):
    """EO interface for Valve (explicit, all backends supported)."""

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

        ratio = self.pressure_ratio  # type: ignore[attr-defined]

        m.add_component(
            f"{self.name}_eq_P",  # type: ignore[attr-defined]
            pyo.Constraint(expr=outlet_var.pyo_P == ratio * inlet_var.pyo_P),
        )
        m.add_component(
            f"{self.name}_eq_T",  # type: ignore[attr-defined]
            pyo.Constraint(expr=outlet_var.pyo_T == inlet_var.pyo_T),
        )
        m.add_component(
            f"{self.name}_eq_F",  # type: ignore[attr-defined]
            pyo.Constraint(expr=outlet_var.pyo_F == inlet_var.pyo_F),
        )
        for c in inlet_var.components:
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
        ratio = self.pressure_ratio  # type: ignore[attr-defined]
        residuals = [
            outlet_syms["P"] - ratio * inlet_syms["P"],
            outlet_syms["T"] - inlet_syms["T"],
            outlet_syms["F"] - inlet_syms["F"],
        ]
        for c in outlet_syms["z"]:
            residuals.append(outlet_syms["z"][c] - inlet_syms["z"][c])
        return residuals
