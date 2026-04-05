"""HeaterEOMixin: EO residuals for the Heater unit (Scipy backend only)."""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ..mixin import EOUnitModelMixin

if TYPE_CHECKING:
    from ..stream_var import StreamVar

_NOT_IMPL_MSG = (
    "HeaterEOMixin only supports ScipyBackend. Heater uses CoolProp for enthalpy "
    "calculations which are not symbolically differentiable. "
    "Use backend='scipy' for flowsheets containing Heater units."
)


class HeaterEOMixin(EOUnitModelMixin):
    """
    EO interface for Heater (Scipy backend only).

    Implements the energy balance residual directly:
        H(T_out, P, z) - H(T_in, P, z) - duty/F = 0

    Pyomo and CasADi backends are not supported because CoolProp
    enthalpy calls are not symbolically differentiable.
    """

    def build_pyomo_block(
        self,
        m: Any,
        inlet_var: "StreamVar",
        outlet_var: "StreamVar",
        config: dict,
    ) -> None:
        raise NotImplementedError(_NOT_IMPL_MSG)

    def get_casadi_residuals(
        self,
        inlet_syms: dict,
        outlet_syms: dict,
        config: dict,
    ) -> list:
        raise NotImplementedError(_NOT_IMPL_MSG)

    def get_scipy_residuals(
        self,
        inlet_vals: dict,
        outlet_vals: dict,
        config: dict,
    ) -> list[float]:
        """
        Energy balance residuals for the heater.

        Equations:
            H(T_out, P, z) - H(T_in, P, z) - duty/F = 0   (energy balance)
            P_out - P_in = 0                               (isobaric)
            F_out - F_in = 0                               (mass balance)
            z_out[c] - z_in[c] = 0  for each component    (iso-composition)
        """
        from processforge.thermo import get_enthalpy_molar

        t_in = inlet_vals["T"]
        p_in = inlet_vals["P"]
        f_in = inlet_vals["flowrate"]
        z_in = inlet_vals["z"]

        t_out = outlet_vals["T"]
        p_out = outlet_vals["P"]
        f_out = outlet_vals["flowrate"]
        z_out = outlet_vals["z"]

        duty = self.duty  # type: ignore[attr-defined]
        comps = list(z_in.keys())

        # Route through provider when attached; fall back to CoolProp.
        provider = getattr(self, "_provider", None)
        if provider is not None:
            h_in = provider.get_thermo_properties(
                {"z": z_in, "T": t_in, "P": p_in}
            )["H"]
            h_out = provider.get_thermo_properties(
                {"z": z_out, "T": t_out, "P": p_out}
            )["H"]
        else:
            h_in = get_enthalpy_molar(z_in, t_in, p_in)
            h_out = get_enthalpy_molar(z_out, t_out, p_out)

        f_eff = f_in if abs(f_in) > 1e-20 else 1e-20
        residuals = [
            h_out - h_in - duty / f_eff,
            p_out - p_in,
            f_out - f_in,
        ]
        for c in comps:
            residuals.append(z_out.get(c, 0.0) - z_in.get(c, 0.0))
        return residuals
