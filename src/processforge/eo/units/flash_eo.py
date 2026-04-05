"""FlashEOMixin: EO residuals for the Flash unit (Scipy backend only, 2 outlets)."""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ..mixin import EOUnitModelMixin

if TYPE_CHECKING:
    from ..stream_var import StreamVar

_NOT_IMPL_MSG = (
    "FlashEOMixin only supports ScipyBackend. Flash uses CoolProp K-values "
    "which are not symbolically differentiable. "
    "Use backend='scipy' for flowsheets containing Flash units."
)


class FlashEOMixin(EOUnitModelMixin):
    """
    EO interface for Flash (Scipy backend only, 2 outlet streams).

    The Flash unit exposes residuals for two outlets (liquid + vapor):
        Liquid: T_liq = T_in, P_liq = P_flash, F_liq = F_in*(1-β), z_liq = x[c]
        Vapor:  T_vap = T_in, P_vap = P_flash, F_vap = F_in*β,     z_vap = y[c]

    where β, x, y are computed internally by solving the Rachford-Rice equation
    with CoolProp K-values (same as the SM ``run()`` method).
    """

    eo_n_outlet_streams: int = 2

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

    def get_scipy_residuals_multi_outlet(
        self,
        inlet_vals: dict,
        outlet_vals_list: list[dict],
        config: dict,
    ) -> list[float]:
        """
        Residuals for both liquid and vapor outlets.

        Args:
            inlet_vals: Feed stream dict.
            outlet_vals_list: [liq_vals, vap_vals] — current EO variable values.
            config: Unit config (must contain ``out_liq`` and ``out_vap``).

        Returns:
            Flat list: ``[liq_residuals..., vap_residuals...]``.
        """
        from processforge import thermo

        t_in = inlet_vals["T"]
        p_flash = self.P  # type: ignore[attr-defined]
        f_in = inlet_vals["flowrate"]
        z_in = inlet_vals["z"]
        comps = list(z_in.keys())

        # Route K-values through provider when attached; fall back to CoolProp.
        provider = getattr(self, "_provider", None)
        if provider is not None:
            props = provider.get_thermo_properties(
                {"z": z_in, "T": t_in, "P": p_flash}
            )
            k_vals = props.get("K_values", thermo.get_K_values(comps, t_in, p_flash))
        else:
            k_vals = thermo.get_K_values(comps, t_in, p_flash)
        beta = thermo.rachford_rice(z_in, k_vals)

        x = {c: z_in[c] / (1 + beta * (k_vals[c] - 1)) for c in z_in}
        y = {c: k_vals[c] * x[c] for c in z_in}
        x_tot = sum(x.values()) or 1.0
        y_tot = sum(y.values()) or 1.0
        x = {c: x[c] / x_tot for c in x}
        y = {c: y[c] / y_tot for c in y}

        f_liq = f_in * (1.0 - beta)
        f_vap = f_in * beta

        liq_vals, vap_vals = outlet_vals_list[0], outlet_vals_list[1]

        residuals: list[float] = []

        # Liquid outlet residuals
        residuals += [
            liq_vals["T"] - t_in,
            liq_vals["P"] - p_flash,
            liq_vals["flowrate"] - f_liq,
        ]
        for c in comps:
            residuals.append(liq_vals["z"].get(c, 0.0) - x.get(c, 0.0))

        # Vapor outlet residuals
        residuals += [
            vap_vals["T"] - t_in,
            vap_vals["P"] - p_flash,
            vap_vals["flowrate"] - f_vap,
        ]
        for c in comps:
            residuals.append(vap_vals["z"].get(c, 0.0) - y.get(c, 0.0))

        return residuals

    def get_scipy_residuals(
        self,
        inlet_vals: dict,
        outlet_vals: dict,
        config: dict,
    ) -> list[float]:
        # Should not be called for Flash — use get_scipy_residuals_multi_outlet
        raise RuntimeError(
            "Flash has 2 outlets. Call get_scipy_residuals_multi_outlet() instead."
        )
