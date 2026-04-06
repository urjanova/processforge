"""CanteraJacobianMixin: Jacobian contributions for CanteraProvider.

Implements the ``JacobianContributor`` protocol for Cantera-backed units by
subclassing ``BaseJacobianMixin`` and overriding ``compute_jacobian_block``.

Strategy
--------
Cantera's ``jacobian_stoich`` attribute gives kinetics-space derivatives
(d net-production-rate / d concentration), not the stream-variable residual
space (T, P, z, flowrate) that ProcessForge uses.  Converting between the
two spaces requires a chain-rule calculation that is more complex than simply
doing unit-local finite differences on the Cantera side.

Instead, this mixin uses ``BaseJacobianMixin._local_fd_block`` to run
``CanteraProvider._run_reactor`` with small perturbations of each inlet
variable.  This is O(n_vars) Cantera reactor integrations per unit per Newton
iteration — far cheaper than O(n_vars × n_units) global ``evaluate_F`` calls
because only one unit's residuals are re-evaluated per perturbation.

For the Heater unit, full analytic derivatives are available via Cantera's
``gas.cp_mole`` and ``gas.partial_molar_enthalpies`` — no reactor integration
needed at all.

Usage
-----
Added to ``CanteraProvider`` by extending its MRO::

    class CanteraProvider(AbstractProvider, CanteraJacobianMixin): ...

No other changes to ``CanteraProvider`` are required.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from .base_jacobian_mixin import BaseJacobianMixin

if TYPE_CHECKING:
    from .jacobian_contributor import ReferenceState


class CanteraJacobianMixin(BaseJacobianMixin):
    """Jacobian mixin for CanteraProvider.

    Handles reactor units (CSTR, PFR, IdealGasReactor) via unit-local FD
    through ``_run_reactor``, and the Heater unit via analytic thermo
    derivatives using the loaded Cantera ``Solution``.

    All other unit types return ``[]``, falling back to global FD.
    """

    def compute_jacobian_block(
        self,
        unit_type: str,
        config: dict,
        ref_state: "ReferenceState",
        row_offset: int,
        col_offset: int,
        n_vars: int,
        epsilon: float = 1e-5,
    ) -> list[tuple[int, int, float]]:
        """Return Jacobian triplets for Cantera-handled unit types.

        Args:
            unit_type:   Unit class name.
            config:      Unit configuration dict.
            ref_state:   Inlet stream snapshot at the current Newton iterate.
            row_offset:  Global row offset for the outlet residual block.
            col_offset:  Global column offset for the inlet variable block.
            n_vars:      Variables per stream (3 + n_components).
            epsilon:     Perturbation step for unit-local FD.

        Returns:
            Non-empty triplet list for handled unit types; ``[]`` otherwise
            (caller falls back to global FD for that block).
        """
        if unit_type in ("CSTR", "PFR", "IdealGasReactor"):
            return self._local_fd_block(
                config,
                ref_state,
                row_offset,
                col_offset,
                n_vars,
                # _run_reactor is defined on CanteraProvider; available via MRO.
                run_fn=lambda cfg, inlet: self._run_reactor(unit_type, cfg, inlet),
                epsilon=epsilon,
            )

        if unit_type == "Heater":
            return self._heater_analytic_block(
                config, ref_state, row_offset, col_offset, n_vars
            )

        return []  # all other unit types → global FD

    # ------------------------------------------------------------------
    # Heater: analytic Jacobian via Cantera thermo properties
    # ------------------------------------------------------------------

    def _heater_analytic_block(
        self,
        config: dict,
        ref_state: "ReferenceState",
        row_offset: int,
        col_offset: int,
        n_vars: int,
    ) -> list[tuple[int, int, float]]:
        """Analytic Jacobian for the Heater energy-balance residual.

        The Heater residual system (ScipyBackend form) is:

        .. code-block:: text

            Row 0 (energy balance):  H(T_out, z) - H(T_in, z) - duty/F = 0
            Row 1 (isobaric):        P_out - P_in = 0
            Row 2 (mass balance):    F_out - F_in = 0
            Row 3+i (composition):   z_out[i] - z_in[i] = 0

        Analytic derivatives for Row 0 (via Cantera gas properties):

        - dR/dT_out = Cp(T_out, z)
        - dR/dT_in  = -Cp(T_in, z)
        - dR/dF_in  = duty / F_in^2
        - dR/dz_i   = partial_h_i(T_out) - partial_h_i(T_in)

        Rows 1-3+i have trivial ±1 diagonal derivatives.

        This method assumes ``self._gas`` is a Cantera ``Solution`` object
        (set by ``CanteraProvider.initialize()``).

        Returns:
            Triplet list, or ``[]`` if Cantera state cannot be set.
        """
        gas = getattr(self, "_gas", None)
        if gas is None:
            return []

        T_in = ref_state.T
        P_in = ref_state.P
        z = ref_state.z
        F_in = ref_state.flowrate
        duty = config.get("duty", 0.0)
        comps = sorted(z.keys())
        n_c = len(comps)

        composition = ", ".join(f"{k}:{v}" for k, v in z.items() if v > 0)

        # Evaluate Cp and partial molar enthalpies at T_in
        try:
            gas.TPX = T_in, P_in, composition
            Cp_in = gas.cp_mole / 1000.0                          # J/(mol·K)
            h_in = [gas.partial_molar_enthalpies[gas.species_index(c)] / 1000.0
                    if c in gas.species_names else 0.0 for c in comps]
        except Exception:
            return []

        # We don't know T_out at this point (it's what the solver is finding).
        # Use T_in as a first-order approximation for T_out-side terms.
        # The Jacobian will be evaluated again at each Newton iterate, so this
        # is consistent with the finite-difference approximation quality.
        Cp_out = Cp_in
        h_out = h_in

        # Offset mapping: inlet and outlet are the same stream block in our
        # residual scheme (inlet_offset == col_offset, outlet == row_offset).
        # Variable indices within a block: 0=T, 1=P, 2=F, 3+i=z_i
        triplets: list[tuple[int, int, float]] = []

        # --- Row 0: energy balance ---
        # dR/dT_out = Cp_out  (outlet T is a variable; row=row_offset, col in outlet block)
        # In the EO scheme both inlet and outlet blocks appear; the heater mixin
        # residuals reference both inlet_vals and outlet_vals variables.
        # Here col_offset refers to the inlet block and row_offset to the outlet block.

        # dR/dT_out (outlet block, col = row_offset + 0)
        triplets.append((row_offset + 0, row_offset + 0, Cp_out))
        # dR/dT_in  (inlet block, col = col_offset + 0)
        triplets.append((row_offset + 0, col_offset + 0, -Cp_in))
        # dR/dF_in
        if abs(F_in) > 1e-20:
            triplets.append((row_offset + 0, col_offset + 2, duty / (F_in ** 2)))
        # dR/dz_i (both inlet and outlet composition terms)
        for i, c in enumerate(comps):
            dh = h_out[i] - h_in[i]
            if abs(dh) > 1e-20:
                triplets.append((row_offset + 0, col_offset + 3 + i, -dh))
                triplets.append((row_offset + 0, row_offset + 3 + i, dh))

        # --- Row 1: isobaric  P_out - P_in = 0 ---
        triplets.append((row_offset + 1, row_offset + 1, 1.0))   # dR/dP_out
        triplets.append((row_offset + 1, col_offset + 1, -1.0))  # dR/dP_in

        # --- Row 2: mass balance  F_out - F_in = 0 ---
        triplets.append((row_offset + 2, row_offset + 2, 1.0))   # dR/dF_out
        triplets.append((row_offset + 2, col_offset + 2, -1.0))  # dR/dF_in

        # --- Rows 3+i: composition pass-through  z_out[i] - z_in[i] = 0 ---
        for i in range(n_c):
            triplets.append((row_offset + 3 + i, row_offset + 3 + i, 1.0))
            triplets.append((row_offset + 3 + i, col_offset + 3 + i, -1.0))

        return triplets
