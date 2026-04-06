"""BaseJacobianMixin: shared utilities for provider Jacobian contributions.

This module provides the base mixin class that concrete provider Jacobian
implementations should subclass.  It supplies:

* A safe default ``compute_jacobian_block`` returning ``[]`` (FD fallback).
* ``_local_fd_block`` — unit-local finite differences through an arbitrary
  ``run_fn``, far cheaper than global FD on the full flowsheet.
* ``_ref_state_to_inlet`` — convert a ReferenceState snapshot to the standard
  stream inlet dict format.
* ``_pack_triplets`` — map local (outlet_j, inlet_i) indices to global
  (row, col) coordinates.

Extension pattern
-----------------
To add Jacobian support to a new provider engine::

    # my_engine_jacobian.py
    from processforge.providers.base_jacobian_mixin import BaseJacobianMixin

    class MyEngineJacobianMixin(BaseJacobianMixin):
        def compute_jacobian_block(self, unit_type, config, ref_state,
                                   row_offset, col_offset, n_vars):
            if unit_type not in ("MyUnit", "MyOtherUnit"):
                return []   # global FD handles all other unit types
            return self._local_fd_block(
                config, ref_state, row_offset, col_offset, n_vars,
                run_fn=lambda cfg, inlet: self._run_my_unit(cfg, inlet),
            )

    # my_provider.py
    from processforge.providers.base import AbstractProvider
    from .my_engine_jacobian import MyEngineJacobianMixin
    from .registry import register_provider

    class MyProvider(AbstractProvider, MyEngineJacobianMixin):
        ...   # implement initialize, get_thermo_properties, compute_unit, teardown

    register_provider("my_engine", MyProvider)

Providers that have no cheap analytic Jacobian (e.g. Monte Carlo engines like
OpenMC) can subclass ``BaseJacobianMixin`` without overriding anything: the
default ``compute_jacobian_block`` returns ``[]`` for all unit types, so the
global FD path is used transparently with no code change anywhere else.
"""
from __future__ import annotations

from copy import deepcopy
from typing import Callable, Optional

from .jacobian_contributor import ReferenceState


class BaseJacobianMixin:
    """Base mixin supplying shared Jacobian utilities for provider implementations.

    Subclass alongside ``AbstractProvider``::

        class MyProvider(AbstractProvider, MyJacobianMixin): ...

    where ``MyJacobianMixin(BaseJacobianMixin)`` overrides
    ``compute_jacobian_block`` for the unit types the engine handles.
    """

    # ------------------------------------------------------------------
    # JacobianContributor protocol method (default: FD fallback)
    # ------------------------------------------------------------------

    def compute_jacobian_block(
        self,
        unit_type: str,
        config: dict,
        ref_state: ReferenceState,
        row_offset: int,
        col_offset: int,
        n_vars: int,
    ) -> list[tuple[int, int, float]]:
        """Default: return [] for all unit types — caller uses global FD.

        Override in concrete mixins to supply analytic or semi-analytic
        Jacobian triplets for specific unit types.
        """
        return []

    # ------------------------------------------------------------------
    # Shared utilities
    # ------------------------------------------------------------------

    def _local_fd_block(
        self,
        config: dict,
        ref_state: ReferenceState,
        row_offset: int,
        col_offset: int,
        n_vars: int,
        run_fn: Callable[[dict, dict], Optional[dict]],
        epsilon: float = 1e-5,
    ) -> list[tuple[int, int, float]]:
        """Unit-local finite differences through ``run_fn``.

        Perturbs each inlet variable (T, P, flowrate, z_0 … z_{n_c-1}) one at
        a time and calls ``run_fn(config, perturbed_inlet)`` to get the
        perturbed outlet.  Returns sparse ``(row, col, value)`` triplets for
        all non-negligible derivatives.

        This is O(n_vars) calls to *one specific engine function* rather than
        O(n_vars × n_units) calls to the global ``evaluate_F`` — typically a
        k× speedup for a flowsheet with k units.

        Args:
            config:     Unit configuration dict.
            ref_state:  Inlet stream snapshot at the current Newton iterate.
            row_offset: Global row index of this unit's first outlet residual.
            col_offset: Global column index of the inlet stream's first variable.
            n_vars:     Number of variables per stream (3 + n_components).
            run_fn:     Callable ``(config, inlet_dict) → outlet_dict | None``.
                        Must return a dict with keys T, P, flowrate, z:{comp:frac}
                        or ``None`` if the run failed (block is skipped).
            epsilon:    Finite-difference step size.

        Returns:
            List of ``(global_row, global_col, value)`` triplets.
            Empty list if ``run_fn`` returns ``None`` at the base point.
        """
        inlet_base = self._ref_state_to_inlet(ref_state)
        outlet_base = run_fn(config, inlet_base)
        if outlet_base is None:
            return []

        comps = sorted(ref_state.z.keys())
        triplets: list[tuple[int, int, float]] = []

        for i, (key, sub_key) in enumerate(_inlet_var_order(comps)):
            inlet_pert = _perturb_inlet(inlet_base, key, sub_key, epsilon)
            outlet_pert = run_fn(config, inlet_pert)
            if outlet_pert is None:
                continue
            diffs = _outlet_diff(outlet_base, outlet_pert, comps, epsilon)
            for j, val in enumerate(diffs):
                if abs(val) > 1e-20:
                    triplets.append((row_offset + j, col_offset + i, val))

        return triplets

    @staticmethod
    def _ref_state_to_inlet(ref_state: ReferenceState) -> dict:
        """Convert a ReferenceState snapshot to the standard inlet stream dict.

        Returns:
            ``{"T": float, "P": float, "flowrate": float, "z": {comp: frac}}``
        """
        return {
            "T": ref_state.T,
            "P": ref_state.P,
            "flowrate": ref_state.flowrate,
            "z": dict(ref_state.z),
        }

    @staticmethod
    def _pack_triplets(
        outlet_var_idx: int,
        inlet_var_idx: int,
        value: float,
        row_offset: int,
        col_offset: int,
    ) -> tuple[int, int, float]:
        """Map local outlet/inlet indices to a global (row, col, value) triplet.

        Variable index order within a stream block:
            0  → T
            1  → P
            2  → flowrate
            3+ → z[comp_i] (alphabetical)

        Args:
            outlet_var_idx: Index of the outlet variable (0-based within block).
            inlet_var_idx:  Index of the inlet variable (0-based within block).
            value:          Derivative value dF_outlet[j] / dx_inlet[i].
            row_offset:     Global offset of the outlet stream block.
            col_offset:     Global offset of the inlet stream block.

        Returns:
            ``(global_row, global_col, value)`` ready to assemble into the
            global sparse Jacobian.
        """
        return (row_offset + outlet_var_idx, col_offset + inlet_var_idx, value)


# ---------------------------------------------------------------------------
# Module-level helpers (private — used only by _local_fd_block)
# ---------------------------------------------------------------------------

def _inlet_var_order(comps: list[str]):
    """Yield ``(key, sub_key)`` for the inlet variable sequence.

    Order: T, P, flowrate, z[comp_0], z[comp_1], …  (comps already sorted).
    ``sub_key`` is ``None`` for scalar fields and the component name for ``z``.
    """
    yield ("T", None)
    yield ("P", None)
    yield ("flowrate", None)
    for c in comps:
        yield ("z", c)


def _perturb_inlet(inlet: dict, key: str, sub_key: Optional[str], epsilon: float) -> dict:
    """Return a deep copy of ``inlet`` with one variable perturbed by ``epsilon``."""
    p = deepcopy(inlet)
    if sub_key is None:
        p[key] = p[key] + epsilon
    else:
        p[key][sub_key] = p[key][sub_key] + epsilon
    return p


def _outlet_diff(
    base: dict, pert: dict, comps: list[str], epsilon: float
) -> list[float]:
    """Return ``(pert - base) / epsilon`` in stream variable order [T, P, F, z…]."""
    diffs: list[float] = [
        (pert["T"] - base["T"]) / epsilon,
        (pert["P"] - base["P"]) / epsilon,
        (pert["flowrate"] - base["flowrate"]) / epsilon,
    ]
    for c in comps:
        diffs.append((pert["z"].get(c, 0.0) - base["z"].get(c, 0.0)) / epsilon)
    return diffs
