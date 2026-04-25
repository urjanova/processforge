"""GlobalJacobianManager: assembles and evaluates the global EO residual and Jacobian."""
from __future__ import annotations

from copy import deepcopy

import numpy as np
import scipy.sparse as sp
from loguru import logger

from .stream_var import StreamVar


class GlobalJacobianManager:
    """
    Owns the global variable indexing table and evaluates F(x) and J(x).

    Variable layout in x:
        For each registered stream s (in registration order):
            x[s.global_offset + 0] = T_s
            x[s.global_offset + 1] = P_s
            x[s.global_offset + 2] = F_s
            x[s.global_offset + 3+i] = z_s[comp_i]

    Residual layout in F(x) mirrors the variable layout:
        Feed streams:  F[offset+j] = x[offset+j] - feed_value[j]  (BC)
        Unit outlets:  F[offset+j] = unit.get_scipy_residuals(...)[j]

    Stream identity:  a stream name appearing as ``out`` on unit A and ``in`` on
    unit B resolves to the **same** StreamVar (same global_offset), so connection
    equations are automatically satisfied — no additional equality constraints needed.
    """

    def __init__(self, components: list[str]) -> None:
        self.components: list[str] = components
        self.n_c: int = len(components)
        self.vars_per_stream: int = 3 + self.n_c  # T, P, F, z[0..N_c-1]

        self._streams: dict[str, StreamVar] = {}
        self._feed_values: dict[str, dict] = {}
        self._packed_feed_cache: dict[str, np.ndarray] = {}
        self._unit_entries: list[dict] = []
        self._n_vars: int = 0

        # Cached derived structures — invalidated when streams are registered
        self._cached_column_groups: list[list[int]] | None = None
        self._stream_by_offset: dict[int, str] | None = None
        self._results: dict[str, dict] | None = None

    # ------------------------------------------------------------------
    # Registration API (called by EOFlowsheet.build())
    # ------------------------------------------------------------------

    def register_stream(self, name: str) -> StreamVar:
        """Register a process stream and assign it a unique global offset."""
        if name in self._streams:
            return self._streams[name]
        sv = StreamVar(
            name=name,
            components=self.components,
            global_offset=self._n_vars,
        )
        self._streams[name] = sv
        self._n_vars += self.vars_per_stream
        # Invalidate caches that depend on the stream set
        self._cached_column_groups = None
        self._stream_by_offset = None
        logger.debug(f"Registered stream '{name}' at offset {sv.global_offset}")
        return sv

    def register_feed(self, name: str, feed_dict: dict) -> StreamVar:
        """Register a feed stream with fixed boundary-condition values."""
        sv = self.register_stream(name)
        self._feed_values[name] = deepcopy(feed_dict)
        # Pre-compute packed numpy array; feed values are constant across the solve.
        self._packed_feed_cache[name] = sv.from_stream_dict(feed_dict)
        return sv

    def register_unit(
        self,
        unit,
        inlet_names: list[str],
        outlet_names: list[str],
        config: dict,
    ) -> None:
        """Register a unit and ensure all its streams are registered."""
        for name in inlet_names + outlet_names:
            self.register_stream(name)
        self._unit_entries.append(
            {
                "unit": unit,
                "inlet_names": inlet_names,
                "outlet_names": outlet_names,
                "config": config,
                # Provider reference stored here so evaluate_J() can route
                # JacobianContributor-capable providers to the analytic path.
                "provider": getattr(unit, "_provider", None),
            }
        )

    # ------------------------------------------------------------------
    # Initial guess
    # ------------------------------------------------------------------

    def get_x0(self, stream_initial_values: dict[str, dict]) -> np.ndarray:
        """
        Build the initial guess vector x0 from a dict of stream initial values.

        Streams not present in ``stream_initial_values`` default to
        T=298.15 K, P=101325 Pa, F=1.0 mol/s, z=0 for all components.
        """
        x0 = np.zeros(self._n_vars)
        defaults = {"T": 298.15, "P": 101325.0, "flowrate": 1.0, "z": {}}
        for name, sv in self._streams.items():
            vals = stream_initial_values.get(name, defaults)
            sv.write_to_vector(x0, vals)
        return x0

    # ------------------------------------------------------------------
    # Core evaluation
    # ------------------------------------------------------------------

    def evaluate_F(self, x: np.ndarray) -> np.ndarray:
        """
        Evaluate the global residual vector F(x).

        Returns a numpy array of length ``n_vars``.  The system is solved when
        ``||F(x)||`` is below the solver tolerance.
        """
        F = np.zeros(self._n_vars)

        # 1. Feed (boundary condition) residuals: x[stream] - feed_value = 0
        for name in self._feed_values:
            sv = self._streams[name]
            packed_feed = self._packed_feed_cache[name]
            o = sv.global_offset
            F[o: o + sv.n_vars] = x[o: o + sv.n_vars] - packed_feed

        # 2. Unit model residuals
        for entry in self._unit_entries:
            unit = entry["unit"]
            config = entry["config"]
            inlet_names = entry["inlet_names"]
            outlet_names = entry["outlet_names"]

            inlet_vals = self._merged_inlet_vals(x, inlet_names)

            n_outlets = getattr(unit, "eo_n_outlet_streams", 1)
            if n_outlets > 1:
                # Multi-outlet unit (e.g. Flash): call combined residual method
                all_outlet_vals = [
                    self._streams[n].to_vector_slice(x) for n in outlet_names
                ]
                residuals_all = unit.get_scipy_residuals_multi_outlet(
                    inlet_vals, all_outlet_vals, config
                )
                n_per = len(residuals_all) // len(outlet_names)
                for k, outlet_name in enumerate(outlet_names):
                    sv_out = self._streams[outlet_name]
                    o = sv_out.global_offset
                    for j, r in enumerate(residuals_all[k * n_per:(k + 1) * n_per]):
                        F[o + j] = r
            else:
                for outlet_name in outlet_names:
                    sv_out = self._streams[outlet_name]
                    outlet_vals = sv_out.to_vector_slice(x)
                    residuals = unit.get_scipy_residuals(inlet_vals, outlet_vals, config)
                    o = sv_out.global_offset
                    for j, r in enumerate(residuals):
                        F[o + j] = r

        return F

    def evaluate_J(
        self, x: np.ndarray, epsilon: float = 1e-7
    ) -> sp.csr_matrix:
        """
        Evaluate the global sparse Jacobian J = dF/dx.

        For units whose provider implements the ``JacobianContributor`` protocol,
        analytic (or semi-analytic) triplets are collected directly from the
        provider's ``compute_jacobian_block`` method.  All other stream column
        groups fall back to the existing block-sparse finite-difference path.

        The ``ReferenceStateRegistry`` ensures every provider that contributes
        to the same stream sees the same (T, P, z) snapshot from the current x
        vector, enforcing reference-state consistency across physics engines.
        """
        from processforge.providers.jacobian_contributor import (
            JacobianContributor,
            ReferenceState,
        )
        from processforge.providers.reference_state_registry import ReferenceStateRegistry

        n = self._n_vars
        F0 = self.evaluate_F(x)

        # --- 1. Snapshot every stream into the registry ---
        registry = ReferenceStateRegistry()
        for name, sv in self._streams.items():
            vals = sv.to_vector_slice(x)
            registry.snapshot_stream(
                name,
                vals["T"],
                vals["P"],
                vals["z"],
                vals["flowrate"],
                sv.global_offset,
            )

        rows_list: list[int] = []
        cols_list: list[int] = []
        data_list: list[float] = []

        # --- 2. Analytic path for JacobianContributor providers ---
        # partial_analytic tracks whether EVERY unit fed by a given inlet
        # stream returned non-empty triplets.  Only fully-covered inlets are
        # exempted from the FD pass below.
        partial_analytic: dict[str, bool] = {}

        for entry in self._unit_entries:
            provider = entry.get("provider")
            if provider is None or not isinstance(provider, JacobianContributor):
                # This unit needs FD; mark all its inlets as not fully covered.
                for iname in entry["inlet_names"]:
                    partial_analytic[iname] = False
                continue

            inlet_names = entry["inlet_names"]
            config = entry["config"]

            if len(inlet_names) == 1:
                ref_state = registry.get(inlet_names[0])
                col_offset = self._streams[inlet_names[0]].global_offset
            else:
                # Flow-weighted merged inlet → build synthetic ReferenceState.
                merged = self._merged_inlet_vals(x, inlet_names)
                ref_state = ReferenceState(
                    stream_name="_merged",
                    T=merged["T"],
                    P=merged["P"],
                    z=merged["z"],
                    flowrate=merged["flowrate"],
                    global_offset=-1,
                )
                col_offset = self._streams[inlet_names[0]].global_offset

            for outlet_name in entry["outlet_names"]:
                sv_out = self._streams[outlet_name]
                triplets = provider.compute_jacobian_block(
                    type(entry["unit"]).__name__,
                    config,
                    ref_state,
                    sv_out.global_offset,
                    col_offset,
                    sv_out.n_vars,
                )
                if triplets:
                    for r, c, v in triplets:
                        rows_list.append(r)
                        cols_list.append(c)
                        data_list.append(v)
                    # Mark inlets as covered only if not already negated.
                    for iname in inlet_names:
                        if iname not in partial_analytic:
                            partial_analytic[iname] = True
                else:
                    # Provider returned [] → FD needed for this inlet's columns.
                    for iname in inlet_names:
                        partial_analytic[iname] = False

        # Streams whose inlets are fully covered analytically skip the FD pass.
        analytic_inlet_streams: set[str] = {
            k for k, v in partial_analytic.items() if v
        }

        # --- 3. Finite-difference pass for non-analytic stream groups ---
        if self._stream_by_offset is None:
            self._stream_by_offset = {sv.global_offset: name for name, sv in self._streams.items()}
        if self._cached_column_groups is None:
            self._cached_column_groups = self._sparsity_column_groups()

        for col_group in self._cached_column_groups:
            if not col_group:
                continue
            stream_name = self._stream_by_offset.get(col_group[0])
            if stream_name in analytic_inlet_streams:
                continue  # already covered by analytic path

            # Perturb all columns in the group simultaneously (they don't interact).
            x_pert = x.copy()
            eps_vec = np.zeros(n)
            for col in col_group:
                eps_vec[col] = epsilon
            x_pert += eps_vec

            F_pert = self.evaluate_F(x_pert)
            # Compute dF once per group — identical for every column in the group.
            dF = (F_pert - F0) / epsilon
            nz_rows = np.nonzero(np.abs(dF) > 1e-20)[0].tolist()
            dF_nz = [float(dF[i]) for i in nz_rows]
            for col in col_group:
                rows_list.extend(nz_rows)
                cols_list.extend([col] * len(nz_rows))
                data_list.extend(dF_nz)

        return sp.csr_matrix(
            (data_list, (rows_list, cols_list)), shape=(n, n)
        )

    # ------------------------------------------------------------------
    # Result extraction
    # ------------------------------------------------------------------

    def extract_results(self, x: np.ndarray) -> dict[str, dict]:
        """
        Extract all stream values from the solution vector.

        Returns the same ``{stream_name: {"T": ..., "P": ..., "flowrate": ...,
        "z": {...}}}`` format as ``Flowsheet.run()``.
        """
        results = {
            name: sv.to_vector_slice(x)
            for name, sv in self._streams.items()
        }
        self._results = results
        return results

    @property
    def scalars(self) -> dict[str, dict]:
        """Return scalars collected from SolverUnit runs."""
        if self._results is None:
            return {}
        return self.scalars_from_results(self._results)

    def scalars_from_results(self, results: dict[str, dict]) -> dict[str, dict]:
        """Extract scalars from unit results dict."""
        return {
            k: v
            for k, v in results.items()
            if isinstance(v, dict) and "status" in v
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _merged_inlet_vals(self, x: np.ndarray, inlet_names: list[str]) -> dict:
        """Flow-weighted merge of one or more inlet streams from x vector."""
        if len(inlet_names) == 1:
            return self._streams[inlet_names[0]].to_vector_slice(x)

        merged: dict = {
            "T": 0.0, "P": 0.0, "flowrate": 0.0,
            "z": {c: 0.0 for c in self.components},
        }
        total_F = 0.0
        for name in inlet_names:
            sv = self._streams[name]
            vals = sv.to_vector_slice(x)
            f = vals["flowrate"]
            total_F += f
            merged["T"] += vals["T"] * f
            merged["P"] += vals["P"] * f
            for c in self.components:
                merged["z"][c] += vals["z"].get(c, 0.0) * f

        if total_F > 1e-20:
            merged["T"] /= total_F
            merged["P"] /= total_F
            for c in self.components:
                merged["z"][c] /= total_F
        merged["flowrate"] = total_F
        return merged

    def _sparsity_column_groups(self) -> list[list[int]]:
        """
        Return column groups that can be perturbed simultaneously.

        Two columns can be grouped if they never appear together in the same
        residual block.  Here we use the conservative grouping: one group per
        stream, where columns within the same stream block are perturbed together.
        """
        groups: list[list[int]] = []
        for sv in self._streams.values():
            o = sv.global_offset
            groups.append(list(range(o, o + sv.n_vars)))
        return groups

    @property
    def n_vars(self) -> int:
        """Total number of scalar variables in the global system."""
        return self._n_vars
