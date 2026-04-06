"""ReferenceStateRegistry: per-Jacobian-evaluation stream snapshot store.

The registry is created once at the start of each ``GlobalJacobianManager.
evaluate_J()`` call and discarded when the call returns.  Its job is to
guarantee T/P consistency across providers:

1. Before any provider is asked for Jacobian contributions, the manager calls
   ``snapshot_stream()`` for every registered stream, recording (T, P, z,
   flowrate) from the current x vector.

2. When a provider's ``compute_jacobian_block()`` needs the inlet state, it
   receives the same ``ReferenceState`` object that every other provider for
   that stream receives.  No provider can observe a stale or different T/P.

3. If two providers independently report different T or P for the same stream
   (e.g. because they each maintain internal state), ``validate_agreement()``
   emits a ``UserWarning`` so the inconsistency is surfaced.

Lifecycle
---------
::

    # Inside GlobalJacobianManager.evaluate_J():
    registry = ReferenceStateRegistry()

    for name, sv in self._streams.items():
        vals = sv.to_vector_slice(x)
        registry.snapshot_stream(name, vals["T"], vals["P"], vals["z"],
                                 vals["flowrate"], sv.global_offset)

    # Pass registry into provider calls:
    ref_state = registry.get(inlet_stream_name)
    triplets = provider.compute_jacobian_block(..., ref_state, ...)

    # Optionally cross-check two providers:
    registry.validate_agreement(stream_name,
                                "cantera", T_cantera, P_cantera,
                                "modelica", T_modelica, P_modelica)
"""
from __future__ import annotations

import warnings

from .jacobian_contributor import ReferenceState


class ReferenceStateRegistry:
    """Per-``evaluate_J()`` snapshot store for stream reference states.

    All providers evaluating a stream receive the same ``ReferenceState``
    object, ensuring they agree on the T/P operating point.

    Args:
        tol_T: Temperature agreement tolerance [K].  Differences larger than
               this trigger a ``UserWarning`` from ``validate_agreement()``.
        tol_P: Pressure agreement tolerance [Pa].
    """

    def __init__(self, tol_T: float = 0.1, tol_P: float = 1.0) -> None:
        self._states: dict[str, ReferenceState] = {}
        self._tol_T = tol_T
        self._tol_P = tol_P

    # ------------------------------------------------------------------
    # Snapshot API
    # ------------------------------------------------------------------

    def snapshot_stream(
        self,
        stream_name: str,
        T: float,
        P: float,
        z: dict[str, float],
        flowrate: float,
        global_offset: int,
    ) -> ReferenceState:
        """Record the current state of a stream from the x vector.

        Calling this method a second time for the same ``stream_name``
        overwrites the previous snapshot (idempotent when x has not changed).

        Args:
            stream_name:   Name of the stream in the flowsheet topology.
            T:             Temperature [K] at the current Newton iterate.
            P:             Pressure [Pa] at the current Newton iterate.
            z:             Mole-fraction dict ``{component: fraction}``.
            flowrate:      Molar flowrate [mol/s].
            global_offset: ``StreamVar.global_offset`` for this stream.

        Returns:
            The newly created ``ReferenceState`` snapshot.
        """
        state = ReferenceState(
            stream_name=stream_name,
            T=T,
            P=P,
            z=dict(z),
            flowrate=flowrate,
            global_offset=global_offset,
        )
        self._states[stream_name] = state
        return state

    def get(self, stream_name: str) -> ReferenceState:
        """Retrieve a previously snapshotted stream state.

        Args:
            stream_name: Name of the stream.

        Returns:
            The ``ReferenceState`` recorded by the most recent
            ``snapshot_stream()`` call for this stream.

        Raises:
            KeyError: If ``snapshot_stream`` was never called for this stream.
        """
        return self._states[stream_name]

    # ------------------------------------------------------------------
    # Consistency check
    # ------------------------------------------------------------------

    def validate_agreement(
        self,
        stream_name: str,
        prov_a: str,
        T_a: float,
        P_a: float,
        prov_b: str,
        T_b: float,
        P_b: float,
    ) -> None:
        """Warn if two providers report different T or P for the same stream.

        Call this after obtaining T/P values from two providers that both
        claim to evaluate at ``stream_name``.  If the values disagree beyond
        the tolerances set in ``__init__``, a ``UserWarning`` is emitted
        describing the mismatch so the user can diagnose provider state drift.

        Args:
            stream_name: Name of the stream being evaluated.
            prov_a:      Human-readable name of the first provider.
            T_a:         Temperature [K] as reported by provider A.
            P_a:         Pressure [Pa] as reported by provider A.
            prov_b:      Human-readable name of the second provider.
            T_b:         Temperature [K] as reported by provider B.
            P_b:         Pressure [Pa] as reported by provider B.
        """
        if abs(T_a - T_b) > self._tol_T:
            warnings.warn(
                f"ProviderBridge: T mismatch on stream '{stream_name}': "
                f"{prov_a}={T_a:.3f} K vs {prov_b}={T_b:.3f} K "
                f"(delta={abs(T_a - T_b):.3f} K, tol={self._tol_T} K). "
                "Both providers must evaluate at the same reference state.",
                stacklevel=3,
            )
        if abs(P_a - P_b) > self._tol_P:
            warnings.warn(
                f"ProviderBridge: P mismatch on stream '{stream_name}': "
                f"{prov_a}={P_a:.1f} Pa vs {prov_b}={P_b:.1f} Pa "
                f"(delta={abs(P_a - P_b):.1f} Pa, tol={self._tol_P} Pa). "
                "Both providers must evaluate at the same reference state.",
                stacklevel=3,
            )
