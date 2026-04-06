"""JacobianContributor: opt-in AD endpoint for ProcessForge physics providers.

Any provider (Cantera, Modelica/FMU, PathSim, OpenMC, etc.) that can supply
its own Jacobian contributions — rather than relying on global finite differences
— implements the JacobianContributor protocol.

The AD endpoint contract
------------------------
A provider announces Jacobian capability by having a ``compute_jacobian_block``
method that matches the protocol below.  The check is::

    isinstance(provider, JacobianContributor)   # True if method is present

Returning an empty list from ``compute_jacobian_block`` is the escape hatch:
it signals "I don't handle this unit type; fall back to global FD for this block."

ReferenceState
--------------
Passed into ``compute_jacobian_block`` by the GlobalJacobianManager.  It carries
the (T, P, z, flowrate) snapshot for the inlet stream *at the current Newton
iterate x*, created by ReferenceStateRegistry before any provider is called.
All providers evaluating the same stream therefore see the same reference point —
this is the T/P consistency guarantee.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable


@dataclass
class ReferenceState:
    """Frozen snapshot of a stream at the current Newton iterate.

    Created once per ``evaluate_J()`` call per stream by
    ``ReferenceStateRegistry.snapshot_stream()``, then passed unchanged
    to every provider that contributes to units fed by that stream.

    Attributes:
        stream_name:   Name of the stream in the flowsheet topology.
        T:             Temperature [K].
        P:             Pressure [Pa].
        z:             Mole-fraction dict ``{component: fraction}``.
        flowrate:      Molar flowrate [mol/s].
        global_offset: Position of this stream's variable block in the
                       global x vector (== StreamVar.global_offset).
    """

    stream_name: str
    T: float
    P: float
    z: dict[str, float] = field(default_factory=dict)
    flowrate: float = 1.0
    global_offset: int = -1


@runtime_checkable
class JacobianContributor(Protocol):
    """Protocol for providers that supply analytic Jacobian contributions.

    Implementing this protocol is **opt-in**.  Providers that do not implement
    it (e.g. CoolPropProvider) are never checked against it and always use the
    global finite-difference path.

    A provider implementing this protocol must define ``compute_jacobian_block``
    with the exact signature below.  No base class is required — the
    ``@runtime_checkable`` decorator allows ``isinstance`` checks without
    inheritance.

    For convenience, providers should subclass ``BaseJacobianMixin`` (from
    ``processforge.providers.base_jacobian_mixin``) which:
    - Provides a default ``compute_jacobian_block`` returning ``[]`` (safe FD
      fallback).
    - Provides ``_local_fd_block``, ``_ref_state_to_inlet``, and
      ``_pack_triplets`` shared utilities to avoid reimplementing indexing math.

    Extension pattern for a new provider::

        # my_engine_jacobian.py
        from processforge.providers.base_jacobian_mixin import BaseJacobianMixin

        class MyEngineJacobianMixin(BaseJacobianMixin):
            def compute_jacobian_block(self, unit_type, config, ref_state,
                                       row_offset, col_offset, n_vars):
                if unit_type not in ("MyUnit",):
                    return []
                return self._local_fd_block(
                    config, ref_state, row_offset, col_offset, n_vars,
                    run_fn=lambda cfg, inlet: self._run_my_unit(cfg, inlet),
                )

        # my_engine_provider.py
        from processforge.providers.base import AbstractProvider
        from .my_engine_jacobian import MyEngineJacobianMixin
        from .registry import register_provider

        class MyEngineProvider(AbstractProvider, MyEngineJacobianMixin):
            ...

        register_provider("my_engine", MyEngineProvider)
    """

    def compute_jacobian_block(
        self,
        unit_type: str,
        config: dict,
        ref_state: ReferenceState,
        row_offset: int,
        col_offset: int,
        n_vars: int,
    ) -> list[tuple[int, int, float]]:
        """Return sparse Jacobian contributions for one unit block.

        Args:
            unit_type:   Unit class name (e.g. ``"CSTR"``, ``"Heater"``).
            config:      Unit params dict from the flowsheet JSON.
            ref_state:   Snapshot of the inlet stream at the current x.
                         The same object is passed to every provider evaluating
                         this stream — do not mutate it.
            row_offset:  Global row offset for this unit's outlet residuals
                         (== outlet ``StreamVar.global_offset``).
            col_offset:  Global column offset for the inlet variables
                         (== inlet ``StreamVar.global_offset``).
            n_vars:      Variables per stream: ``3 + n_components``
                         (T, P, flowrate, z_0 … z_{n_c-1}).

        Returns:
            A list of ``(global_row, global_col, value)`` triplets that will
            be assembled into the global sparse Jacobian matrix.

            Return an **empty list** to signal that this provider does not
            handle the given ``unit_type`` — the caller falls back to global
            finite differences for this block.
        """
        ...
