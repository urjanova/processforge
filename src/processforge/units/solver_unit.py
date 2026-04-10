"""Generic provider-driven simulation unit (SolverUnit).

SolverUnit is a standalone FEM/neutronics/etc. unit that delegates all
execution to its attached provider via ``provider.run_simulation()``.
It has no mandatory inlet or outlet process stream — it is not a
stream-transforming unit.

Any provider that implements ``run_simulation()`` can be paired with this
unit type (Festim, OpenMC, future solvers).  Provider-specific configuration
lives in the flowsheet JSON ``solver_config`` key and is passed through to
the provider opaquely as part of a :class:`~processforge.types.UnitConfig`.

Flowsheet JSON example::

    "units": {
        "fem_2d": {
            "type": "SolverUnit",
            "provider": "festim",
            "sim_type": "heat_2d",
            "material": 10,
            "solver_config": {
                "mesh": {"nx": 20, "ny": 20},
                "boundary_conditions": [...],
                "traps": [...],
                "exports": [...]
            }
        }
    }
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from loguru import logger

from .provider_mixin import ProviderMixin
from processforge.types import UnitConfig

if TYPE_CHECKING:
    from processforge.types import SimulationResult
    from processforge.providers.base import AbstractProvider


class SolverUnit(ProviderMixin):
    """Generic provider-driven simulation unit.

    Delegates all execution to ``provider.run_simulation()``.  Has no
    mandatory inlet or outlet stream — it is a standalone FEM/neutronics
    solver, not a stream-transforming unit.

    Any provider that implements ``run_simulation()`` can be used here
    (Festim, OpenMC, future solvers).  Provider-specific configuration
    lives in the JSON ``solver_config`` key and is passed through to the
    provider opaquely as part of a typed ``UnitConfig``.

    The last result is stored on ``self._last_result`` for post-run
    introspection.

    JSON config fields:
      sim_type:      str — provider-interpreted simulation mode (e.g. "heat_2d")
      material:      int — ``id`` of the structural material
      solver_config: dict — all provider-specific settings (mesh, BCs, traps, …)
    """

    def __init__(self, name: str, **params):
        self.name = name
        self.params = params
        self._last_result: Optional["SimulationResult"] = None
        # Injected by build_units via ProviderMixin — declared here for type checkers
        self._provider: Optional["AbstractProvider"] = None

    def _run_impl(self, inlet: dict) -> dict:
        """Convert params to a typed ``UnitConfig`` and call ``provider.run_simulation()``.

        Args:
            inlet: Current inlet stream state dict (may be empty for standalone sims).

        Returns:
            Flat dict suitable for flowsheet result storage:
            ``{"status": "completed", "sim_type": ..., <scalar_key>: <value>, …}``

        Raises:
            RuntimeError: If no provider is attached (missing ``"provider"`` key
                in unit config or provider not declared in the flowsheet).
        """
        if not hasattr(self, "_provider") or self._provider is None:
            raise RuntimeError(
                f"SolverUnit '{self.name}' has no provider attached. "
                "Declare a 'provider' key in the unit config pointing to a "
                "provider declared in the flowsheet 'providers' block."
            )

        # build_units filters 'type' and 'material' from **params before passing
        # them to __init__, then sets unit.material separately.  Reconstruct a
        # complete dict so UnitConfig gets all required fields.
        full_params = {
            "type": "SolverUnit",
            **self.params,
        }
        if hasattr(self, "material"):
            full_params.setdefault("material", self.material)

        unit_cfg = UnitConfig.from_dict(full_params)
        logger.info(
            f"SolverUnit '{self.name}': dispatching sim_type='{unit_cfg.sim_type}' "
            f"to provider {type(self._provider).__name__}"
        )

        result = self._provider.run_simulation(unit_cfg, inlet)
        self._last_result = result
        return result.as_dict()
