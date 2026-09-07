"""BaseSimulationProvider — shared scaffolding for engine-style providers.

Engine providers (OpenMC, FESTIM, and future providers that run standalone
simulations via ``SolverUnit``) share a common lifecycle:

* they do not provide per-stream thermodynamics,
* they do not intercept ``compute_unit()`` calls,
* they keep a material registry built from the flowsheet,
* they resolve a writable run directory with a temp-dir fallback.

This base class captures those defaults so concrete providers only implement
``initialize()`` and ``run_simulation()``.
"""
from __future__ import annotations

import os
import pathlib
import tempfile
from abc import ABC
from typing import TYPE_CHECKING, Optional

from .base import AbstractProvider

if TYPE_CHECKING:
    from processforge.types import FlowsheetConfig, ProviderConfig


class BaseSimulationProvider(AbstractProvider, ABC):
    """Base class for providers that run standalone engine simulations.

    Provides default implementations of ``get_thermo_properties``,
    ``compute_unit``, and ``teardown``, plus a helper for resolving a writable
    run directory.
    """

    def __init__(self) -> None:
        self._materials: dict = {}
        self._provider_output_dir: str = "outputs"
        self._initialized: bool = False

    def initialize(
        self,
        provider_config: "ProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        """Set up the provider. Must be implemented by subclasses."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement initialize()"
        )

    def get_thermo_properties(self, stream: dict) -> dict:
        """Stream thermodynamics are not supported by engine providers."""
        raise NotImplementedError(
            f"{type(self).__name__} does not support stream thermodynamics."
        )

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        """Engine providers use ``run_simulation()`` via ``SolverUnit``."""
        return None

    def teardown(self) -> None:
        """Release provider state."""
        self._initialized = False

    def _resolve_run_dir(self) -> pathlib.Path:
        """Create the run output directory, falling back to a temp dir on failure.

        The configured output dir (often a mounted volume) may not be writable by
        the current user. Fall back to a temp dir so the simulation can still run;
        only the persisted artifacts are lost.
        """
        run_dir = pathlib.Path(self._provider_output_dir)
        try:
            run_dir.mkdir(parents=True, exist_ok=True)
        except (PermissionError, OSError) as exc:
            run_dir = pathlib.Path(tempfile.mkdtemp(prefix="processforge_"))
            from loguru import logger

            logger.warning(
                f"{type(self).__name__}: cannot use output dir "
                f"'{self._provider_output_dir}' ({exc}); "
                f"falling back to '{run_dir}'."
            )
        return run_dir

    def _expand_output_dir(self, output_dir: str) -> str:
        """Expand env vars and anchor relative paths under the run output root."""
        out_dir = os.path.expandvars(output_dir)
        if not os.path.isabs(out_dir):
            root = os.environ.get("PROCESSFORGE_OUTPUT_DIR", "outputs")
            out_dir = os.path.join(root, out_dir)
        return out_dir
