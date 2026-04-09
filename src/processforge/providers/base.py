"""AbstractProvider: base class for all Processforge provider implementations."""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from processforge.types import MaterialDef, SimulationResult, UnitConfig


class AbstractProvider(ABC):
    """
    Bridge between Processforge flowsheet declarations and calculation engines.

    A provider translates stream/unit resource definitions into calls against a
    specific backend (CoolProp, Cantera, OpenModelica FMU, etc.).

    Lifecycle
    ---------
    1. ``initialize(provider_config, flowsheet_config)`` — called once at
       flowsheet startup to load mechanisms, compile FMUs, etc.
    2. ``get_thermo_properties(stream)`` / ``compute_unit(...)`` — called per
       simulation step.
    3. ``teardown()`` — called after the simulation to release resources.

    Context manager protocol is supported so providers can be used with
    ``with`` statements::

        with provider:
            results = flowsheet.run()
    """

    @abstractmethod
    def initialize(self, provider_config: dict, flowsheet_config: dict) -> None:
        """Set up provider resources from the flowsheet configuration.

        Args:
            provider_config: The provider's own config block from the JSON
                ``providers`` map (e.g. ``{"type": "cantera", "mechanism": "gri30.yaml"}``).
            flowsheet_config: The full validated flowsheet config dict.
        """

    @abstractmethod
    def get_thermo_properties(self, stream: dict) -> dict:
        """Compute thermodynamic properties for a stream snapshot.

        Args:
            stream: ``{"T": float, "P": float, "z": {comp: frac}}``

        Returns:
            ``{"H": float [J/mol], "Cp": float [J/(mol·K)], "K_values": {comp: float}}``
        """

    @abstractmethod
    def compute_unit(
        self,
        unit_type: str,
        config: dict,
        inlet: dict,
    ) -> Optional[dict]:
        """Optionally compute a unit operation outlet.

        Return a fully populated outlet stream dict if this provider handles
        the given unit type, or ``None`` to fall through to the default
        sequential-modular unit logic.

        Args:
            unit_type: Unit class name (``"Pump"``, ``"CSTR"``, etc.).
            config:    Unit params dict (``self.params`` on the unit object).
            inlet:     Inlet stream ``{"T", "P", "flowrate", "z"}``.

        Returns:
            Outlet stream dict, or ``None`` to delegate to default logic.
        """

    @abstractmethod
    def teardown(self) -> None:
        """Release managed resources (FMU instances, Cantera solutions, etc.)."""

    # ------------------------------------------------------------------
    # Simulation provider interface (optional — override for FEM/neutronics)
    # ------------------------------------------------------------------

    def run_simulation(
        self,
        unit_config: "UnitConfig",
        inlet: dict,
    ) -> "SimulationResult":
        """Run a provider-native standalone simulation (FEM, neutronics, etc.).

        Called by ``SolverUnit._run_impl()``. Receives a typed ``UnitConfig``
        (including ``sim_type`` and ``solver_config``) and the inlet stream
        state dict. Returns a typed ``SimulationResult``.

        Providers that only handle stream-transforming units (CoolProp, Cantera)
        need not override — the default raises ``NotImplementedError`` so the
        error is clear if a ``SolverUnit`` is accidentally paired with them.

        Args:
            unit_config: Typed unit definition (sim_type, material, solver_config, …).
            inlet:       Current inlet stream state dict.

        Returns:
            A :class:`~processforge.types.SimulationResult` with status and any scalar outputs.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement run_simulation. "
            "Use SolverUnit only with providers that support standalone simulations "
            "(e.g. festim, openmc)."
        )

    @classmethod
    def validate_material(
        cls,
        mat_name: str,
        mat_def: "MaterialDef",
        unit_cfg: "UnitConfig",
    ) -> list:
        """Return validation error strings for a material used by this provider.

        Called by the flowsheet validator for every unit backed by this provider.
        Returning an empty list means the material is valid for this provider.

        Override in each provider to enforce provider-specific required fields.
        For example, Festim requires D_0 and E_D; OpenMC requires nuclides.

        Args:
            mat_name: Key from the flowsheet ``materials`` section.
            mat_def:  Typed :class:`~processforge.types.MaterialDef` for this material.
            unit_cfg: Typed :class:`~processforge.types.UnitConfig` for the unit referencing it.

        Returns:
            List of human-readable error strings. Empty list = material is valid.
        """
        return []

    # ------------------------------------------------------------------
    # Context manager support
    # ------------------------------------------------------------------

    def __enter__(self) -> "AbstractProvider":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.teardown()
        return False
