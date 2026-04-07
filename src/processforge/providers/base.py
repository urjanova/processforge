"""AbstractProvider: base class for all Processforge provider implementations."""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional


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
    # Context manager support
    # ------------------------------------------------------------------

    def __enter__(self) -> "AbstractProvider":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.teardown()
        return False
