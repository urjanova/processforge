"""Storage abstraction for CAE cases.

The CLI only ever talks to a :class:`CAERepository`. The default implementation
(:class:`~processforge.cae.file_repository.FileCAERepository`) reads a local
``cae/`` directory, but the interface leaves room for a future remote-backed
repository (database or HTTP service) without touching call sites.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from .types import CAECase, CAECaseRef, CAECaseLink


class CAERepository(ABC):
    """Read-only access to CAE cases and flowsheet links."""

    @abstractmethod
    def list_cases(self) -> list[CAECaseRef]:
        """Return every case version available in this repository."""

    @abstractmethod
    def load_case(self, case_id: str, version: str | None = None) -> CAECase:
        """Load one case, selecting the highest version when *version* is None."""

    @abstractmethod
    def links_for_flowsheet(self, flowsheet_path: str) -> list[CAECaseLink]:
        """Return the case links declared for *flowsheet_path*.

        Implementations should raise ``CAEError`` when the flowsheet is not in
        the index.
        """
