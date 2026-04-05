"""CSTR (Continuous Stirred Tank Reactor) unit — requires CanteraProvider."""
from __future__ import annotations

from loguru import logger

from .provider_mixin import ProviderMixin


class CSTR(ProviderMixin):
    """Continuous Stirred Tank Reactor unit operation.

    Actual reactor integration is performed by the attached ``CanteraProvider``
    via ``provider.compute_unit("CSTR", self.params, inlet)``.

    If no Cantera provider is attached (e.g. the unit is used without declaring
    a ``cantera`` provider), ``_run_impl`` issues a warning and returns the
    inlet stream unchanged — useful for topology validation without Cantera
    installed.

    JSON configuration example::

        "R-101": {
            "type": "CSTR",
            "provider": "cantera",
            "in": "feed",
            "out": "products",
            "residence_time": 10.0,
            "material": 1
        }
    """

    def __init__(self, name: str, residence_time: float = 1.0, **kwargs):
        self.name = name
        self.params = {"residence_time": residence_time, **kwargs}

    def _run_impl(self, inlet: dict) -> dict:
        logger.warning(
            f"CSTR '{self.name}': no CanteraProvider attached. "
            "Returning inlet stream unchanged. "
            "Declare a 'cantera' provider and set \"provider\": \"cantera\" on this unit."
        )
        return dict(inlet)
