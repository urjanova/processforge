"""PFR (Plug Flow Reactor) unit — requires CanteraProvider."""
from __future__ import annotations

from loguru import logger

from .provider_mixin import ProviderMixin


class PFR(ProviderMixin):
    """Plug Flow Reactor unit operation.

    Actual reactor integration is performed by the attached ``CanteraProvider``
    via ``provider.compute_unit("PFR", self.params, inlet)``.

    If no Cantera provider is attached, ``_run_impl`` issues a warning and
    returns the inlet stream unchanged.

    JSON configuration example::

        "R-201": {
            "type": "PFR",
            "provider": "cantera",
            "in": "feed",
            "out": "products",
            "residence_time": 5.0,
            "material": 1
        }
    """

    def __init__(self, name: str, residence_time: float = 1.0, **kwargs):
        self.name = name
        self.params = {"residence_time": residence_time, **kwargs}

    def _run_impl(self, inlet: dict) -> dict:
        logger.warning(
            f"PFR '{self.name}': no CanteraProvider attached. "
            "Returning inlet stream unchanged. "
            "Declare a 'cantera' provider and set \"provider\": \"cantera\" on this unit."
        )
        return dict(inlet)
