"""Provider abstraction layer for ProcessForge.

Providers are the 'how' behind flowsheet declarations — they bridge ProcessForge
stream/unit definitions to underlying calculation engines (CoolProp, Cantera,
OpenModelica, etc.).

Declare providers in the flowsheet JSON::

    {
      "providers": {
        "cantera": {"type": "cantera", "mechanism": "gri30.yaml"}
      },
      "default_provider": "cantera",
      "units": {
        "R-101": {"type": "CSTR", "provider": "cantera", ...}
      }
    }

Omitting the ``providers`` block (or a unit-level ``"provider"`` key) falls
back to the built-in CoolProp provider — existing flowsheets are unaffected.
"""

from .base import AbstractProvider
from .coolprop_provider import CoolPropProvider
from .registry import get_provider_class, register_provider

__all__ = [
    "AbstractProvider",
    "CoolPropProvider",
    "get_provider_class",
    "register_provider",
]
