"""Provider lifecycle manager for Processforge flowsheets.

Builds, initialises, and tears down the provider map derived from the
``providers`` block in a flowsheet JSON.  Mirrors the two-level default
pattern used for solver backends:

* No ``providers`` block → CoolProp for everything (zero behavior change).
* ``providers`` block present, unit has no ``"provider"`` key → CoolProp.
* ``"default_provider": "<name>"`` overrides the per-unit fallback.
* ``"provider": "<name>"`` on a unit selects that named provider explicitly.
"""
from __future__ import annotations

from typing import Optional

from loguru import logger
from pydantic import BaseModel, ConfigDict

from processforge.types import CoolPropProviderConfig, FlowsheetConfig, ProviderConfig

from .coolprop_provider import CoolPropProvider
from .registry import get_provider_class

_BUILTIN_DEFAULT_KEY = "__coolprop__"
_DEFAULT_KEY = "__default__"


class UnitProviderConfig(BaseModel):
    """Thin wrapper around a unit's config dict for provider resolution.

    Extracts only the ``provider`` key, which is all
    :meth:`ProviderMap.resolve` needs.
    """

    model_config = ConfigDict(extra="allow")

    provider: Optional[str] = None


class ProviderMap(BaseModel):
    """Typed, dict-compatible container for initialised providers.

    Supports the dict-like operations that existing call sites rely on
    (``__getitem__``, ``__contains__``, ``__bool__``, ``.values()``)
    while adding a :meth:`resolve` method for provider lookup.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    _providers: dict[str, object] = {}
    _default: Optional[object] = None

    def __init__(self, **data):
        super().__init__(**data)
        object.__setattr__(self, "_providers", data.get("_providers", {}))
        object.__setattr__(self, "_default", data.get("_default"))

    # -- dict-like access ---------------------------------------------------

    def __getitem__(self, key: str):
        return self._providers[key]

    def __contains__(self, key: str) -> bool:
        return key in self._providers

    def __bool__(self) -> bool:
        return bool(self._providers)

    def values(self):
        return self._providers.values()

    def keys(self):
        return self._providers.keys()

    def items(self):
        return self._providers.items()

    def get(self, key: str, default=None):
        return self._providers.get(key, default)

    # -- provider resolution ------------------------------------------------

    def resolve(self, unit_config: UnitProviderConfig) -> object:
        """Return the provider for a unit, falling back to the default.

        Raises:
            ValueError: If the unit references an undeclared provider name.
        """
        pname = unit_config.provider
        if pname is not None:
            if pname not in self._providers:
                raise ValueError(
                    f"Unit references provider '{pname}', which is not declared in "
                    f"the 'providers' block. "
                    f"Declared providers: {[k for k in self._providers if not k.startswith('__')]}"
                )
            return self._providers[pname]
        if self._default is not None:
            return self._default
        raise ValueError(
            "No provider specified and no default provider configured."
        )


def build_provider_map(
    providers_config: dict[str, ProviderConfig],
    flowsheet_config: FlowsheetConfig,
) -> ProviderMap:
    """Parse the ``providers`` block and return a ready-to-use provider map.

    The returned :class:`ProviderMap` always contains:

    * ``"__coolprop__"`` — the built-in CoolProp fallback (always present).
    * One entry per named provider declared in ``providers_config``.
    * ``"__default__"`` — whichever provider units fall back to when they
      omit the ``"provider"`` key.  Equals ``__coolprop__`` unless the
      flowsheet sets ``"default_provider"``.

    Args:
        providers_config: Mapping of provider name → typed provider config.
        flowsheet_config: Typed representation of the full flowsheet config.

    Returns:
        Initialised :class:`ProviderMap`.
    """
    # Step 1: always try to seed with the built-in CoolProp fallback.
    coolprop = CoolPropProvider()
    try:
        coolprop.initialize(CoolPropProviderConfig(), flowsheet_config)
        providers: dict = {_BUILTIN_DEFAULT_KEY: coolprop}
    except ImportError as e:
        logger.debug(f"Skipping implicit CoolProp fallback: {e}")
        providers = {}
        coolprop = None

    # Step 2: instantiate every declared provider.
    for name, cfg in providers_config.items():
        ptype = cfg.type
        cls = get_provider_class(ptype)
        instance = cls()
        instance.initialize(cfg, flowsheet_config)
        providers[name] = instance
        logger.info(f"Initialized provider '{name}' (type='{ptype}')")

    # Step 3: resolve the active default.
    default = None
    user_default = flowsheet_config.default_provider
    if user_default is not None:
        if user_default not in providers:
            raise ValueError(
                f"'default_provider' references '{user_default}', which is not "
                f"declared in the 'providers' block. "
                f"Declared providers: {[k for k in providers if not k.startswith('__')]}"
            )
        default = providers[user_default]
        logger.info(f"Default provider set to '{user_default}'")
    elif coolprop is not None:
        default = coolprop

    return ProviderMap(_providers=providers, _default=default)


def teardown_providers(provider_map: ProviderMap | None) -> None:
    """Call ``teardown()`` on every provider in the map.

    Errors during teardown are logged as warnings rather than raised so that
    all providers get a chance to clean up.
    """
    if provider_map is None:
        return
    seen: set[int] = set()
    for name, provider in provider_map.items():
        # Avoid calling teardown twice when __default__ aliases another entry.
        pid = id(provider)
        if pid in seen:
            continue
        seen.add(pid)
        try:
            provider.teardown()
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"Provider '{name}' teardown raised: {exc}")
