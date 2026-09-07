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

from collections.abc import MutableMapping
from typing import Optional

from loguru import logger
from pydantic import BaseModel, ConfigDict

from processforge.types import CoolPropProviderConfig, FlowsheetConfig, ProviderConfig

from .base import AbstractProvider
from .coolprop_provider import CoolPropProvider
from .errors import ProviderNotAvailableError
from .registry import get_provider_class, is_containerized

_BUILTIN_DEFAULT_KEY = "__coolprop__"
_DEFAULT_KEY = "__default__"


class UnitProviderConfig(BaseModel):
    """Thin wrapper around a unit's config dict for provider resolution.

    Extracts only the ``provider`` key, which is all
    :meth:`ProviderMap.resolve` needs.
    """

    model_config = ConfigDict(extra="allow")

    provider: Optional[str] = None


class ProviderMap(MutableMapping[str, AbstractProvider]):
    """Typed, dict-compatible container for initialised providers.

    Implements the :class:`~collections.abc.MutableMapping` protocol so existing
    call sites can use ``__getitem__``, ``__contains__``, ``.values()``,
    iteration, ``len()``, etc.  Adds :meth:`resolve` for provider lookup and a
    ``_default`` fallback.
    """

    def __init__(
        self,
        providers: Optional[dict[str, AbstractProvider]] = None,
        default: Optional[AbstractProvider] = None,
    ) -> None:
        self._providers: dict[str, AbstractProvider] = dict(providers or {})
        self._default: Optional[AbstractProvider] = default

    # -- MutableMapping protocol --------------------------------------------

    def __getitem__(self, key: str) -> AbstractProvider:
        return self._providers[key]

    def __setitem__(self, key: str, value: AbstractProvider) -> None:
        self._providers[key] = value

    def __delitem__(self, key: str) -> None:
        del self._providers[key]

    def __iter__(self):
        return iter(self._providers)

    def __len__(self) -> int:
        return len(self._providers)

    def __repr__(self) -> str:
        names = list(self._providers.keys())
        default_name = self._default_name()
        return f"ProviderMap({names!r}, default={default_name!r})"

    # -- provider resolution ------------------------------------------------

    def resolve(self, unit_config: UnitProviderConfig) -> AbstractProvider:
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

    def _default_name(self) -> Optional[str]:
        """Return the key under which the default provider is stored, if any."""
        if self._default is None:
            return None
        for name, provider in self._providers.items():
            if provider is self._default:
                return name
        return "<external>"


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
        providers: dict[str, AbstractProvider] = {_BUILTIN_DEFAULT_KEY: coolprop}
    except ProviderNotAvailableError as exc:
        logger.debug(f"Skipping implicit CoolProp fallback: {exc}")
        providers = {}
        coolprop = None

    # Step 2: instantiate every declared provider.
    for name, cfg in providers_config.items():
        ptype = cfg.type
        if is_containerized(ptype):
            # Containerized providers (OpenMC, FESTIM, …) run their backend
            # inside a Docker container; the CLI only needs a thin HTTP client
            # that talks to the container's provider_server.py.
            from .container_client import ContainerProviderClient

            instance: AbstractProvider = ContainerProviderClient(ptype)
            logger.info(
                f"Initialized provider '{name}' (type='{ptype}') via container client"
            )
        else:
            cls = get_provider_class(ptype)
            instance = cls()
            logger.info(f"Initialized provider '{name}' (type='{ptype}')")
        instance.initialize(cfg, flowsheet_config)
        providers[name] = instance

    # Step 3: resolve the active default.
    default: Optional[AbstractProvider] = None
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

    return ProviderMap(providers=providers, default=default)


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
