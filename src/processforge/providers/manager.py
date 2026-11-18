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

from processforge.types import CoolPropProviderConfig, FlowsheetConfig, ProviderConfig, RunContext

from .base import AbstractProvider
from .coolprop_provider import CoolPropProvider
from .errors import ProviderNotAvailableError
from .registry import get_provider_class, is_containerized

_BUILTIN_DEFAULT_KEY = "__coolprop__"
_DEFAULT_KEY = "__default__"


class UnitProviderConfig(BaseModel):
    """Thin wrapper around a unit's config dict for provider resolution.

    Extracts the ``provider`` and ``thermo_provider`` keys, which are used by
    :meth:`ProviderMap.resolve` and :meth:`ProviderMap.thermo_resolve` respectively.
    """

    model_config = ConfigDict(extra="allow")

    provider: Optional[str] = None
    thermo_provider: Optional[str] = None


class ProviderMap(MutableMapping[str, AbstractProvider]):
    """Typed, dict-compatible container for initialised providers.

    Implements the :class:`~collections.abc.MutableMapping` protocol so existing
    call sites can use ``__getitem__``, ``__contains__``, ``.values()``,
    iteration, ``len()``, etc.  Adds :meth:`resolve` for provider lookup and a
    ``_default`` fallback.

    Additionally supports a separate ``thermo_default`` provider used for stream
    thermo property computation (``H``, ``Cp``, ``K_values``).  This is always
    CoolProp when declared; units/streams may override with ``thermo_provider``.
    """

    def __init__(
        self,
        providers: Optional[dict[str, AbstractProvider]] = None,
        default: Optional[AbstractProvider] = None,
        thermo_default: Optional[AbstractProvider] = None,
    ) -> None:
        self._providers: dict[str, AbstractProvider] = dict(providers or {})
        self._default: Optional[AbstractProvider] = default
        self._thermo_default: Optional[AbstractProvider] = thermo_default

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
        """Return the provider for a unit's operation, falling back to the default.

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
            f"Unit '{unit_config.unit_name}' omits 'provider' and no default provider "
            f"is configured. Declare 'provider' on the unit or set "
            f"'default_provider' in the flowsheet."
        )

    def thermo_resolve(self, unit_config: UnitProviderConfig) -> AbstractProvider:
        """Return the thermo provider for a unit's outlet stream.

        Falls back to the thermo_default (CoolProp if declared) if no
        explicit ``thermo_provider`` is set on the unit.
        """
        tname = unit_config.thermo_provider
        if tname is not None:
            if tname not in self._providers:
                raise ValueError(
                    f"Unit '{unit_config.unit_name}' references thermo_provider "
                    f"'{tname}', which is not declared in 'providers'. "
                    f"Declared providers: {[k for k in self._providers if not k.startswith('__')]}"
                )
            return self._providers[tname]
        if self._thermo_default is not None:
            return self._thermo_default
        raise ValueError(
            f"Unit '{unit_config.unit_name}' omits 'thermo_provider' and no default "
            f"thermo provider is configured. Declare 'thermo_provider' on the unit "
            f"or ensure 'coolprop' is declared in 'providers'."
        )
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

    The returned :class:`ProviderMap` contains only entries that appear in
    ``providers_config``.  There is no implicit CoolProp fallback; if CoolProp
    is needed it must be declared in the ``providers`` block.

    The provider map has two purposes:

    1. **Unit operation** — ``unit.provider`` selects which backend performs
       hydraulic unit computations (Pump, Valve, Strainer, Pipes).  The default
       is CoolProp when ``default_provider`` is set; otherwise every unit must
       explicitly set ``provider``.
    2. **Stream thermo** — ``unit.thermo_provider`` / ``stream.thermo_provider``
       selects which backend computes ``{H, Cp, K_values}`` for that stream.
       The default is always CoolProp when declared; otherwise every unit/stream
       must explicitly set ``thermo_provider``.

    Args:
        providers_config: Mapping of provider name → typed provider config.
        flowsheet_config: Typed representation of the full flowsheet config.

    Returns:
        Initialised :class:`ProviderMap`.
    """
    # Step 1: instantiate every declared provider.
    providers: dict[str, AbstractProvider] = {}
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

    # Step 2: resolve the unit default.
    # If the user set default_provider, use it; otherwise the default is
    # CoolProp only if it appears in the declared providers.  No silent fallback.
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
    elif "coolprop" in providers:
        default = providers["coolprop"]
        logger.info(f"Default provider set to 'coolprop'")
    else:
        default = None
        logger.debug("No default provider declared; units must set provider explicitly")

    # Step 3: resolve the thermo default.
    # CoolProp is the default thermo backend when declared; otherwise units/
    # streams must set thermo_provider explicitly.
    if "coolprop" in providers:
        thermo_default = providers["coolprop"]
        logger.info("Default thermo provider set to 'coolprop'")
    else:
        thermo_default = None
        logger.debug("No default thermo provider declared; units/streams must set thermo_provider explicitly")

    return ProviderMap(providers=providers, default=default,
                       thermo_default=thermo_default)


def apply_run_context(
    provider_map: ProviderMap | None,
    context: RunContext | None,
) -> None:
    """Propagate run-level metadata to every provider in *provider_map*.

    Providers that do not consume run context (e.g. CoolProp) inherit a no-op
    implementation from ``AbstractProvider``.
    """
    if provider_map is None or context is None:
        return
    for provider in provider_map.values():
        provider.set_run_context(context)


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
