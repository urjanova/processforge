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

from loguru import logger

from .coolprop_provider import CoolPropProvider
from .registry import get_provider_class

_BUILTIN_DEFAULT_KEY = "__coolprop__"
_DEFAULT_KEY = "__default__"


def build_provider_map(
    providers_config: dict,
    flowsheet_config: dict,
) -> dict[str, "AbstractProvider"]:  # noqa: F821
    """Parse the ``providers`` block and return a ready-to-use provider map.

    The returned dict always contains:

    * ``"__coolprop__"`` — the built-in CoolProp fallback (always present).
    * One entry per named provider declared in ``providers_config``.
    * ``"__default__"`` — whichever provider units fall back to when they
      omit the ``"provider"`` key.  Equals ``__coolprop__`` unless the
      flowsheet sets ``"default_provider"``.

    Args:
        providers_config: The ``config["providers"]`` dict (may be empty).
        flowsheet_config: The full validated flowsheet config dict.

    Returns:
        Mapping of provider name → initialised ``AbstractProvider`` instance.
    """
    # Step 1: always seed with the built-in CoolProp fallback.
    coolprop = CoolPropProvider()
    coolprop.initialize({}, flowsheet_config)
    provider_map: dict = {_BUILTIN_DEFAULT_KEY: coolprop}

    # Step 2: instantiate every declared provider.
    for name, cfg in providers_config.items():
        ptype = cfg.get("type")
        if not ptype:
            raise ValueError(
                f"Provider '{name}' is missing the required 'type' field."
            )
        _maybe_import_provider(ptype)
        cls = get_provider_class(ptype)
        instance = cls()
        instance.initialize(cfg, flowsheet_config)
        provider_map[name] = instance
        logger.info(f"Initialized provider '{name}' (type='{ptype}')")

    # Step 3: resolve the active default.
    user_default = flowsheet_config.get("default_provider")
    if user_default is not None:
        if user_default not in provider_map:
            raise ValueError(
                f"'default_provider' references '{user_default}', which is not "
                f"declared in the 'providers' block. "
                f"Declared providers: {[k for k in provider_map if not k.startswith('__')]}"
            )
        provider_map[_DEFAULT_KEY] = provider_map[user_default]
        logger.info(f"Default provider set to '{user_default}'")
    else:
        provider_map[_DEFAULT_KEY] = coolprop

    return provider_map


def get_unit_provider(
    provider_map: dict,
    unit_config: dict,
) -> "AbstractProvider":  # noqa: F821
    """Resolve the provider for a specific unit.

    Falls back to ``__default__`` (CoolProp, or user-declared default) when
    the unit config contains no ``"provider"`` key.

    Args:
        provider_map: Map returned by :func:`build_provider_map`.
        unit_config:  The unit's config dict from the flowsheet JSON.

    Raises:
        ValueError: If the unit references an undeclared provider name.
    """
    pname = unit_config.get("provider", _DEFAULT_KEY)
    if pname not in provider_map:
        raise ValueError(
            f"Unit references provider '{pname}', which is not declared in "
            f"the 'providers' block. "
            f"Declared providers: {[k for k in provider_map if not k.startswith('__')]}"
        )
    return provider_map[pname]


def teardown_providers(provider_map: dict) -> None:
    """Call ``teardown()`` on every provider in the map.

    Errors during teardown are logged as warnings rather than raised so that
    all providers get a chance to clean up.
    """
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


def _maybe_import_provider(ptype: str) -> None:
    """Lazily import optional provider modules so they self-register."""
    if ptype == "cantera":
        from . import cantera_provider  # noqa: F401
    elif ptype == "modelica":
        from . import modelica_provider  # noqa: F401
    # "coolprop" is already in the registry from registry._seed_registry().
