"""ProviderMixin — intercepts unit.run() to allow provider-computed outlets."""
from __future__ import annotations


class ProviderMixin:
    """Mixin that routes ``run(inlet)`` through the attached provider.

    Before calling the unit's own ``_run_impl()``, this mixin asks the
    attached ``_provider`` whether it wants to fully handle this unit type
    (e.g. ``CanteraProvider`` handling a ``CSTR``).  If the provider returns
    a non-``None`` outlet dict, that result is used and ``_run_impl`` is
    skipped entirely.  If it returns ``None``, ``_run_impl`` runs as normal.

    No provider attached (``_provider`` absent or ``None``) → always calls
    ``_run_impl`` directly, preserving 100 % backward compatibility.
    """

    def run(self, inlet: dict) -> dict:
        provider = getattr(self, "_provider", None)
        if provider is not None:
            result = provider.compute_unit(
                type(self).__name__,
                getattr(self, "params", {}),
                inlet,
            )
            if result is not None:
                return result
        return self._run_impl(inlet)

    def _run_impl(self, inlet: dict) -> dict:  # pragma: no cover
        raise NotImplementedError(
            f"{type(self).__name__} must implement _run_impl()"
        )
