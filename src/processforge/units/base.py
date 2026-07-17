"""Base mixin for all process unit models.

Provides the standard ``run()`` / ``_run_impl()`` contract used by the
sequential-modular (SM) flowsheet engine, including automatic provider
delegation when a provider (e.g. ``CanteraProvider``) is attached.
"""


class BaseUnitMixin:
    """Base mixin for all SM unit models.

    ``run(inlet)`` is the entry point called by the flowsheet.  If a
    ``_provider`` is attached, it gets first chance to handle the
    computation via ``provider.compute_unit()``.  If it returns
    ``None`` (declines), falls through to ``_run_impl(inlet)``.

    Subclasses must override ``_run_impl(inlet)`` with their own SM
    logic.
    """

    def run(self, inlet: dict) -> dict:
        """Process *inlet* stream, delegating to provider if attached."""
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
        """Override in subclass with unit-specific SM logic."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement _run_impl()"
        )
