"""Provider-side error classification.

This module provides:

* A generic exception hierarchy for provider failures:
  :class:`ProviderError`, :class:`ProviderInitError`,
  :class:`ProviderNotAvailableError`, :class:`ProviderRuntimeError`,
  :class:`ProviderCleanupError`, :class:`ProviderValidationError`, and
  :class:`ProviderConfigError`.
* A structured record of a run-time engine failure:
  :class:`ProviderRunError`.
* :func:`classify_run_error` — turn an exception (plus any captured
  stdout/stderr from the engine) into a :class:`ProviderRunError`, picking a
  category and a concrete remediation hint where one is known.
* :func:`make_failed_output` — build the ``EngineOutput(status="failed")`` that
  providers return, populated with the classification.

Categories are intentionally broad so a new engine signature can be added here
without touching the providers.
"""
from __future__ import annotations

import re
from typing import Optional

from pydantic import BaseModel

# ---------------------------------------------------------------------------
# Generic exception hierarchy
# ---------------------------------------------------------------------------


class ProviderError(Exception):
    """Base class for all provider-related failures."""


class ProviderInitError(ProviderError):
    """Raised when a provider fails to initialize.

    Covers missing optional dependencies, unreachable container services,
    bad configuration, and missing data files detected during ``initialize()``.
    """


class ProviderNotAvailableError(ProviderInitError):
    """Raised when a provider backend is not installed or unreachable."""


class ProviderRuntimeError(ProviderError):
    """Raised when a provider fails during a simulation run."""


class ProviderCleanupError(ProviderError):
    """Raised when a provider fails to release resources during teardown.

    These errors are usually logged rather than re-raised so that all
    providers get a chance to clean up.
    """


class ProviderValidationError(ProviderError):
    """Raised when provider-specific material or unit config validation fails."""


class ProviderConfigError(ProviderError):
    """Raised when a provider's ``provider_config`` block is invalid."""


# ---------------------------------------------------------------------------
# Error categories. Strings (not an Enum) so they serialize cleanly in JSON
# through the provider HTTP API and are easy to match on in tests/docs.
# ---------------------------------------------------------------------------

# Generic / shared categories
CONVERGENCE = "convergence"
INPUT_VALIDATION = "input_validation"
ENVIRONMENT = "environment"
UNKNOWN = "unknown"

# OpenMC-specific categories
NUCLEAR_DATA = "nuclear_data"
CROSS_SECTIONS = "cross_sections"
MPI_ABORT = "mpi_abort"
GEOMETRY = "geometry"
TALLY = "tally"

# FESTIM-specific categories
MESH_QUALITY = "mesh_quality"
BC_SETUP = "bc_setup"
SOLVER_CONVERGENCE = "solver_convergence"
MATERIAL_PROPERTY = "material_property"

# ---------------------------------------------------------------------------
# Engine-scoped error signatures
# ---------------------------------------------------------------------------
# Each entry is ``(category, compiled-regex, remediation hint)``. First match
# wins. Engine-specific signatures are checked first, then generic signatures.

_ERROR_SIGNATURES: dict[str, list[tuple[str, re.Pattern[str], str]]] = {
    "openmc": [
        (
            NUCLEAR_DATA,
            re.compile(
                r"nuclear data library does not contain cross sections",
                re.IGNORECASE,
            ),
            (
                "The cross-section library lacks data at the requested "
                "temperature. "
                "Either set the material temperature to a value present in the library "
                "(e.g. 300 K), or enable openmc.Settings.temperature handling "
                "(temperature_method / multipole interpolation) so intermediate "
                "temperatures are treated."
            ),
        ),
        (
            CROSS_SECTIONS,
            re.compile(
                r"cross section[s]? (file|library|data|xml)|"
                r"could not (find|read|open).*cross_section|"
                r"no cross sections (available|found)|"
                r"cross_sections\.xml",
                re.IGNORECASE,
            ),
            (
                "Cross-section data could not be located or parsed. Verify the "
                "provider 'cross_sections' path points at a valid cross_sections.xml "
                "and that the data directory is mounted into the container."
            ),
        ),
        (
            MPI_ABORT,
            re.compile(r"MPI_ABORT|mpi abort|orted: \w+ rank", re.IGNORECASE),
            (
                "The solver process crashed (MPI_ABORT). The root cause is usually the "
                "preceding engine error (e.g. nuclear-data or geometry failure); fix "
                "that first. If it occurs with no preceding error, the container may "
                "be out of memory or MPI may be misconfigured."
            ),
        ),
        (
            GEOMETRY,
            re.compile(
                r"geometry ?error|cannot find cell|no cell found|"
                r"particle (got )?lost|surface .*? not (found|defined)|"
                r"universe .*? not (found|defined)|geometry does not contain",
                re.IGNORECASE,
            ),
            (
                "The problem geometry is invalid or the source/particles fall outside "
                "it. Check geometry_config dimensions and the source_point location, "
                "and ensure every referenced material region is filled."
            ),
        ),
        (
            TALLY,
            re.compile(
                r"tally ?error|tally .*? (not|does not) (exist|have)|invalid filter",
                re.IGNORECASE,
            ),
            (
                "A tally or its filter is invalid. Check mesh_tallies / tally scores, "
                "filter IDs, and that referenced cells/meshes exist in the geometry."
            ),
        ),
    ],
    "festim": [
        (
            MESH_QUALITY,
            re.compile(
                r"mesh.*(too coarse|negative volume|invalid|quality)|"
                r"non-matching.*mesh|mesh generation failed",
                re.IGNORECASE,
            ),
            (
                "The FESTIM mesh is invalid or too coarse for the problem. "
                "Check mesh_config vertices/segments and ensure subdomain borders "
                "align with mesh nodes."
            ),
        ),
        (
            BC_SETUP,
            re.compile(
                r"boundary condition.*(not|does not|undefined)|unknown boundary|"
                r"bc.*(missing|undefined)|no bc",
                re.IGNORECASE,
            ),
            (
                "A boundary condition is missing, duplicated, or references an "
                "undefined surface/subdomain. Verify boundary_conditions entries "
                "and their species/subdomain references."
            ),
        ),
        (
            SOLVER_CONVERGENCE,
            re.compile(
                r"solver did not converge|newton.*(diverged|failed)|"
                r"nonlinear solver|divergence|solver reached maximum|"
                r"time step.*failed",
                re.IGNORECASE,
            ),
            (
                "The FESTIM solver did not converge. Increase max_iterations, "
                "relax atol/rtol, or refine the mesh in solver_config."
            ),
        ),
        (
            MATERIAL_PROPERTY,
            re.compile(
                r"material.*property|missing D_0|missing E_D|"
                r"solubility_law|pre-exponential|activation energy",
                re.IGNORECASE,
            ),
            (
                "A FESTIM material property is missing or invalid. Ensure every "
                "material defines extra.D_0, extra.E_D, and any solubility-related "
                "fields required by its boundary conditions."
            ),
        ),
    ],
    "generic": [
        (
            CONVERGENCE,
            re.compile(
                r"did not converge|maximum number of (iterations|resampling)|"
                r"failed to converge|stagnat",
                re.IGNORECASE,
            ),
            (
                "The solve did not converge. Increase batches/iterations or relax "
                "tolerances in solver_config; verify the model is well-posed."
            ),
        ),
        (
            INPUT_VALIDATION,
            re.compile(
                r"valueerror|keyerror|typeerror|validation ?error|"
                r"expected .*? (got|found)|missing .*? argument",
                re.IGNORECASE,
            ),
            (
                "The engine rejected the resolved configuration. This usually means a "
                "value in solver_config/geometry_config is out of range or mistyped."
            ),
        ),
        (
            ENVIRONMENT,
            re.compile(
                r"permission denied|no such file or directory|"
                r"cannot (write|create) .*? directory|"
                r"shared library|lib[a-z0-9]+\.so",
                re.IGNORECASE,
            ),
            (
                "An environment/container issue (file permissions, missing library, "
                "or unwritable output dir). Check the container mounts and that the "
                "output directory is writable."
            ),
        ),
    ],
}

# Flattened category -> hint lookup built from all signature sets. Engine-specific
# hints override generic hints for the same category.
_CATEGORY_HINTS: dict[str, str] = {}
for _signatures in _ERROR_SIGNATURES.values():
    for _cat, _pat, _hint in _signatures:
        _CATEGORY_HINTS[_cat] = _hint

# Default hint used when no category matched.
_DEFAULT_HINT = (
    "Inspect the engine log / run directory for the root cause; this is a "
    "run-time error raised by the provider, not a flowsheet configuration error."
)


# ---------------------------------------------------------------------------
# Structured run-error record
# ---------------------------------------------------------------------------


class ProviderRunError(BaseModel):
    """Structured record of a run-time failure that originated in a provider/engine.

    Attributes:
        category: One of the ``*`` constants above (``nuclear_data``, ``mpi_abort``,
            …). ``unknown`` when no signature matched.
        source: Always ``"provider"`` — distinguishes engine run-time failures
            from flowsheet/setup validation errors.
        message: A concise, human-readable one-line summary of the failure.
        type: The Python exception class name (e.g. ``RuntimeError``).
        detail: The full captured error text (engine stderr/stdout + traceback
            tail). May be long; this is what gets logged.
        hint: A concrete remediation suggestion where a category is recognised,
            else a generic pointer to check the engine log.
    """

    category: str = UNKNOWN
    source: str = "provider"
    message: str = ""
    type: str = "Exception"
    detail: str = ""
    hint: str = _DEFAULT_HINT

    @classmethod
    def from_exception(
        cls,
        exc: BaseException,
        captured: str = "",
        category: Optional[str] = None,
    ) -> "ProviderRunError":
        """Build a :class:`ProviderRunError` from an exception (category overridable)."""
        text = str(exc)
        if captured:
            text = f"{text}\n\n{captured}".strip()
        # Prefer an explicitly supplied category (e.g. a provider already knows).
        if category is None:
            category = _classify_text(text)
        return cls(
            category=category,
            message=_summarize(text) or f"{type(exc).__name__}: {text}",
            type=type(exc).__name__,
            detail=text,
            hint=_hint_for(category),
        )


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def classify_run_error(
    engine: str,
    exc: BaseException,
    captured: str = "",
) -> ProviderRunError:
    """Classify a provider run-time exception into a :class:`ProviderRunError`.

    Args:
        engine: Engine name (``"openmc"``, ``"festim"``, …). Used to select
            engine-specific error signatures; generic signatures are used as
            a fallback.
        exc: The exception raised by the engine run.
        captured: Optional captured engine stdout/stderr, used to improve the
            classification when the raised message is terse.

    Returns:
        A populated :class:`ProviderRunError` (``source="provider"``).
    """
    err = ProviderRunError.from_exception(
        exc, captured=captured, category=_classify_text(f"{exc}\n\n{captured}".strip(), engine)
    )
    err.message = f"[{engine}] {err.message}"
    return err


def make_failed_output(
    engine: str,
    sim_type: str,
    run_dir,
    err: ProviderRunError,
    unit: str = "",
) -> "EngineOutput":  # type: ignore[name-defined]  # imported lazily to avoid cycle
    """Build the ``EngineOutput(status="failed")`` providers return on a run error."""
    from processforge.types import EngineOutput

    return EngineOutput(
        status="failed",
        engine=engine,
        sim_type=sim_type,
        unit=unit,
        error=err,
        run_dir=str(getattr(run_dir, "resolve", lambda: run_dir)()),
        diagnostics={
            "run_dir": str(getattr(run_dir, "resolve", lambda: run_dir)()),
            "error": err.detail,
            "error_category": err.category,
            "error_source": err.source,
        },
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _classify_text(text: str, engine: str = "generic") -> str:
    """Return the first matching category for ``text``.

    Engine-specific signatures are checked first, then generic signatures.
    Unknown engines fall back to generic signatures only.
    """
    signatures: list[tuple[str, re.Pattern[str], str]] = []
    if engine in _ERROR_SIGNATURES and engine != "generic":
        signatures.extend(_ERROR_SIGNATURES[engine])
    signatures.extend(_ERROR_SIGNATURES["generic"])

    for category, pattern, _hint in signatures:
        if pattern.search(text):
            return category
    return UNKNOWN


def _hint_for(category: str) -> str:
    return _CATEGORY_HINTS.get(category, _DEFAULT_HINT)


def _summarize(text: str, limit: int = 280) -> str:
    """Return the most diagnostic single line of ``text`` (engine ERROR lines win)."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return ""
    for ln in lines:
        if ln.upper().startswith("ERROR") or "ERROR:" in ln.upper():
            return ln[:limit]
    # Fall back to the first line that looks like a message (not a banner/traceback).
    for ln in lines:
        if ln.startswith(("Traceback", "File ", "raise ", "Proc:", "NOTE:")):
            continue
        return ln[:limit]
    return lines[0][:limit]


__all__ = [
    # Exception hierarchy
    "ProviderError",
    "ProviderInitError",
    "ProviderNotAvailableError",
    "ProviderRuntimeError",
    "ProviderCleanupError",
    "ProviderValidationError",
    "ProviderConfigError",
    # Classification data model / helpers
    "ProviderRunError",
    "classify_run_error",
    "make_failed_output",
    # Categories
    "BC_SETUP",
    "CONVERGENCE",
    "CROSS_SECTIONS",
    "ENVIRONMENT",
    "GEOMETRY",
    "INPUT_VALIDATION",
    "MATERIAL_PROPERTY",
    "MESH_QUALITY",
    "MPI_ABORT",
    "NUCLEAR_DATA",
    "SOLVER_CONVERGENCE",
    "TALLY",
    "UNKNOWN",
]
