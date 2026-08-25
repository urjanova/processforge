"""Provider-side run-error classification.

When a simulation is executed *inside* a provider/engine (e.g. OpenMC or FESTIM
running in a Docker container) the run can fail for reasons that live entirely
in the engine — missing nuclear-data temperatures, an MPI process abort, a
geometry that excludes a source point, a tally that references an unknown
filter, … . These are distinct from *flowsheet/setup* errors (bad JSON, unknown
material, schema validation), which are caught earlier during initialization or
flowsheet validation.

This module gives every provider a single, structured way to capture such a
run-time failure:

* :class:`ProviderRunError` — a typed, engine-agnostic record of what failed,
  attributed to ``source="provider"`` so Processforge (and the user) can tell a
  runtime engine error apart from a flowsheet configuration error.
* :func:`classify_run_error` — turn an exception (plus any captured
  stdout/stderr from the engine) into a :class:`ProviderRunError`, picking a
  category and a
  concrete remediation hint where one is known.
* :func:`make_failed_output` — build the ``EngineOutput(status="failed")`` that
  providers return, populated with the classification.

Categories are intentionally broad so a new engine signature can be added here
without touching the providers.
"""
from __future__ import annotations

import re
from typing import Optional

from pydantic import BaseModel

# Error categories. Strings (not an Enum) so they serialize cleanly in JSON
# through the provider HTTP API and are easy to match on in tests/docs.
NUCLEAR_DATA = "nuclear_data"
CROSS_SECTIONS = "cross_sections"
MPI_ABORT = "mpi_abort"
GEOMETRY = "geometry"
TALLY = "tally"
CONVERGENCE = "convergence"
INPUT_VALIDATION = "input_validation"
ENVIRONMENT = "environment"
UNKNOWN = "unknown"

# (category, compiled-regex, remediation hint). First match wins.
# Order matters: more specific (and more actionable) patterns come first.
_ERROR_SIGNATURES: list[tuple[str, re.Pattern[str], str]] = [
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
        re.compile(r"tally ?error|tally .*? (not|does not) (exist|have)|invalid filter", re.IGNORECASE),
        (
            "A tally or its filter is invalid. Check mesh_tallies / tally scores, "
            "filter IDs, and that referenced cells/meshes exist in the geometry."
        ),
    ),
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
]


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
    hint: str = (
        "Inspect the engine log / run directory for the root cause; this is a "
        "run-time error raised by the provider, not a flowsheet configuration error."
    )

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


def classify_run_error(
    engine: str,
    exc: BaseException,
    captured: str = "",
) -> ProviderRunError:
    """Classify a provider run-time exception into a :class:`ProviderRunError`.

    Args:
        engine: Engine name (``"openmc"``, ``"festim"``, …) — recorded for context
            only (it does not change classification today).
        exc: The exception raised by the engine run.
        captured: Optional captured engine stdout/stderr, used to improve the
            classification when the raised message is terse.

    Returns:
        A populated :class:`ProviderRunError` (``source="provider"``).
    """
    err = ProviderRunError.from_exception(exc, captured=captured)
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
        diagnostics={
            "run_dir": str(getattr(run_dir, "resolve", lambda: run_dir)()),
            "error": err.detail,
            "error_category": err.category,
            "error_source": err.source,
        },
    )


def _classify_text(text: str) -> str:
    for category, pattern, _hint in _ERROR_SIGNATURES:
        if pattern.search(text):
            return category
    return UNKNOWN


def _hint_for(category: str) -> str:
    for cat, _pattern, hint in _ERROR_SIGNATURES:
        if cat == category:
            return hint
    return ProviderRunError().hint


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
