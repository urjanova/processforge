"""Tests for the provider run-error classification utility."""
from __future__ import annotations

import pathlib

from processforge.providers.errors import (
    BC_SETUP,
    CROSS_SECTIONS,
    ENVIRONMENT,
    GEOMETRY,
    INPUT_VALIDATION,
    MATERIAL_PROPERTY,
    MESH_QUALITY,
    MPI_ABORT,
    NUCLEAR_DATA,
    ProviderCleanupError,
    ProviderConfigError,
    ProviderError,
    ProviderInitError,
    ProviderNotAvailableError,
    ProviderRunError,
    ProviderRuntimeError,
    ProviderValidationError,
    SOLVER_CONVERGENCE,
    UNKNOWN,
    classify_run_error,
    make_failed_output,
)
from processforge.types import EngineOutput


def _classify(msg: str, engine: str = "openmc", captured: str = "") -> ProviderRunError:
    return classify_run_error(engine, RuntimeError(msg), captured=captured)


# ---------------------------------------------------------------------------
# Exception hierarchy
# ---------------------------------------------------------------------------


def test_provider_error_hierarchy():
    assert issubclass(ProviderNotAvailableError, ProviderInitError)
    assert issubclass(ProviderInitError, ProviderError)
    assert issubclass(ProviderRuntimeError, ProviderError)
    assert issubclass(ProviderCleanupError, ProviderError)
    assert issubclass(ProviderValidationError, ProviderError)
    assert issubclass(ProviderConfigError, ProviderError)


def test_provider_not_available_error_carries_message():
    exc = ProviderNotAvailableError("backend missing")
    assert str(exc) == "backend missing"
    assert isinstance(exc, ProviderInitError)


# ---------------------------------------------------------------------------
# OpenMC categories
# ---------------------------------------------------------------------------


def test_nuclear_data_category_and_hint():
    err = _classify(
        "Nuclear data library does not contain cross sections for Ni58 at or near "
        "400.000000 K. Available temperatures are 300 K."
    )
    assert err.category == NUCLEAR_DATA
    assert err.source == "provider"
    assert "temperature" in err.hint.lower()
    assert err.type == "RuntimeError"


def test_mpi_abort_category():
    err = _classify(
        "RuntimeError: something bad\n"
        "MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD"
    )
    assert err.category == MPI_ABORT


def test_cross_sections_category():
    err = _classify("Could not find cross_sections.xml at the configured path")
    assert err.category == CROSS_SECTIONS


def test_geometry_category():
    err = _classify("GeometryError: Could not find cell containing particle")
    assert err.category == GEOMETRY


# ---------------------------------------------------------------------------
# FESTIM categories
# ---------------------------------------------------------------------------


def test_festim_mesh_quality_category():
    err = _classify(
        "Mesh generation failed: negative volume cells detected",
        engine="festim",
    )
    assert err.category == MESH_QUALITY
    assert "mesh" in err.hint.lower()


def test_festim_bc_setup_category():
    err = _classify(
        "Boundary condition undefined for surface subdomain id=2",
        engine="festim",
    )
    assert err.category == BC_SETUP


def test_festim_solver_convergence_category():
    err = _classify(
        "Newton solver diverged after 25 iterations",
        engine="festim",
    )
    assert err.category == SOLVER_CONVERGENCE


def test_festim_material_property_category():
    err = _classify(
        "Material 'tungsten' is missing property pre-exponential D_0",
        engine="festim",
    )
    assert err.category == MATERIAL_PROPERTY


# ---------------------------------------------------------------------------
# Engine scoping / fallback behavior
# ---------------------------------------------------------------------------


def test_unknown_fallback():
    err = _classify("some totally opaque engine complaint with no known signature")
    assert err.category == UNKNOWN
    assert err.source == "provider"
    assert "run-time" in err.hint


def test_message_prefixed_with_engine():
    err = _classify("Nuclear data library does not contain cross sections for Ni58")
    assert err.message.startswith("[openmc]")


def test_captured_text_improves_classification():
    # Raised message is terse; the real signal is in the captured engine stderr.
    err = _classify(
        "RuntimeError",
        captured="ERROR: Nuclear data library does not contain cross sections for U235",
    )
    assert err.category == NUCLEAR_DATA


def test_openmc_specific_signature_does_not_match_festim_engine():
    # The same nuclear-data text should be UNKNOWN for festim because it has no
    # cross-section concept; generic signatures should still apply.
    err = _classify(
        "Nuclear data library does not contain cross sections for Ni58",
        engine="festim",
    )
    assert err.category == UNKNOWN


def test_generic_environment_category_matches_any_engine():
    err = _classify(
        "Permission denied when writing to outputs/festim",
        engine="festim",
    )
    assert err.category == ENVIRONMENT


def test_generic_input_validation_category_matches_any_engine():
    err = _classify(
        "ValueError: solver_config expects float got str",
        engine="my_custom_engine",
    )
    assert err.category == INPUT_VALIDATION


# ---------------------------------------------------------------------------
# ProviderRunError record / output builders
# ---------------------------------------------------------------------------


def test_from_exception_sets_fields():
    err = ProviderRunError.from_exception(ValueError("expected float got str"))
    assert err.type == "ValueError"
    assert err.category == INPUT_VALIDATION
    assert err.detail


def test_make_failed_output_shape(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    err = _classify("Nuclear data library does not contain cross sections for Ni58")
    out = make_failed_output("openmc", "eigenvalue_reactor", run_dir, err, unit="msre")
    assert isinstance(out, EngineOutput)
    assert out.status == "failed"
    assert out.engine == "openmc"
    assert out.sim_type == "eigenvalue_reactor"
    assert out.unit == "msre"
    assert out.error is err
    assert out.diagnostics["error_category"] == NUCLEAR_DATA
    assert out.diagnostics["error_source"] == "provider"
    assert out.diagnostics["run_dir"] == str(run_dir.resolve())
    assert out.diagnostics["error"]


def test_engine_output_roundtrip_serialization(tmp_path):
    err = _classify("Nuclear data library does not contain cross sections for Ni58")
    out = make_failed_output("openmc", "eigenvalue_reactor", tmp_path, err)
    dumped = out.model_dump()
    restored = EngineOutput.model_validate(dumped)
    assert restored.status == "failed"
    assert restored.error is not None
    assert restored.error.category == NUCLEAR_DATA
