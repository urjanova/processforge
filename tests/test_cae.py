"""Tests for CAE case loading, evidence resolution and ``pf runs --cae``."""

import io
import json

import pytest
import typer
import yaml
from loguru import logger
from typer.testing import CliRunner

from processforge.cae import CAEError, FileCAERepository, generate_cae_reports
from processforge.cli.runs import runs
from processforge.result import save_results_zarr
from processforge.types import (
    EngineOutput,
    OutputArtifact,
    OutputField,
    Quantity,
    RunManifest,
)


@pytest.fixture
def logs():
    """Capture loguru error output (its default sink bypasses pytest capture)."""
    buffer = io.StringIO()
    sink_id = logger.add(buffer, level="ERROR")
    try:
        yield buffer
    finally:
        logger.remove(sink_id)


def _engine_output(tmp_path, *, keff=1.04, std=0.002) -> EngineOutput:
    statepoint = tmp_path / "statepoint.h5"
    statepoint.write_bytes(b"fake-statepoint")
    return EngineOutput(
        engine="openmc",
        sim_type="eigenvalue_reactor",
        status="completed",
        fields=[
            OutputField(
                name="k_eff",
                quantity=Quantity(value=keff, unit="", std_dev=std),
                kind="scalar",
                source="keff",
            )
        ],
        artifacts=[
            OutputArtifact(
                name="statepoint",
                kind="statepoint",
                local_path=str(statepoint),
                source="local",
            )
        ],
    )


def _build_run(tmp_path, monkeypatch, *, base="cae_flowsheet"):
    monkeypatch.setenv("PROCESSFORGE_OUTPUT_DIR", str(tmp_path))
    archive_path = tmp_path / f"{base}.pfarchive"
    runs_dir = archive_path / "runs"
    runs_dir.mkdir(parents=True)

    run_id = "20260101T000000Z_aaaaaa"
    output = _engine_output(tmp_path)
    stream_result = {"T": 300.0, "P": 101325.0, "flowrate": 1.0}
    save_results_zarr(
        {"openmc_solver": output, "feed": stream_result},
        str(archive_path / "results" / run_id / "results.zarr"),
        None,
    )
    manifest = RunManifest(
        run_id=run_id,
        timestamp="2026-01-01T00:00:00+00:00",
        mode="steady",
        units={"openmc_solver": output},
        streams={
            "feed": {
                "stream": "feed",
                "engine": "coolprop",
                "fields": [
                    OutputField(name="P", quantity=Quantity(value=101325.0, unit="Pa")),
                    OutputField(name="T", quantity=Quantity(value=300.0, unit="K")),
                ],
            }
        },
    )
    (runs_dir / f"{run_id}.json").write_text(manifest.model_dump_json())
    (archive_path / "latest_run").write_text(run_id)

    flowsheet_path = tmp_path / f"{base}.json"
    flowsheet_path.write_text(
        json.dumps(
            {
                "units": {
                    "openmc_solver": {"type": "SolverUnit", "provider": "openmc"}
                },
                "streams": {"feed": {"T": 300.0, "P": 101325.0, "flowrate": 1.0}},
            }
        )
    )
    return flowsheet_path, run_id


def _write_repo(root, *, cases: dict[str, dict], links: list[dict]):
    cases_dir = root / "cases"
    cases_dir.mkdir(parents=True)
    for filename, case in cases.items():
        (cases_dir / filename).write_text(yaml.safe_dump(case, sort_keys=False))
    index = {"schema_version": "1.0", "flowsheets": links}
    (root / "index.yaml").write_text(yaml.safe_dump(index, sort_keys=False))


def _case(*, case_id="sample-cae", version="1.0.0", claims=None, args=None, evidence=None):
    return {
        "schema_version": "1.0",
        "case": {
            "id": case_id,
            "title": "Sample case",
            "version": version,
            "jurisdictions": ["UK-GDA"],
        },
        "requirements": [
            {"id": "REQ-001", "statement": "The system shall be safe."},
        ],
        "claims": claims
        or [
            {
                "id": "CLAIM-001",
                "statement": "The reactor is subcritical.",
                "jurisdiction": "UK-GDA",
                "requirement_refs": ["REQ-001"],
                "argument_refs": ["ARG-001"],
            }
        ],
        "arguments": args
        or [
            {
                "id": "ARG-001",
                "statement": "k_eff is below unity.",
                "evidence_refs": ["EVD-001"],
            }
        ],
        "evidence": evidence
        or [
            {
                "id": "EVD-001",
                "description": "k_eff estimate",
                "ref": {
                    "type": "unit",
                    "name": "openmc_solver",
                    "field": "k_eff",
                    "extract": "value",
                },
            }
        ],
    }


def _invoke(flowsheet_path, *args):
    app = typer.Typer()
    app.command()(runs)
    return CliRunner().invoke(app, [str(flowsheet_path), *args])


def test_cae_report_resolves_unit_and_artifact_evidence(tmp_path, monkeypatch):
    flowsheet_path, run_id = _build_run(tmp_path, monkeypatch)
    repo = tmp_path / "cae"
    _write_repo(
        repo,
        cases={"sample-cae.ca.yaml": _case()},
        links=[{"path": str(flowsheet_path), "cases": [{"case_id": "sample-cae"}]}],
    )

    result = _invoke(flowsheet_path, run_id, "--cae", "--cae-dir", str(repo))

    assert result.exit_code == 0, result.output
    assert "CAE ASSURANCE REPORT" in result.output
    assert "CLAIM-001" in result.output
    assert "k_eff estimate" in result.output
    assert "1.04" in result.output
    assert run_id in result.output
    assert "Flowsheet hash" in result.output


def test_cae_report_artifact_evidence(tmp_path, monkeypatch):
    flowsheet_path, run_id = _build_run(tmp_path, monkeypatch)
    repo = tmp_path / "cae"
    evidence = [
        {
            "id": "EVD-001",
            "description": "statepoint retained",
            "ref": {"type": "artifact", "name": "statepoint", "unit": "openmc_solver"},
        }
    ]
    case = _case(evidence=evidence)
    case["arguments"][0]["evidence_refs"] = ["EVD-001"]
    _write_repo(
        repo,
        cases={"sample-cae.ca.yaml": case},
        links=[{"path": str(flowsheet_path), "cases": [{"case_id": "sample-cae"}]}],
    )

    result = _invoke(flowsheet_path, run_id, "--cae", "--cae-dir", str(repo))
    assert result.exit_code == 0, result.output
    assert "artifact 'statepoint'" in result.output
    assert "statepoint.h5" in result.output
    assert "resolved" in result.output


def test_cae_stream_evidence_supports_extract_ops(tmp_path, monkeypatch):
    flowsheet_path, run_id = _build_run(tmp_path, monkeypatch)
    repo = tmp_path / "cae"
    evidence = [
        {
            "id": "EVD-001",
            "description": "mean feed pressure",
            "ref": {"type": "stream", "name": "feed", "field": "P", "extract": "mean"},
        }
    ]
    case = _case(evidence=evidence)
    case["arguments"][0]["evidence_refs"] = ["EVD-001"]
    _write_repo(
        repo,
        cases={"sample-cae.ca.yaml": case},
        links=[{"path": str(flowsheet_path), "cases": [{"case_id": "sample-cae"}]}],
    )

    result = _invoke(flowsheet_path, run_id, "--cae", "--cae-dir", str(repo))
    assert result.exit_code == 0, result.output
    assert "stream 'feed' field 'P' (mean)" in result.output
    assert "101325" in result.output


def test_cae_fails_when_flowsheet_has_no_link(tmp_path, monkeypatch, logs):
    flowsheet_path, run_id = _build_run(tmp_path, monkeypatch)
    repo = tmp_path / "cae"
    _write_repo(repo, cases={}, links=[])

    result = _invoke(flowsheet_path, run_id, "--cae", "--cae-dir", str(repo))
    assert result.exit_code == 1
    assert "No CAE case linked" in logs.getvalue()


def test_cae_fails_on_unresolved_evidence(tmp_path, monkeypatch, logs):
    flowsheet_path, run_id = _build_run(tmp_path, monkeypatch)
    repo = tmp_path / "cae"
    evidence = [
        {
            "id": "EVD-001",
            "description": "missing field",
            "ref": {
                "type": "unit",
                "name": "openmc_solver",
                "field": "does_not_exist",
                "extract": "value",
            },
        }
    ]
    case = _case(evidence=evidence)
    _write_repo(
        repo,
        cases={"sample-cae.ca.yaml": case},
        links=[{"path": str(flowsheet_path), "cases": [{"case_id": "sample-cae"}]}],
    )

    result = _invoke(flowsheet_path, run_id, "--cae", "--cae-dir", str(repo))
    assert result.exit_code == 1
    err = logs.getvalue()
    assert "Cannot resolve CAE evidence" in err
    assert "does_not_exist" in err


def test_cae_fails_when_unit_not_in_flowsheet(tmp_path, monkeypatch, logs):
    flowsheet_path, run_id = _build_run(tmp_path, monkeypatch)
    repo = tmp_path / "cae"
    evidence = [
        {
            "id": "EVD-001",
            "description": "wrong unit",
            "ref": {
                "type": "unit",
                "name": "not_in_flowsheet",
                "field": "k_eff",
                "extract": "value",
            },
        }
    ]
    case = _case(evidence=evidence)
    _write_repo(
        repo,
        cases={"sample-cae.ca.yaml": case},
        links=[{"path": str(flowsheet_path), "cases": [{"case_id": "sample-cae"}]}],
    )

    result = _invoke(flowsheet_path, run_id, "--cae", "--cae-dir", str(repo))
    assert result.exit_code == 1
    assert "not declared in the flowsheet" in logs.getvalue()


def test_cae_jurisdiction_filter(tmp_path, monkeypatch):
    flowsheet_path, run_id = _build_run(tmp_path, monkeypatch)
    repo = tmp_path / "cae"
    claims = [
        {
            "id": "CLAIM-UK",
            "statement": "UK claim",
            "jurisdiction": "UK-GDA",
            "requirement_refs": ["REQ-001"],
            "argument_refs": ["ARG-001"],
        },
        {
            "id": "CLAIM-US",
            "statement": "US claim",
            "jurisdiction": "US-NRC",
            "requirement_refs": ["REQ-001"],
            "argument_refs": ["ARG-US"],
        },
    ]
    args = [
        {"id": "ARG-001", "statement": "uk", "evidence_refs": ["EVD-001"]},
        {"id": "ARG-US", "statement": "us", "evidence_refs": ["EVD-002"]},
    ]
    evidence = [
        {
            "id": "EVD-001",
            "description": "uk evidence",
            "ref": {"type": "unit", "name": "openmc_solver", "field": "k_eff"},
        },
        {
            "id": "EVD-002",
            "description": "us evidence",
            "ref": {"type": "unit", "name": "openmc_solver", "field": "k_eff"},
        },
    ]
    _write_repo(
        repo,
        cases={"sample-cae.ca.yaml": _case(claims=claims, args=args, evidence=evidence)},
        links=[{"path": str(flowsheet_path), "cases": [{"case_id": "sample-cae"}]}],
    )

    result = _invoke(
        flowsheet_path, run_id, "--cae", "--cae-dir", str(repo), "--jurisdiction", "US-NRC"
    )
    assert result.exit_code == 0, result.output
    assert "CLAIM-US" in result.output
    assert "CLAIM-UK" not in result.output
    assert "uk evidence" not in result.output
    assert "Jurisdiction    : US-NRC" in result.output


def test_cae_fails_for_unknown_jurisdiction(tmp_path, monkeypatch, logs):
    flowsheet_path, run_id = _build_run(tmp_path, monkeypatch)
    repo = tmp_path / "cae"
    _write_repo(
        repo,
        cases={"sample-cae.ca.yaml": _case()},
        links=[{"path": str(flowsheet_path), "cases": [{"case_id": "sample-cae"}]}],
    )

    result = _invoke(
        flowsheet_path, run_id, "--cae", "--cae-dir", str(repo), "--jurisdiction", "EU"
    )
    assert result.exit_code == 1
    assert "No CAE claims found for jurisdiction 'EU'" in logs.getvalue()


def test_cae_requires_run_id(tmp_path, monkeypatch, logs):
    flowsheet_path, _ = _build_run(tmp_path, monkeypatch)
    result = _invoke(flowsheet_path, "--cae")
    assert result.exit_code == 1
    assert "requires a run id" in logs.getvalue()


def test_repository_selects_highest_version(tmp_path):
    repo_root = tmp_path / "cae"
    _write_repo(
        repo_root,
        cases={
            "sample-cae-1.0.0.ca.yaml": _case(version="1.0.0"),
            "sample-cae-1.2.0.ca.yaml": _case(version="1.2.0"),
        },
        links=[],
    )
    repo = FileCAERepository(str(repo_root))
    assert repo.load_case("sample-cae").case.version == "1.2.0"
    assert repo.load_case("sample-cae", "1.0.0").case.version == "1.0.0"
    with_error = None
    try:
        repo.load_case("sample-cae", "9.9.9")
    except CAEError as exc:
        with_error = str(exc)
    assert with_error and "no version '9.9.9'" in with_error


def test_case_validation_rejects_unknown_evidence(tmp_path):
    repo_root = tmp_path / "cae"
    bad = _case()
    bad["arguments"][0]["evidence_refs"] = ["EVD-MISSING"]
    _write_repo(
        repo_root,
        cases={"sample-cae.ca.yaml": bad},
        links=[],
    )
    repo = FileCAERepository(str(repo_root))
    try:
        repo.list_cases()
    except CAEError as exc:
        assert "unknown evidence" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("expected CAEError")


def test_generate_cae_reports_without_repo_raises(tmp_path, monkeypatch):
    flowsheet = tmp_path / "x.json"
    flowsheet.write_text("{}")
    monkeypatch.delenv("PROCESSFORGE_CAE_DIR", raising=False)
    monkeypatch.chdir(tmp_path)
    try:
        generate_cae_reports(
            flowsheet_path=str(flowsheet),
            flowsheet_config={},
            run_id="r",
            run_summary={"present": True},
            manifest=None,
            zarr_path="",
            cae_root=str(tmp_path / "missing"),
        )
    except CAEError as exc:
        assert "No CAE repository found" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("expected CAEError")
