"""Tests for the CLI layer — simulate.py, cli/common.py, cli/display.py, and command modules."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from loguru import logger


TESTS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = TESTS_DIR.parent


# ---------------------------------------------------------------------------
# Loguru capture fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def log_output():
    """Capture loguru output into a list of strings."""
    messages: list[str] = []
    handler_id = logger.add(lambda m: messages.append(str(m)), format="{message}")
    yield messages
    logger.remove(handler_id)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def coolprop_flowsheet(tmp_path):
    """Write a minimal valid flowsheet with only coolprop."""
    flowsheet = {
        "providers": {"coolprop": {"type": "coolprop"}},
        "streams": {},
        "units": {},
        "simulation": {"mode": "steady"},
    }
    path = tmp_path / "test_flowsheet.json"
    with open(path, "w") as f:
        json.dump(flowsheet, f)
    return path


@pytest.fixture
def empty_flowsheet(tmp_path):
    """Write a flowsheet with no providers."""
    flowsheet = {
        "providers": {},
        "streams": {},
        "units": {},
        "simulation": {"mode": "steady"},
    }
    path = tmp_path / "empty.json"
    with open(path, "w") as f:
        json.dump(flowsheet, f)
    return path


@pytest.fixture
def invalid_json(tmp_path):
    """Write an invalid JSON file."""
    path = tmp_path / "bad.json"
    path.write_text("{invalid json content")
    return path


# ---------------------------------------------------------------------------
# main() / argparse dispatch
# ---------------------------------------------------------------------------

class TestMain:
    def test_no_command_prints_help_and_exits(self):
        from processforge.simulate import main

        with patch("sys.argv", ["processforge"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 1

    def test_help_exits_zero(self):
        from processforge.simulate import main

        with patch("sys.argv", ["processforge", "--help"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0

    def test_unknown_command_exits_nonzero(self):
        from processforge.simulate import main

        with patch("sys.argv", ["processforge", "nonexistent-cmd"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code != 0

    def test_debug_flag_on_unhandled_exception(self, tmp_path):
        from processforge.simulate import main

        bad_path = tmp_path / "no_such_file.json"
        with patch("sys.argv", ["processforge", "--debug", "validate", str(bad_path)]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 1


# ---------------------------------------------------------------------------
# cli/common.py — output_root
# ---------------------------------------------------------------------------

class TestOutputRoot:
    def test_default(self):
        from processforge.cli.common import output_root

        with patch.dict(os.environ, {}, clear=True):
            assert output_root() == "outputs"

    def test_env_override(self):
        from processforge.cli.common import output_root

        with patch.dict(os.environ, {"PROCESSFORGE_OUTPUT_DIR": "/data"}):
            assert output_root() == "/data"


# ---------------------------------------------------------------------------
# cli/common.py — require_existing_file
# ---------------------------------------------------------------------------

class TestRequireExistingFile:
    def test_missing_file_exits(self):
        from processforge.cli.common import require_existing_file

        with pytest.raises(SystemExit) as exc_info:
            require_existing_file("/no/such/file.json")
        assert exc_info.value.code == 1

    def test_existing_file_ok(self, tmp_path):
        from processforge.cli.common import require_existing_file

        f = tmp_path / "exists.json"
        f.write_text("{}")
        require_existing_file(str(f))


# ---------------------------------------------------------------------------
# cli/common.py — extract_providers
# ---------------------------------------------------------------------------

class TestExtractProviders:
    def test_valid_file(self, coolprop_flowsheet):
        from processforge.cli.common import extract_providers

        providers = extract_providers(str(coolprop_flowsheet))
        assert "coolprop" in providers

    def test_empty_providers_warns(self, empty_flowsheet, log_output):
        from processforge.cli.common import extract_providers

        providers = extract_providers(str(empty_flowsheet))
        assert providers == {}
        assert any("No providers declared" in m for m in log_output)

    def test_invalid_json_exits(self, invalid_json):
        from processforge.cli.common import extract_providers

        with pytest.raises(SystemExit) as exc_info:
            extract_providers(str(invalid_json))
        assert exc_info.value.code == 1

    def test_missing_file_exits(self):
        from processforge.cli.common import extract_providers

        with pytest.raises(SystemExit) as exc_info:
            extract_providers("/no/such/file.json")
        assert exc_info.value.code == 1


# ---------------------------------------------------------------------------
# cli/common.py — build_run_metadata
# ---------------------------------------------------------------------------

class TestBuildRunMetadata:
    def test_basic(self):
        from processforge.cli.common import build_run_metadata

        meta = build_run_metadata({"streams": {}}, 1e-6, 50, "scipy")
        assert meta["version"]
        assert meta["flowsheet_hash"]
        assert meta["solver_settings"]["tol"] == 1e-6
        assert meta["solver_settings"]["max_iter"] == 50
        assert meta["solver_settings"]["backend"] == "scipy"

    def test_deterministic_hash(self):
        from processforge.cli.common import build_run_metadata

        m1 = build_run_metadata({"a": 1}, 1e-6, 50, "scipy")
        m2 = build_run_metadata({"a": 1}, 1e-6, 50, "scipy")
        assert m1["flowsheet_hash"] == m2["flowsheet_hash"]

    def test_different_configs_different_hash(self):
        from processforge.cli.common import build_run_metadata

        m1 = build_run_metadata({"a": 1}, 1e-6, 50, "scipy")
        m2 = build_run_metadata({"a": 2}, 1e-6, 50, "scipy")
        assert m1["flowsheet_hash"] != m2["flowsheet_hash"]


# ---------------------------------------------------------------------------
# cli/common.py — build_divergence_report
# ---------------------------------------------------------------------------

class TestBuildDivergenceReport:
    def test_basic(self):
        from processforge.cli.common import build_divergence_report

        report = build_divergence_report(
            drifted_params=["streams.feed.T"],
            solver_stats={"final_norm": 1e-3, "iterations": 10},
            x_last=[1.0, 2.0, 3.0],
            var_names=["a", "b", "c"],
            breakdown=[],
        )
        assert report["drifted_params"] == ["streams.feed.T"]
        assert report["final_norm"] == 1e-3
        assert report["x_last"] == [1.0, 2.0, 3.0]
        assert report["var_names"] == ["a", "b", "c"]
        assert report["timestamp"]

    def test_empty_x_last_and_var_names(self):
        from processforge.cli.common import build_divergence_report

        report = build_divergence_report(
            drifted_params=[],
            solver_stats={},
            x_last=[],
            var_names=[],
            breakdown=[],
        )
        assert report["x_last"] == []
        assert report["var_names"] == []

    def test_numpy_array_x_last(self):
        import numpy as np
        from processforge.cli.common import build_divergence_report

        report = build_divergence_report(
            drifted_params=[],
            solver_stats={},
            x_last=np.array([1.0, 2.0]),
            var_names=["a", "b"],
            breakdown=[],
        )
        assert report["x_last"] == [1.0, 2.0]


# ---------------------------------------------------------------------------
# cli/display.py
# ---------------------------------------------------------------------------

class TestDisplay:
    def test_print_structural_diff_empty(self, log_output):
        from processforge.cli.display import print_structural_diff

        print_structural_diff({"added": {}, "modified": {}, "removed": {}})
        assert any("no structural changes" in m for m in log_output)

    def test_print_structural_diff_with_changes(self, log_output):
        from processforge.cli.display import print_structural_diff

        diff = {
            "added": {"mixer": "Mixer"},
            "modified": {},
            "removed": {"old_heater": "Heater"},
        }
        print_structural_diff(diff)
        assert any("+ mixer" in m for m in log_output)
        assert any("- old_heater" in m for m in log_output)

    def test_print_unit_mismatches_empty(self, log_output):
        from processforge.cli.display import print_unit_mismatches

        print_unit_mismatches([])
        assert any("No unit mismatches" in m for m in log_output)


# ---------------------------------------------------------------------------
# cli/init.py — scaffold-only mode
# ---------------------------------------------------------------------------

class TestCmdInit:
    def test_scaffold_only(self, tmp_path, log_output):
        from processforge.cli.init import cmd_init

        args = argparse.Namespace(flowsheet=None, path=str(tmp_path))
        cmd_init(args)
        assert (tmp_path / ".processforge" / "config.json").exists()
        assert (tmp_path / "outputs").is_dir()
        assert any("initialised successfully" in m for m in log_output)

    def test_stale_snapshots_removed(self, tmp_path, log_output):
        from processforge.cli.init import cmd_init

        outputs = tmp_path / "outputs"
        outputs.mkdir()
        stale = outputs / "old.pfstate"
        stale.mkdir()
        (stale / "data.bin").write_bytes(b"stale")

        args = argparse.Namespace(flowsheet=None, path=str(tmp_path))
        cmd_init(args)
        assert not stale.exists()
        assert any("1 stale snapshot" in m for m in log_output)

    def test_missing_flowsheet_exits(self, tmp_path):
        from processforge.cli.init import cmd_init

        args = argparse.Namespace(flowsheet="nonexistent.json", path=str(tmp_path))
        with pytest.raises(SystemExit) as exc_info:
            cmd_init(args)
        assert exc_info.value.code == 1

    def test_empty_providers_flowsheet(self, tmp_path, log_output):
        from processforge.cli.init import cmd_init

        fs = tmp_path / "empty.json"
        fs.write_text(json.dumps({"providers": {}, "streams": {}, "units": {}}))
        args = argparse.Namespace(flowsheet=str(fs), path=str(tmp_path))
        cmd_init(args)
        assert any("No providers declared" in m for m in log_output)
        assert any("initialised successfully" in m for m in log_output)


# ---------------------------------------------------------------------------
# cli/validate.py
# ---------------------------------------------------------------------------

class TestCmdValidate:
    def test_missing_file_exits(self, tmp_path):
        from processforge.cli.validate import cmd_validate

        args = argparse.Namespace(flowsheet=str(tmp_path / "nope.json"))
        with pytest.raises(SystemExit) as exc_info:
            cmd_validate(args)
        assert exc_info.value.code == 1

    def test_valid_coolprop_flowsheet(self, coolprop_flowsheet, log_output):
        from processforge.cli.validate import cmd_validate

        args = argparse.Namespace(flowsheet=str(coolprop_flowsheet))
        cmd_validate(args)
        assert any("All providers ready" in m for m in log_output)


# ---------------------------------------------------------------------------
# cli/run.py
# ---------------------------------------------------------------------------

class TestCmdRun:
    def test_missing_file_exits(self, tmp_path):
        from processforge.cli.run import cmd_run

        args = argparse.Namespace(flowsheet=str(tmp_path / "nope.json"), export_images=False)
        with pytest.raises(SystemExit) as exc_info:
            cmd_run(args)
        assert exc_info.value.code == 1


# ---------------------------------------------------------------------------
# cli/diagram.py
# ---------------------------------------------------------------------------

class TestCmdDiagram:
    def test_missing_file_exits(self, tmp_path):
        from processforge.cli.diagram import cmd_diagram

        args = argparse.Namespace(
            flowsheet=str(tmp_path / "nope.json"),
            output_dir="diagrams",
            format="png",
        )
        with pytest.raises(SystemExit) as exc_info:
            cmd_diagram(args)
        assert exc_info.value.code == 1

    def test_invalid_json_exits(self, invalid_json):
        from processforge.cli.diagram import cmd_diagram

        args = argparse.Namespace(
            flowsheet=str(invalid_json),
            output_dir="diagrams",
            format="png",
        )
        with pytest.raises(SystemExit) as exc_info:
            cmd_diagram(args)
        assert exc_info.value.code == 1
