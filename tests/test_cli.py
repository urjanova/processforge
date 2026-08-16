"""Tests for the CLI layer — simulate.py, cli/common.py, cli/display.py, and command modules."""
from __future__ import annotations

import json
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner
from loguru import logger

from processforge.simulate import app


TESTS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = TESTS_DIR.parent

runner = CliRunner()


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
# main() / typer dispatch
# ---------------------------------------------------------------------------

class TestMain:
    def test_no_command_prints_help_and_exits(self):
        result = runner.invoke(app, [])
        assert result.exit_code == 1

    def test_help_exits_zero(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "Processforge" in result.output

    def test_unknown_command_exits_nonzero(self):
        result = runner.invoke(app, ["nonexistent-cmd"])
        assert result.exit_code != 0

    def test_debug_flag_on_unhandled_exception(self, tmp_path):
        bad_path = tmp_path / "no_such_file.json"
        result = runner.invoke(app, ["--debug", "validate", str(bad_path)])
        assert result.exit_code == 1


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
        from processforge.cli.init import init

        init(flowsheet=None, path=str(tmp_path))
        assert (tmp_path / ".processforge" / "config.json").exists()
        assert (tmp_path / "outputs").is_dir()
        assert any("initialised successfully" in m for m in log_output)

    def test_stale_snapshots_removed(self, tmp_path, log_output):
        from processforge.cli.init import init

        outputs = tmp_path / "outputs"
        outputs.mkdir()
        stale = outputs / "old.pfstate"
        stale.mkdir()
        (stale / "data.bin").write_bytes(b"stale")

        init(flowsheet=None, path=str(tmp_path))
        assert not stale.exists()
        assert any("1 stale snapshot" in m for m in log_output)

    def test_missing_flowsheet_exits(self, tmp_path):
        from processforge.cli.init import init

        with pytest.raises(SystemExit) as exc_info:
            init(flowsheet="nonexistent.json", path=str(tmp_path))
        assert exc_info.value.code == 1

    def test_empty_providers_flowsheet(self, tmp_path, log_output):
        from processforge.cli.init import init

        fs = tmp_path / "empty.json"
        fs.write_text(json.dumps({"providers": {}, "streams": {}, "units": {}}))
        init(flowsheet=str(fs), path=str(tmp_path))
        assert any("No providers declared" in m for m in log_output)
        assert any("initialised successfully" in m for m in log_output)


# ---------------------------------------------------------------------------
# cli/validate.py
# ---------------------------------------------------------------------------

class TestCmdValidate:
    def test_missing_file_exits(self, tmp_path):
        from processforge.cli.validate import validate

        with pytest.raises(SystemExit) as exc_info:
            validate(flowsheet=str(tmp_path / "nope.json"))
        assert exc_info.value.code == 1

    def test_valid_coolprop_flowsheet(self, coolprop_flowsheet, log_output):
        from processforge.cli.validate import validate

        validate(flowsheet=str(coolprop_flowsheet))
        assert any("All providers ready" in m for m in log_output)


# ---------------------------------------------------------------------------
# cli/run.py
# ---------------------------------------------------------------------------

class TestCmdRun:
    def test_missing_file_exits(self, tmp_path):
        from processforge.cli.run import run

        with pytest.raises(SystemExit) as exc_info:
            run(flowsheet=str(tmp_path / "nope.json"), export_images=False)
        assert exc_info.value.code == 1

    def _containerized_config(self, flowsheet_path):
        return {
            "providers": {"openmc": {"type": "openmc"}},
            "streams": {},
            "units": {},
            "simulation": {"mode": "steady"},
            "_config_path": str(flowsheet_path),
        }

    def _patch_run_os(self, run_mod):
        """Patch run's ``os`` so the compose-file existence check passes
        while keeping real ``makedirs`` behaviour."""
        real_os = run_mod.os
        mock_os = MagicMock()
        mock_os.path.exists.return_value = True
        mock_os.makedirs = real_os.makedirs
        mock_os.path.join = real_os.path.join
        return mock_os

    def test_containerized_flowsheet_runs_without_lifecycle(self, tmp_path, log_output):
        import processforge.cli.run as run_mod
        from processforge.cli.run import run

        flowsheet = tmp_path / "fs.json"
        flowsheet.write_text(json.dumps(self._containerized_config(flowsheet)["providers"]))

        mock_fs = MagicMock()
        mock_fs.run.return_value = {}
        mock_fs.converged = True
        mock_fs.x0 = []
        mock_fs.var_names = []

        with patch.object(run_mod, "os", self._patch_run_os(run_mod)), \
             patch.object(run_mod, "validate_runtime_flowsheet",
                          return_value=self._containerized_config(flowsheet)), \
             patch.object(run_mod, "check_providers"), \
             patch.object(run_mod, "EOFlowsheet", return_value=mock_fs), \
             patch.object(run_mod, "ProcessStateArchive", MagicMock()):
            run(flowsheet=str(flowsheet), export_images=False)

    def test_coolprop_flowsheet_no_lifecycle(self, tmp_path):
        import processforge.cli.run as run_mod
        from processforge.cli.run import run

        flowsheet = tmp_path / "fs.json"
        flowsheet.write_text(json.dumps({"providers": {"coolprop": {"type": "coolprop"}}}))

        mock_fs = MagicMock()
        mock_fs.run.return_value = {}
        mock_fs.converged = True
        mock_fs.x0 = []
        mock_fs.var_names = []

        # A stale docker-compose.yml is present, but the flowsheet only uses
        # pip providers — Docker must never be touched.
        with patch.object(run_mod, "os", self._patch_run_os(run_mod)), \
             patch.object(run_mod, "validate_runtime_flowsheet",
                          return_value={"providers": {"coolprop": {"type": "coolprop"}},
                                        "streams": {}, "units": {},
                                        "simulation": {"mode": "steady"},
                                        "_config_path": str(flowsheet)}), \
             patch.object(run_mod, "check_providers"), \
             patch.object(run_mod, "EOFlowsheet", return_value=mock_fs), \
             patch.object(run_mod, "ProcessStateArchive", MagicMock()):
            run(flowsheet=str(flowsheet), export_images=False)

    def test_no_compose_no_lifecycle(self, tmp_path):
        import processforge.cli.run as run_mod
        from processforge.cli.run import run

        flowsheet = tmp_path / "fs.json"
        flowsheet.write_text(json.dumps(self._containerized_config(flowsheet)["providers"]))

        mock_fs = MagicMock()
        mock_fs.run.return_value = {}
        mock_fs.converged = True
        mock_fs.x0 = []
        mock_fs.var_names = []

        # Containerized flowsheet but no compose file — run still proceeds
        # (containers are assumed already running).
        with patch.object(run_mod, "os") as mock_os, \
             patch.object(run_mod, "validate_runtime_flowsheet",
                          return_value=self._containerized_config(flowsheet)), \
             patch.object(run_mod, "check_providers"), \
             patch.object(run_mod, "EOFlowsheet", return_value=mock_fs), \
             patch.object(run_mod, "ProcessStateArchive", MagicMock()):
            mock_os.path.exists.return_value = False  # no compose file
            run(flowsheet=str(flowsheet), export_images=False)


class TestCmdPlan:
    @pytest.fixture
    def container_flowsheet(self, tmp_path):
        """Write a minimal valid flowsheet using a Docker-containerized provider."""
        flowsheet = {
            "providers": {"openmc": {"type": "openmc"}},
            "streams": {},
            "units": {},
            "simulation": {"mode": "steady"},
        }
        path = tmp_path / "container.json"
        with open(path, "w") as f:
            json.dump(flowsheet, f)
        return path

    def test_container_provider_health_ok(self, container_flowsheet, log_output):
        import processforge.cli.common as common_mod
        from processforge.cli.plan import plan

        health = {"status": "ready", "provider_type": "openmc"}
        with patch.object(
            common_mod, "_ping_provider_health", return_value=(True, health)
        ), patch.object(common_mod, "_resolve_provider_url", return_value="http://localhost:9001"):
            plan(flowsheet=str(container_flowsheet), no_diagram=True)

        assert any("Provider / Container Health" in m for m in log_output)
        assert any("[OK] openmc" in m and "provider_type=openmc" in m for m in log_output)

    def test_container_provider_health_failure_exits(self, container_flowsheet, log_output):
        import processforge.cli.common as common_mod
        from processforge.cli.plan import plan

        with patch.object(
            common_mod, "_ping_provider_health", return_value=(False, "Connection refused")
        ), patch.object(common_mod, "_resolve_provider_url", return_value="http://localhost:9001"), \
                pytest.raises(SystemExit) as exc_info:
            plan(flowsheet=str(container_flowsheet), no_diagram=True)

        assert exc_info.value.code == 1
        assert any("[ERR] openmc" in m and "unreachable" in m for m in log_output)


class TestCmdApply:
    def _containerized_config(self, flowsheet_path):
        return {
            "providers": {"openmc": {"type": "openmc"}},
            "streams": {},
            "units": {},
            "simulation": {"mode": "steady"},
            "_config_path": str(flowsheet_path),
        }

    def _patch_apply_os(self, apply_mod):
        real_os = apply_mod.os
        mock_os = MagicMock()
        mock_os.path.exists.return_value = True
        mock_os.makedirs = real_os.makedirs
        mock_os.path.join = real_os.path.join
        return mock_os

    def test_containerized_flowsheet_runs_without_lifecycle(self, tmp_path, log_output):
        import processforge.cli.apply as apply_mod
        from processforge.cli.apply import apply

        flowsheet = tmp_path / "fs.json"
        flowsheet.write_text(json.dumps(self._containerized_config(flowsheet)["providers"]))

        mock_fs = MagicMock()
        mock_fs.run.return_value = {}
        mock_fs.converged = True
        mock_fs.x0 = []
        mock_fs.x_converged = []
        mock_fs.var_names = []
        mock_fs.solver_stats = {"final_norm": 0.0}

        with patch.object(apply_mod, "os", self._patch_apply_os(apply_mod)), \
             patch.object(apply_mod, "validate_runtime_flowsheet",
                          return_value=self._containerized_config(flowsheet)), \
             patch.object(apply_mod, "check_providers"), \
             patch.object(apply_mod, "EOFlowsheet", return_value=mock_fs), \
             patch.object(apply_mod, "save_snapshot", return_value="snap-1"), \
             patch.object(apply_mod, "load_state_manager", return_value=(MagicMock(), None)):
            apply(flowsheet=str(flowsheet))

    def test_coolprop_flowsheet_no_lifecycle(self, tmp_path):
        import processforge.cli.apply as apply_mod
        from processforge.cli.apply import apply

        flowsheet = tmp_path / "fs.json"
        flowsheet.write_text(json.dumps({"providers": {"coolprop": {"type": "coolprop"}}}))

        mock_fs = MagicMock()
        mock_fs.run.return_value = {}
        mock_fs.converged = True
        mock_fs.x0 = []
        mock_fs.x_converged = []
        mock_fs.var_names = []
        mock_fs.solver_stats = {"final_norm": 0.0}

        # A stale docker-compose.yml is present, but the flowsheet only uses
        # pip providers — Docker must never be touched.
        with patch.object(apply_mod, "os", self._patch_apply_os(apply_mod)), \
             patch.object(apply_mod, "validate_runtime_flowsheet",
                          return_value={"providers": {"coolprop": {"type": "coolprop"}},
                                        "streams": {}, "units": {},
                                        "simulation": {"mode": "steady"},
                                        "_config_path": str(flowsheet)}), \
             patch.object(apply_mod, "check_providers"), \
             patch.object(apply_mod, "EOFlowsheet", return_value=mock_fs), \
             patch.object(apply_mod, "save_snapshot", return_value="snap-1"), \
             patch.object(apply_mod, "load_state_manager", return_value=(MagicMock(), None)):
            apply(flowsheet=str(flowsheet))


# ---------------------------------------------------------------------------
# cli/diagram.py
# ---------------------------------------------------------------------------

class TestCmdDiagram:
    def test_missing_file_exits(self, tmp_path):
        from processforge.cli.diagram import diagram

        with pytest.raises(SystemExit) as exc_info:
            diagram(
                flowsheet=str(tmp_path / "nope.json"),
                output_dir="diagrams",
                format="png",
            )
        assert exc_info.value.code == 1

    def test_invalid_json_exits(self, invalid_json):
        from processforge.cli.diagram import diagram

        with pytest.raises(SystemExit) as exc_info:
            diagram(
                flowsheet=str(invalid_json),
                output_dir="diagrams",
                format="png",
            )
        assert exc_info.value.code == 1
