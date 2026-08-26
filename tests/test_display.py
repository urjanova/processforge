"""Tests for CLI display helpers (drift / structural diff formatting)."""
import loguru

from processforge.cli.display import print_param_drift, _resolve_path


def _capture_logs():
    """Install a loguru sink that records emitted message strings."""
    records: list[str] = []

    def sink(message):
        records.append(str(message))

    loguru.logger.add(sink, format="{message}", level="INFO")
    return records


def test_resolve_path_nested():
    cfg = {"simulation": {"tf": 25.0}, "units": {"p": {"solver_config": {"batches": 20}}}}
    assert _resolve_path(cfg, "simulation.tf") == 25.0
    assert _resolve_path(cfg, "units.p.solver_config.batches") == 20
    assert _resolve_path(cfg, "units.p.missing") == "<missing>"


def test_print_param_drift_no_changes():
    records = _capture_logs()
    print_param_drift([], {}, {})
    assert any("no parameter changes" in r for r in records)


def test_print_param_drift_shows_values():
    records = _capture_logs()
    old = {"simulation": {"tf": 25.0}}
    new = {"simulation": {"tf": 40.0}}
    print_param_drift(["simulation.tf"], old, new)
    assert any("simulation.tf: 25.0 → 40.0" in r for r in records)
