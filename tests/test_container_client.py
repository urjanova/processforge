"""Tests for the CLI-side container provider HTTP client serialization."""

import json
import urllib.request

from processforge.providers.container_client import ContainerProviderClient
from processforge.types import UnitConfig


def test_serialize_unit_config_includes_geometry_config():
    """geometry_config must survive the CLI -> container body serialization.

    Regression test for the bug where geometry_config (a known UnitConfig
    field) was dropped by the hardcoded serialization list, so the container
    received an empty geometry_config and failed validation.
    """
    geometry_config = {
        "type": "reactor_core",
        "core_radius": 72.5,
        "core_height": 160.0,
        "core_material": "salt",
    }
    uc = UnitConfig(
        type="SolverUnit",
        provider="openmc",
        material=3,
        sim_type="eigenvalue_reactor",
        solver_config={"batches": 20},
        geometry_config=geometry_config,
    )

    body = ContainerProviderClient._serialize_unit_config(uc)

    assert body["geometry_config"] == geometry_config
    assert body["solver_config"] == {"batches": 20}
    assert body["sim_type"] == "eigenvalue_reactor"


def test_run_simulation_body_excludes_output_dir(monkeypatch):
    """The CLI must not dictate the container's scratch path via output_dir.

    The server resolves its own run folder under PROCESSFORGE_OUTPUT_DIR, so
    the top-level body field is removed (it was never consumed and caused a
    client/server path mismatch).
    """
    captured = {}

    class _Resp:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def read(self):
            return json.dumps(
                {
                    "status": "completed",
                    "engine": "openmc",
                    "sim_type": "eigenvalue_reactor",
                    "fields": [],
                    "artifacts": [],
                    "diagnostics": {},
                }
            ).encode()

    def _urlopen(req, timeout=0):
        if req.full_url.endswith("/health"):
            resp = _Resp()
            resp.read = lambda: json.dumps({"status": "ready"}).encode()
            return resp
        captured["body"] = json.loads(req.data.decode())
        return _Resp()

    monkeypatch.setattr(urllib.request, "urlopen", _urlopen)

    client = ContainerProviderClient("openmc")
    client.initialize(
        type("Cfg", (), {"url": "http://localhost:9000", "output_dir": "outputs/openmc", "model_dump": lambda self: {"url": "http://localhost:9000", "output_dir": "outputs/openmc"}})(),
        type("FS", (), {"materials": {}})(),
    )
    client.run_simulation(
        UnitConfig(
            type="SolverUnit",
            provider="openmc",
            material=3,
            sim_type="eigenvalue_reactor",
        ),
        inlet={},
    )

    assert "output_dir" not in captured["body"]
    assert "provider_config" in captured["body"]
