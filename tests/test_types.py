from __future__ import annotations

import numpy as np
import pytest

from processforge.types import (
    CanteraProviderConfig,
    CoolPropProviderConfig,
    FestimProviderConfig,
    FlowsheetConfig,
    MaterialDef,
    ModelicaProviderConfig,
    OpenMCProviderConfig,
    RunInfo,
    SimulationResult,
    SnapshotState,
    UnitConfig,
    provider_config_from_dict,
)


# ---------------------------------------------------------------------------
# MaterialDef
# ---------------------------------------------------------------------------

class TestMaterialDef:
    def test_direct_construction(self):
        m = MaterialDef(id=1, density=2.0, extra={"D_0": 1e-4})
        assert m.id == 1
        assert m.density == 2.0
        assert m.extra == {"D_0": 1e-4}

    def test_from_dict_separates_known_keys(self):
        m = MaterialDef.from_dict({"id": 1, "density": 2.0, "D_0": 1e-4})
        assert m.id == 1
        assert m.density == 2.0
        assert m.extra == {"D_0": 1e-4}

    def test_from_dict_flattens_nested_extra(self):
        m = MaterialDef.from_dict({"id": 1, "extra": {"D_0": 1e-4}, "E_act": 0.5})
        assert m.id == 1
        assert m.extra == {"D_0": 1e-4, "E_act": 0.5}

    def test_from_dict_nested_extra_is_overridden_by_flat(self):
        m = MaterialDef.from_dict({"id": 1, "extra": {"D_0": 1e-4}, "D_0": 0.5})
        assert m.extra["D_0"] == 0.5

    def test_from_dict_empty(self):
        m = MaterialDef.from_dict({"id": 1})
        assert m.id == 1
        assert m.density is None
        assert m.extra == {}

    def test_get_typed_field(self):
        m = MaterialDef(id=1, density=2.0)
        assert m.get("id") == 1
        assert m.get("density") == 2.0

    def test_get_extra_field(self):
        m = MaterialDef(id=1, extra={"D_0": 1e-4})
        assert m.get("D_0") == 1e-4

    def test_get_returns_default(self):
        m = MaterialDef(id=1)
        assert m.get("nonexistent", "fallback") == "fallback"
        assert m.get("nonexistent") is None

    def test_defaults(self):
        m = MaterialDef(id=1)
        assert m.depletable is False
        assert m.nuclides == []
        assert m.extra == {}

    def test_nuclides_list(self):
        m = MaterialDef(id=1, nuclides=["H1", "O16"])
        assert m.nuclides == ["H1", "O16"]


# ---------------------------------------------------------------------------
# UnitConfig
# ---------------------------------------------------------------------------

class TestUnitConfig:
    def test_direct_construction(self):
        u = UnitConfig(type="SolverUnit")
        assert u.type == "SolverUnit"

    def test_from_dict_string_in_becomes_single_input(self):
        u = UnitConfig.from_dict({"type": "SolverUnit", "in": "stream_1"})
        assert u.inputs == ["stream_1"]

    def test_from_dict_list_in_passes_through(self):
        u = UnitConfig.from_dict({"type": "SolverUnit", "in": ["a", "b"]})
        assert u.inputs == ["a", "b"]

    def test_from_dict_no_in_defaults_empty(self):
        u = UnitConfig.from_dict({"type": "SolverUnit"})
        assert u.inputs == []

    def test_from_dict_in_is_none(self):
        u = UnitConfig.from_dict({"type": "SolverUnit", "in": None})
        assert u.inputs == []

    def test_from_dict_separates_extra(self):
        u = UnitConfig.from_dict({"type": "SolverUnit", "mesh": "fine", "tolerance": 1e-6})
        assert u.extra == {"mesh": "fine", "tolerance": 1e-6}
        assert u.type == "SolverUnit"

    def test_from_dict_known_fields(self):
        u = UnitConfig.from_dict({
            "type": "SolverUnit",
            "material": 1,
            "provider": "openmc",
            "out": "out_stream",
            "sim_type": "heat_2d",
            "solver_config": {"n_steps": 100},
        })
        assert u.material == 1
        assert u.provider == "openmc"
        assert u.out == "out_stream"
        assert u.sim_type == "heat_2d"
        assert u.solver_config == {"n_steps": 100}

    def test_retentate_and_permeate(self):
        u = UnitConfig.from_dict({
            "type": "SolverUnit",
            "retentate_out": "ret_stream",
            "permeate_out": "perm_stream",
        })
        assert u.retentate_out == "ret_stream"
        assert u.permeate_out == "perm_stream"

    def test_in_not_in_extra(self):
        u = UnitConfig.from_dict({"type": "SolverUnit", "in": ["s1"], "foo": "bar"})
        assert "in" not in u.extra
        assert u.extra == {"foo": "bar"}


# ---------------------------------------------------------------------------
# Provider config classes
# ---------------------------------------------------------------------------

class TestCoolPropProviderConfig:
    def test_type(self):
        assert CoolPropProviderConfig().type == "coolprop"

    def test_from_dict_ignores_input(self):
        c = CoolPropProviderConfig.from_dict({"type": "coolprop", "anything": 42})
        assert c.type == "coolprop"


class TestCanteraProviderConfig:
    def test_defaults(self):
        c = CanteraProviderConfig()
        assert c.type == "cantera"
        assert c.mechanism == "gri30.yaml"
        assert c.phase is None

    def test_from_dict(self):
        c = CanteraProviderConfig.from_dict({"mechanism": "h2o2.yaml", "phase": "gas"})
        assert c.mechanism == "h2o2.yaml"
        assert c.phase == "gas"

    def test_from_dict_partial(self):
        c = CanteraProviderConfig.from_dict({})
        assert c.mechanism == "gri30.yaml"
        assert c.phase is None


class TestModelicaProviderConfig:
    def test_defaults(self):
        m = ModelicaProviderConfig()
        assert m.type == "modelica"
        assert m.output_dir == "outputs"

    def test_from_dict(self):
        m = ModelicaProviderConfig.from_dict({"output_dir": "custom_out"})
        assert m.output_dir == "custom_out"

    def test_from_dict_default(self):
        m = ModelicaProviderConfig.from_dict({})
        assert m.output_dir == "outputs"


class TestOpenMCProviderConfig:
    def test_defaults(self):
        o = OpenMCProviderConfig()
        assert o.type == "openmc"
        assert o.output_dir == "outputs/openmc"
        assert o.cross_sections is None
        assert o.docker_image is None

    def test_from_dict(self):
        o = OpenMCProviderConfig.from_dict({
            "output_dir": "custom",
            "cross_sections": "/path/to/xs.xml",
            "docker_image": "my/img",
        })
        assert o.output_dir == "custom"
        assert o.cross_sections == "/path/to/xs.xml"
        assert o.docker_image == "my/img"

    def test_from_dict_partial(self):
        o = OpenMCProviderConfig.from_dict({})
        assert o.output_dir == "outputs/openmc"
        assert o.cross_sections is None


class TestFestimProviderConfig:
    def test_defaults(self):
        f = FestimProviderConfig()
        assert f.type == "festim"
        assert f.url is None
        assert f.output_dir == "outputs/festim"
        assert f.docker_image is None

    def test_direct_construction(self):
        f = FestimProviderConfig(url="http://localhost:9002")
        assert f.url == "http://localhost:9002"

    def test_from_dict(self):
        f = FestimProviderConfig.from_dict({
            "url": "http://example.com",
            "output_dir": "custom",
            "docker_image": "my/img",
        })
        assert f.url == "http://example.com"
        assert f.output_dir == "custom"
        assert f.docker_image == "my/img"

    def test_from_dict_partial(self):
        f = FestimProviderConfig.from_dict({})
        assert f.output_dir == "outputs/festim"
        assert f.url is None


# ---------------------------------------------------------------------------
# provider_config_from_dict
# ---------------------------------------------------------------------------

class TestProviderConfigFromDict:
    def test_coolprop(self):
        c = provider_config_from_dict({"type": "coolprop"})
        assert isinstance(c, CoolPropProviderConfig)

    def test_cantera(self):
        c = provider_config_from_dict({"type": "cantera", "mechanism": "h2o2.yaml"})
        assert isinstance(c, CanteraProviderConfig)
        assert c.mechanism == "h2o2.yaml"

    def test_modelica(self):
        c = provider_config_from_dict({"type": "modelica"})
        assert isinstance(c, ModelicaProviderConfig)

    def test_openmc(self):
        c = provider_config_from_dict({"type": "openmc", "cross_sections": "/x"})
        assert isinstance(c, OpenMCProviderConfig)
        assert c.cross_sections == "/x"

    def test_festim(self):
        c = provider_config_from_dict({"type": "festim", "url": "http://x"})
        assert isinstance(c, FestimProviderConfig)
        assert c.url == "http://x"

    def test_missing_type_raises(self):
        with pytest.raises(ValueError, match="missing"):
            provider_config_from_dict({"foo": "bar"})

    def test_unknown_type_raises(self):
        with pytest.raises(ValueError, match="Unknown"):
            provider_config_from_dict({"type": "nonexistent"})

    def test_all_configs_have_type(self):
        for ptype in ("coolprop", "cantera", "modelica", "openmc", "festim"):
            c = provider_config_from_dict({"type": ptype})
            assert c.type == ptype


# ---------------------------------------------------------------------------
# FlowsheetConfig
# ---------------------------------------------------------------------------

class TestFlowsheetConfig:
    MINIMAL = {
        "providers": {"coolprop": {"type": "coolprop"}},
        "units": {"mixer": {"type": "SolverUnit", "in": ["a"]}},
        "materials": {"tungsten": {"id": 1, "density": 19.3}},
    }

    def test_from_dict_minimal(self):
        cfg = FlowsheetConfig.from_dict(self.MINIMAL)
        assert isinstance(cfg.providers["coolprop"], CoolPropProviderConfig)
        assert isinstance(cfg.units["mixer"], UnitConfig)
        assert isinstance(cfg.materials["tungsten"], MaterialDef)
        assert cfg.materials["tungsten"].density == 19.3

    def test_from_dict_extra(self):
        raw = {**self.MINIMAL, "_config_path": "/tmp/cfg.json"}
        cfg = FlowsheetConfig.from_dict(raw)
        assert cfg.extra["_config_path"] == "/tmp/cfg.json"

    def test_get_known_field(self):
        raw = {**self.MINIMAL, "default_provider": "coolprop"}
        cfg = FlowsheetConfig.from_dict(raw)
        assert cfg.get("default_provider") == "coolprop"

    def test_get_extra_field(self):
        raw = {**self.MINIMAL, "_runtime_key": 42}
        cfg = FlowsheetConfig.from_dict(raw)
        assert cfg.get("_runtime_key") == 42

    def test_get_returns_default(self):
        cfg = FlowsheetConfig.from_dict(self.MINIMAL)
        assert cfg.get("nonexistent", "fallback") == "fallback"
        assert cfg.get("nonexistent") is None

    def test_defaults(self):
        cfg = FlowsheetConfig()
        assert cfg.providers == {}
        assert cfg.units == {}
        assert cfg.materials == {}
        assert cfg.streams == {}
        assert cfg.metadata is None


# ---------------------------------------------------------------------------
# SimulationResult
# ---------------------------------------------------------------------------

class TestSimulationResult:
    def test_construction(self):
        r = SimulationResult(status="completed", sim_type="heat_2d")
        assert r.status == "completed"
        assert r.sim_type == "heat_2d"

    def test_as_dict_flattens_scalars(self):
        r = SimulationResult(
            status="completed",
            sim_type="heat_2d",
            scalars={"T_max": 1200.0, "flux": 0.5},
        )
        d = r.as_dict()
        assert d == {"status": "completed", "sim_type": "heat_2d", "T_max": 1200.0, "flux": 0.5}

    def test_as_dict_excludes_metadata(self):
        r = SimulationResult(
            status="completed",
            sim_type="heat_2d",
            metadata={"xdmf_files": ["out.xdmf"]},
        )
        d = r.as_dict()
        assert "metadata" not in d


# ---------------------------------------------------------------------------
# RunInfo
# ---------------------------------------------------------------------------

class TestRunInfo:
    REQUIRED = dict(
        git_hash="abc123",
        timestamp="2024-01-01T00:00:00",
        python_version="3.12",
        platform="linux",
        processforge_version="0.1.0",
        mode="run",
        backend="local",
    )

    def test_construction(self):
        r = RunInfo(**self.REQUIRED)
        for k, v in self.REQUIRED.items():
            assert getattr(r, k) == v

    def test_as_dict(self):
        r = RunInfo(**self.REQUIRED, x0=[1.0, 2.0], var_names=["a", "b"])
        d = r.as_dict()
        assert d["x0"] == [1.0, 2.0]
        assert d["var_names"] == ["a", "b"]
        assert d["pkg_versions"] == {}

    def test_as_dict_none_fields(self):
        r = RunInfo(**self.REQUIRED)
        d = r.as_dict()
        assert d["x0"] is None
        assert d["var_names"] is None

    def test_isinstance(self):
        r = RunInfo(**self.REQUIRED)
        assert isinstance(r, RunInfo)


# ---------------------------------------------------------------------------
# SnapshotState
# ---------------------------------------------------------------------------

class TestSnapshotState:
    def test_construction_minimal(self):
        s = SnapshotState(config={"a": 1}, x=[1.0, 2.0])
        assert s.config == {"a": 1}
        assert s.x == [1.0, 2.0]
        assert s.var_names == []
        assert s.snapshot_id == ""

    def test_construction_with_x_delta(self):
        arr = np.array([1.0, 2.0, 3.0])
        s = SnapshotState(config={}, x=[], x_delta=arr)
        assert np.array_equal(s.x_delta, arr)

    def test_x_delta_default(self):
        s = SnapshotState(config={}, x=[])
        assert isinstance(s.x_delta, np.ndarray)
        assert len(s.x_delta) == 0

    def test_as_dict(self):
        s = SnapshotState(
            config={"k": "v"},
            x=[1.0],
            var_names=["a"],
            snapshot_id="snap1",
            timestamp="t1",
            metadata={"origin": "test"},
            x_delta=np.array([0.1, 0.2]),
            parent_snapshot_id="parent1",
        )
        d = s.as_dict()
        assert d["config"] == {"k": "v"}
        assert d["x"] == [1.0]
        assert d["var_names"] == ["a"]
        assert d["x_delta"] == [0.1, 0.2]
        assert d["parent_snapshot_id"] == "parent1"

    def test_as_dict_empty_x_delta(self):
        s = SnapshotState(config={}, x=[])
        d = s.as_dict()
        assert d["x_delta"] == []

    def test_isinstance(self):
        s = SnapshotState(config={}, x=[])
        assert isinstance(s, SnapshotState)
