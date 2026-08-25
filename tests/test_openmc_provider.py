"""Tests for the OpenMC provider using a fake ``openmc`` module.

The provider only imports ``openmc`` lazily (inside ``initialize`` /
``run_simulation``), so a stub registered in ``sys.modules`` is enough to
exercise the full build → run → extract pipeline without the real package.
"""

import os
import pathlib
import sys
import types
from types import SimpleNamespace

import pandas as pd
import pytest
from pydantic import ValidationError

from processforge.providers.openmc_provider import OpenMCProvider
from processforge.schemas.openmc.openmc_model import (
    OpenMCSetting,
    ReactorCoreGeometryConfig,
)
from processforge.types import MaterialDef, OpenMCProviderConfig, UnitConfig


# ---------------------------------------------------------------------------
# Fake openmc module
# ---------------------------------------------------------------------------


class _RunControl:
    result = None
    raise_msg = None
    sp_keff = None
    sp_tallies = None
    calls = []


_CTRL = _RunControl()


class _FakeMaterial:
    def __init__(self, material_id, name, **kwargs):
        self.id = material_id
        self.name = name
        self.temperature = None
        self.depletable = False
        self.density = None
        self.density_units = None
        self._nuclides = []
        self._elements = []

    def set_density(self, units, density):
        self.density_units = units
        self.density = density

    def add_nuclide(self, name, percent, percent_type="ao"):
        self._nuclides.append((name, percent, percent_type))

    def add_element(self, element, percent, percent_type="ao"):
        self._elements.append((element, percent, percent_type))


class _FakeMaterials:
    def __init__(self, materials=None):
        self._materials = list(materials or [])


class _FakeRegularMesh:
    pass


class _FakeMeshFilter:
    def __init__(self, mesh):
        self.mesh = mesh


class _FakeTally:
    def __init__(self, tally_id, name=None):
        self.id = tally_id
        self.name = name
        self.filters = []
        self.scores = ["flux"]
        self.nuclides = None
        self.estimator = "tracklength"

    def get_pandas_dataframe(self, scores=None):
        return pd.DataFrame({"mean": [10.0, 20.0], "std. dev.": [1.0, 3.0]})


class _FakeTallies:
    def __init__(self, tallies=None):
        self._tallies = list(tallies or [])


class _FakePoint:
    def __init__(self, xyz):
        self.xyz = xyz


class _FakeIsotropic:
    pass


class _FakeWatt:
    pass


class _FakeDiscrete:
    def __init__(self, values, probabilities):
        self.values = values
        self.probabilities = probabilities


class _FakeIndependentSource:
    pass


class _FakeSettings:
    pass


class _FakeSphere:
    def __init__(self, r=None, boundary_type=None):
        self.r = r
        self.boundary_type = boundary_type

    def __neg__(self):
        return self

    def __pos__(self):
        return self


class _FakeRegion:
    def __init__(self, surf=None, neg=False):
        self.surf = surf
        self.neg = neg

    def __or__(self, other):
        return _FakeRegion()

    def __and__(self, other):
        return _FakeRegion()

    def __pos__(self):
        return self

    def __neg__(self):
        return self


class _FakeSurface:
    def __init__(self, *args, **kwargs):
        pass

    def __pos__(self):
        return _FakeRegion(self)

    def __neg__(self):
        return _FakeRegion(self, neg=True)


class _FakeCell:
    def __init__(self, fill=None, region=None):
        self.fill = fill
        self.region = region


class _FakeGeometry:
    def __init__(self, cells=None):
        self._cells = list(cells or [])


class _FakeStatePoint:
    def __init__(self, path, tallies=None, keff=None):
        self.path = path
        self._tallies = dict(
            tallies if tallies is not None else (_CTRL.sp_tallies or {})
        )
        self.keff = keff if keff is not None else _CTRL.sp_keff

    def get_tally(self, id):
        try:
            return self._tallies[id]
        except KeyError as exc:
            raise KeyError(f"tally {id} not found in statepoint") from exc


class _FakeModel:
    def __init__(self, geometry, materials, settings, tallies=None):
        self.geometry = geometry
        self.materials = materials
        self.settings = settings
        self.tallies = tallies

    def run(self, cwd, export_model_xml=False, **kwargs):
        _CTRL.calls.append(
            {
                "cwd": str(cwd),
                "export_model_xml": export_model_xml,
                "xs_env": os.environ.get("OPENMC_CROSS_SECTIONS"),
                "run_mode": getattr(self.settings, "run_mode", None),
            }
        )
        if _CTRL.raise_msg:
            raise RuntimeError(_CTRL.raise_msg)
        if _CTRL.result is None:
            return None
        return pathlib.Path(_CTRL.result)


def make_fake_openmc() -> types.ModuleType:
    mod = types.ModuleType("openmc")
    mod.Material = _FakeMaterial
    mod.Materials = _FakeMaterials
    mod.RegularMesh = _FakeRegularMesh
    mod.MeshFilter = _FakeMeshFilter
    mod.Tally = _FakeTally
    mod.Tallies = _FakeTallies
    mod.IndependentSource = _FakeIndependentSource
    mod.Settings = _FakeSettings
    mod.Sphere = _FakeSphere
    mod.ZCylinder = _FakeSurface
    mod.XPlane = _FakeSurface
    mod.YPlane = _FakeSurface
    mod.ZPlane = _FakeSurface
    mod.Union = _FakeRegion
    mod.Intersection = _FakeRegion
    mod.Cell = _FakeCell
    mod.Geometry = _FakeGeometry
    mod.StatePoint = _FakeStatePoint
    mod.stats = SimpleNamespace(
        Point=_FakePoint,
        Isotropic=_FakeIsotropic,
        Watt=_FakeWatt,
        Discrete=_FakeDiscrete,
    )
    mod.model = SimpleNamespace(Model=_FakeModel)
    return mod


@pytest.fixture
def fake_openmc(monkeypatch):
    fake = make_fake_openmc()
    monkeypatch.setitem(sys.modules, "openmc", fake)
    _CTRL.calls.clear()
    _CTRL.result = None
    _CTRL.raise_msg = None
    _CTRL.sp_keff = None
    _CTRL.sp_tallies = None
    return fake


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mat(id_, name, **overrides):
    base = {
        "id": id_,
        "density": 2.2,
        "density_units": "g/cm3",
        "nuclides": [{"name": "Li7", "percent": 1.0, "percent_type": "ao"}],
    }
    base.update(overrides)
    return MaterialDef.from_dict(base)


def _point_geometry(material="salt", **overrides):
    geo = {
        "source_point": {"xyz": [0.0, 0.0, 0.0]},
        "sphere_radius": 200.0,
        "sphere_material": material,
    }
    geo.update(overrides)
    return geo


def _unit_config(sim_type="eigenvalue_csg", **sc_overrides):
    sc = {
        "batches": 20,
        "inactive": 5,
        "particles": 1000,
        "mesh_tallies": [
            {
                "tally_id": 1,
                "name": "flux",
                "lower_left": [-1, -1, -1],
                "upper_right": [1, 1, 1],
                "dimension": [2, 2, 2],
                "scores": ["flux"],
            }
        ],
    }
    sc.update(sc_overrides)
    return UnitConfig.from_dict(
        {
            "type": "SolverUnit",
            "sim_type": sim_type,
            "solver_config": sc,
            "geometry_config": _point_geometry(),
        }
    )


def _unit_config_reactor(**sc_overrides):
    sc = {"batches": 20, "inactive": 5, "particles": 1000}
    sc.update(sc_overrides)
    geo = {
        "type": "reactor_core",
        "core_radius": 72.5,
        "core_height": 160.0,
        "reflector_thickness": 50.0,
        "core_material": "salt",
        "reflector_material": "graphite",
        "vessel_material": "inconel",
        "source_point": {"xyz": [0.0, 0.0, 0.0]},
    }
    return UnitConfig.from_dict(
        {
            "type": "SolverUnit",
            "sim_type": "eigenvalue_reactor",
            "solver_config": sc,
            "geometry_config": geo,
        }
    )


def _default_materials():
    return {
        "salt": _mat(3, "salt"),
        "graphite": _mat(1, "graphite"),
        "inconel": _mat(4, "inconel"),
    }


def _init_provider(tmp_path, materials, cross_sections=None):
    provider = OpenMCProvider()
    cfg = OpenMCProviderConfig(
        output_dir=str(tmp_path / "openmc_run"), cross_sections=cross_sections
    )
    flowsheet = SimpleNamespace(materials=materials)
    provider.initialize(cfg, flowsheet)
    return provider


# ---------------------------------------------------------------------------
# run_simulation behaviour
# ---------------------------------------------------------------------------


def test_run_eigenvalue_extracts_keff_and_tallies(tmp_path, fake_openmc):
    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.result = str(tmp_path / "openmc_run" / "statepoint.7.h5")
    _CTRL.sp_keff = SimpleNamespace(n=1.02, s=0.004)
    _CTRL.sp_tallies = {1: _FakeTally(1, "flux")}

    result = provider.run_simulation(_unit_config(), {})

    assert result.status == "completed"
    assert result.field_value("k_eff") == pytest.approx(1.02)
    assert result.get_field("k_eff").quantity.std_dev == pytest.approx(0.004)
    # Fake mesh exposes no geometry -> cell volume falls back to 1.0, so the
    # volume-weighted integral is the sum of the per-cell means (10 + 20 = 30)
    # and the mean is 15.0.
    assert result.field_value("tally_1_flux_mean") == pytest.approx(15.0)
    assert result.field_value("tally_1_flux_integrated") == pytest.approx(30.0)
    assert result.diagnostics["run_dir"] == str((tmp_path / "openmc_run").resolve())
    assert result.diagnostics["statepoint_path"].endswith("statepoint.7.h5")


def test_statepoint_path_comes_from_model_run_return(tmp_path, fake_openmc):
    """A stale batches count must not break statepoint discovery."""
    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.result = str(tmp_path / "openmc_run" / "statepoint.7.h5")
    result = provider.run_simulation(_unit_config(batches=999), {})

    assert result.status == "completed"
    assert result.diagnostics["statepoint_path"].endswith("statepoint.7.h5")


def test_cwd_never_changes(tmp_path, fake_openmc):
    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.result = str(tmp_path / "openmc_run" / "statepoint.20.h5")
    cwd_before = pathlib.Path.cwd()

    provider.run_simulation(_unit_config(), {})

    assert pathlib.Path.cwd() == cwd_before
    assert _CTRL.calls[0]["cwd"] == str(tmp_path / "openmc_run")
    assert _CTRL.calls[0]["export_model_xml"] is False


def test_model_run_uses_configured_run_mode(tmp_path, fake_openmc):
    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.result = str(tmp_path / "openmc_run" / "statepoint.20.h5")

    provider.run_simulation(_unit_config("fixed_source_point"), {})

    assert _CTRL.calls[-1]["run_mode"] == "fixed source"

    provider.run_simulation(_unit_config("eigenvalue_csg"), {})

    assert _CTRL.calls[-1]["run_mode"] == "eigenvalue"


def test_cross_sections_env_set_during_run_and_restored(tmp_path, fake_openmc):
    xs_file = tmp_path / "cross_sections.xml"
    xs_file.write_text("<cross_sections/>")
    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.result = str(tmp_path / "openmc_run" / "statepoint.20.h5")

    prior = os.environ.get("OPENMC_CROSS_SECTIONS")
    provider.run_simulation(_unit_config(cross_sections=str(xs_file)), {})

    assert _CTRL.calls[0]["xs_env"] == str(xs_file)
    assert os.environ.get("OPENMC_CROSS_SECTIONS") == prior


def test_run_exception_returns_failed_result(tmp_path, fake_openmc):
    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.raise_msg = "boom: no cross sections loaded"

    prior = os.environ.get("OPENMC_CROSS_SECTIONS")
    xs_file = tmp_path / "custom" / "xs.xml"
    xs_file.parent.mkdir(parents=True)
    xs_file.write_text("<cross_sections/>")
    result = provider.run_simulation(_unit_config(cross_sections=str(xs_file)), {})

    assert result.status == "failed"
    assert "boom" in result.diagnostics["error"]
    assert result.diagnostics["run_dir"] == str((tmp_path / "openmc_run").resolve())
    assert os.environ.get("OPENMC_CROSS_SECTIONS") == prior


def test_run_exception_is_classified(tmp_path, fake_openmc):
    """A nuclear-data run failure must be classified with source='provider'."""
    from processforge.providers.errors import NUCLEAR_DATA

    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.raise_msg = (
        "Nuclear data library does not contain cross sections for Ni58 at or near "
        "400.000000 K. Available temperatures are 300 K."
    )
    result = provider.run_simulation(_unit_config(), {})

    assert result.status == "failed"
    assert result.error is not None
    assert result.error.category == NUCLEAR_DATA
    assert result.error.source == "provider"
    assert result.diagnostics["error_category"] == NUCLEAR_DATA
    assert "temperature" in result.error.hint.lower()


def test_run_without_statepoint_returns_failed(tmp_path, fake_openmc):
    provider = _init_provider(tmp_path, _default_materials())
    result = provider.run_simulation(_unit_config(), {})

    assert result.status == "failed"
    assert "no statepoint" in result.diagnostics["error"]


def test_missing_tally_recorded_as_warning(tmp_path, fake_openmc):
    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.result = str(tmp_path / "openmc_run" / "statepoint.20.h5")

    # Statepoint exposes no tallies — extraction should degrade with a warning.
    result = provider.run_simulation(
        _unit_config(
            mesh_tallies=[
                {
                    "tally_id": 7,
                    "name": "ghost",
                    "lower_left": [-1, -1, -1],
                    "upper_right": [1, 1, 1],
                    "dimension": [2, 2, 2],
                    "scores": ["flux"],
                }
            ]
        ),
        {},
    )

    assert result.status == "completed"
    assert any("tally id=7" in w for w in result.diagnostics["tally_warnings"])


# ---------------------------------------------------------------------------
# initialize validation
# ---------------------------------------------------------------------------


def test_initialize_rejects_duplicate_material_ids(tmp_path, fake_openmc):
    materials = {
        "salt_a": _mat(3, "salt_a"),
        "salt_b": _mat(3, "salt_b"),
    }
    with pytest.raises(ValueError, match="id=3"):
        _init_provider(tmp_path, materials)


def test_initialize_rejects_missing_cross_sections(tmp_path, fake_openmc):
    missing = str(tmp_path / "does_not_exist" / "cross_sections.xml")
    with pytest.raises(RuntimeError, match="cross-sections file not found"):
        _init_provider(tmp_path, _default_materials(), cross_sections=missing)


# ---------------------------------------------------------------------------
# validate_material
# ---------------------------------------------------------------------------


def test_validate_material_accepts_valid(tmp_path):
    errors = OpenMCProvider.validate_material(
        "salt", _mat(3, "salt", density_units="sum"), None
    )
    assert errors == []


def test_validate_material_requires_density_unless_sum(tmp_path):
    errors = OpenMCProvider.validate_material(
        "salt", _mat(3, "salt", density=None), None
    )
    assert any("density" in e for e in errors)


def test_validate_material_rejects_invalid_density_units(tmp_path):
    errors = OpenMCProvider.validate_material(
        "salt", _mat(3, "salt", density_units="lightyears"), None
    )
    assert any("density_units" in e for e in errors)


def test_validate_material_rejects_no_components(tmp_path):
    errors = OpenMCProvider.validate_material(
        "salt", _mat(3, "salt", nuclides=[], elements=[]), None
    )
    assert any("nuclides or elements" in e for e in errors)


def test_validate_material_rejects_malformed_nuclides(tmp_path):
    mdef = _mat(3, "salt", nuclides=[{"percent": 1.0}])
    errors = OpenMCProvider.validate_material("salt", mdef, None)
    assert any("missing 'name'" in e for e in errors)


def test_validate_material_rejects_bad_percent_type(tmp_path):
    mdef = _mat(
        3, "salt", nuclides=[{"name": "Li7", "percent": 1.0, "percent_type": "wt"}]
    )
    errors = OpenMCProvider.validate_material("salt", mdef, None)
    assert any("percent_type" in e for e in errors)


def test_validate_material_rejects_non_dict_element(tmp_path):
    mdef = _mat(3, "salt", elements=["not-a-dict"])
    errors = OpenMCProvider.validate_material("salt", mdef, None)
    assert any("elements[0]" in e for e in errors)


# ---------------------------------------------------------------------------
# OpenMCSetting schema validation
# ---------------------------------------------------------------------------


def test_openmc_setting_minimal_defaults():
    cfg = OpenMCSetting()
    assert cfg.batches == 20
    assert cfg.inactive == 5
    assert cfg.particles == 1000
    assert cfg.run_mode.value == "eigenvalue"


def test_openmc_setting_rejects_invalid_run_mode():
    with pytest.raises(ValidationError):
        OpenMCSetting(run_mode="fixed_source_with_underscores")


def test_openmc_setting_accepts_enum_string_coercion():
    cfg = OpenMCSetting(run_mode="fixed source")
    assert cfg.run_mode.value == "fixed source"


def test_openmc_setting_resolution_preset_overrides():
    cfg = OpenMCSetting(resolution="high")
    assert cfg.batches == 50
    assert cfg.inactive == 10
    assert cfg.particles == 100000


def test_openmc_setting_rejects_inconsistent_mesh_dims():
    with pytest.raises(ValidationError):
        OpenMCSetting(
            mesh_tallies=[
                {
                    "tally_id": 1,
                    "lower_left": [0, 0, 0],
                    "upper_right": [1, 1],
                    "dimension": [2, 2, 2],
                    "scores": ["flux"],
                }
            ]
        )


def test_openmc_setting_rejects_duplicate_tally_ids():
    tally = {
        "tally_id": 1,
        "lower_left": [0, 0, 0],
        "upper_right": [1, 1, 1],
        "dimension": [2, 2, 2],
        "scores": ["flux"],
    }
    with pytest.raises(ValidationError):
        OpenMCSetting(mesh_tallies=[tally, dict(tally)])


def test_openmc_setting_rejects_empty_scores():
    with pytest.raises(ValidationError):
        OpenMCSetting(
            mesh_tallies=[
                {
                    "tally_id": 1,
                    "lower_left": [0, 0, 0],
                    "upper_right": [1, 1, 1],
                    "dimension": [2, 2, 2],
                    "scores": [],
                }
            ]
        )


# ---------------------------------------------------------------------------
# ReactorCoreGeometryConfig validation
# ---------------------------------------------------------------------------


def test_reactor_core_geometry_rejects_unknown_material():
    with pytest.raises(ValidationError):
        ReactorCoreGeometryConfig.model_validate(
            {
                "core_radius": 72.5,
                "core_height": 160.0,
                "core_material": "nope",
                "reflector_material": "graphite",
                "vessel_material": "inconel",
            },
            context={"materials": {"graphite", "inconel"}},
        )


def test_reactor_core_geometry_accepts_valid():
    cfg = ReactorCoreGeometryConfig.model_validate(
        {
            "core_radius": 72.5,
            "core_height": 160.0,
            "core_material": "salt",
            "reflector_material": "graphite",
            "vessel_material": "inconel",
        },
        context={"materials": {"salt", "graphite", "inconel"}},
    )
    assert cfg.core_material == "salt"


def test_run_reactor_core_builds_and_runs(tmp_path, fake_openmc):
    """The approximate reactor geometry builds and runs in eigenvalue mode."""
    provider = _init_provider(tmp_path, _default_materials())
    _CTRL.result = str(tmp_path / "openmc_run" / "statepoint.20.h5")
    _CTRL.sp_keff = SimpleNamespace(n=1.01, s=0.003)

    result = provider.run_simulation(_unit_config_reactor(), {})

    assert result.status == "completed"
    assert result.field_value("k_eff") == pytest.approx(1.01)
    assert _CTRL.calls[-1]["run_mode"] == "eigenvalue"
