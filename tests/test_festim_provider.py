"""Tests for the FESTIM provider (in-process architecture).

These tests validate:
- The provider registers correctly.
- The provider satisfies the AbstractProvider contract.
- Material validation catches missing/invalid D_0/E_D.
- FestimProviderConfig parses from dict correctly.
- The strategy registry and dispatch work correctly.
- The strategy builds a ``HydrogenTransportProblem`` from the typed
  ``FestimModel`` (via a ``fake_festim`` module injected into ``sys.modules``).
- ``run_simulation`` runs the problem in the resolved run directory and
  extracts scalars + time-series from the exports.
"""

from __future__ import annotations

import importlib
import json
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import pytest

from processforge.providers.base import AbstractProvider
from processforge.providers.registry import _PROVIDERS
from processforge.types import MaterialDef

# Ensure the festim provider module is imported (self-registers on import)
importlib.import_module("processforge.providers.festim_provider")


# ---------------------------------------------------------------------------
# fake_festim — minimal stand-in for the `festim` package
# ---------------------------------------------------------------------------


class _FakeMesh1D:
    def __init__(self, vertices):
        self.vertices = vertices


class _FakeMaterial:
    def __init__(
        self,
        D_0,
        E_D,
        K_S_0=None,
        E_K_S=None,
        thermal_conductivity=None,
        density=None,
        heat_capacity=None,
        name=None,
        solubility_law="none",
    ):
        self.D_0 = D_0
        self.E_D = E_D
        self.K_S_0 = K_S_0
        self.E_K_S = E_K_S
        self.thermal_conductivity = thermal_conductivity
        self.density = density
        self.heat_capacity = heat_capacity
        self.name = name
        self.solubility_law = solubility_law


class _FakeSpecies:
    def __init__(self, name, mobile=True):
        self.name = name
        self.mobile = mobile

    def __str__(self):
        return self.name


class _FakeImplicitSpecies:
    def __init__(self, n, others=None, name=None):
        self.n = n
        self.others = others
        self.name = name


class _FakeVolumeSubdomain1D:
    def __init__(self, id, borders, material):
        self.id = id
        self.borders = borders
        self.material = material


class _FakeSurfaceSubdomain1D:
    def __init__(self, id, x):
        self.id = id
        self.x = x


class _FakeReaction:
    def __init__(self, reactant, k_0, E_k, volume, product=None, p_0=None, E_p=None):
        self.reactant = reactant
        self.k_0 = k_0
        self.E_k = E_k
        self.volume = volume
        self.product = product
        self.p_0 = p_0
        self.E_p = E_p


class _FakeParticleSource:
    def __init__(self, value, volume, species):
        self.value = value
        self.volume = volume
        self.species = species


class _FakeInitialConcentration:
    def __init__(self, value, volume, species):
        self.value = value
        self.volume = volume
        self.species = species


class _FakeFixedConcentrationBC:
    def __init__(self, subdomain, value, species):
        self.subdomain = subdomain
        self.value = value
        self.species = species


class _FakeParticleFluxBC:
    def __init__(self, subdomain, value, species):
        self.subdomain = subdomain
        self.value = value
        self.species = species


class _FakeSievertsBC:
    def __init__(self, subdomain, S_0, E_S, pressure, species):
        self.subdomain = subdomain
        self.S_0 = S_0
        self.E_S = E_S
        self.pressure = pressure
        self.species = species


class _FakeHenrysBC:
    def __init__(self, subdomain, H_0, E_H, pressure, species):
        self.subdomain = subdomain
        self.H_0 = H_0
        self.E_H = E_H
        self.pressure = pressure
        self.species = species


class _FakeSurfaceFlux:
    def __init__(self, field, surface, filename=None):
        self.field = field
        self.surface = surface
        self.filename = filename
        self.t = []
        self.data = []
        self.value = None

    @property
    def title(self):
        return f"{self.field.name} flux surface {self.surface.id}"


class _FakeTotalVolume:
    def __init__(self, field, volume, filename=None):
        self.field = field
        self.volume = volume
        self.filename = filename
        self.t = []
        self.data = []
        self.value = None

    @property
    def title(self):
        return f"Total {self.field.name} volume {self.volume.id}"


class _FakeProfile1DExport:
    def __init__(self, field, subdomain=None, times=None):
        self.field = field
        self.subdomain = subdomain
        self.times = times
        self.t = []
        self.data = []
        self.x = [0.0, 1.0]


class _FakeStepsize:
    def __init__(
        self,
        initial_value,
        growth_factor=None,
        cutback_factor=None,
        target_nb_iterations=None,
        max_stepsize=None,
        milestones=None,
        milestone_tolerance=1e-5,
    ):
        self.initial_value = initial_value
        self.growth_factor = growth_factor
        self.cutback_factor = cutback_factor
        self.target_nb_iterations = target_nb_iterations
        self.max_stepsize = max_stepsize
        self.milestones = milestones or []
        self.milestone_tolerance = milestone_tolerance


class _FakeSettings:
    def __init__(
        self,
        atol,
        rtol,
        max_iterations=30,
        transient=True,
        final_time=None,
        element_degree=1,
        stepsize=None,
        convergence_criterion="residual",
    ):
        self.atol = atol
        self.rtol = rtol
        self.max_iterations = max_iterations
        self.transient = transient
        self.final_time = final_time
        self.element_degree = element_degree
        self.stepsize = stepsize
        self.convergence_criterion = convergence_criterion


class _FakeHydrogenTransportProblem:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def initialise(self):
        pass

    def run(self):
        for export in getattr(self, "exports", []) or []:
            if hasattr(export, "title"):
                export.value = 12.34
            export.t = [0.0, 1.0, 2.0]
            export.data = [1.0, 2.0, 3.0]


def make_fake_festim() -> types.ModuleType:
    mod = types.ModuleType("festim")
    mod.Mesh1D = _FakeMesh1D
    mod.Material = _FakeMaterial
    mod.Species = _FakeSpecies
    mod.ImplicitSpecies = _FakeImplicitSpecies
    mod.VolumeSubdomain1D = _FakeVolumeSubdomain1D
    mod.SurfaceSubdomain1D = _FakeSurfaceSubdomain1D
    mod.Reaction = _FakeReaction
    mod.ParticleSource = _FakeParticleSource
    mod.InitialConcentration = _FakeInitialConcentration
    mod.FixedConcentrationBC = _FakeFixedConcentrationBC
    mod.ParticleFluxBC = _FakeParticleFluxBC
    mod.SievertsBC = _FakeSievertsBC
    mod.HenrysBC = _FakeHenrysBC
    mod.SurfaceFlux = _FakeSurfaceFlux
    mod.TotalVolume = _FakeTotalVolume
    mod.Profile1DExport = _FakeProfile1DExport
    mod.Stepsize = _FakeStepsize
    mod.Settings = _FakeSettings
    mod.HydrogenTransportProblem = _FakeHydrogenTransportProblem
    return mod


@pytest.fixture
def fake_festim(monkeypatch):
    fake = make_fake_festim()
    monkeypatch.setitem(sys.modules, "festim", fake)
    return fake


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mat(id_, name, **overrides):
    base = {"id": id_, "extra": {"D_0": 4.1e-7, "E_D": 0.39}}
    base.update(overrides)
    return MaterialDef.from_dict(base)


def _init_provider(tmp_path, materials, fake_festim=None):
    from processforge.providers.festim_provider import FestimProvider
    from processforge.types import FestimProviderConfig

    provider = FestimProvider()
    cfg = FestimProviderConfig(output_dir=str(tmp_path / "festim_run"))
    flowsheet = SimpleNamespace(materials=materials)
    provider.initialize(cfg, flowsheet)
    return provider


def _unit_config_from_flowsheet(path):
    from processforge.types import UnitConfig

    raw = json.loads(Path(path).read_text())
    unit_name, unit_dict = next(iter(raw["units"].items()))
    return UnitConfig.from_dict(unit_dict)


# ---------------------------------------------------------------------------
# 1. Registration
# ---------------------------------------------------------------------------


class TestFestimRegistration:
    """FESTIM provider must be registered in the global provider registry."""

    def test_festim_registered(self):
        assert "festim" in _PROVIDERS

    def test_festim_class_is_festim_provider(self):
        from processforge.providers.festim_provider import FestimProvider

        assert _PROVIDERS["festim"] is FestimProvider

    def test_festim_subclasses_abstract_provider(self):
        from processforge.providers.festim_provider import FestimProvider

        assert issubclass(FestimProvider, AbstractProvider)

    def test_all_abstract_methods_implemented(self):
        from processforge.providers.festim_provider import FestimProvider

        abstract_methods = set(AbstractProvider.__abstractmethods__)
        missing = abstract_methods - set(dir(FestimProvider))
        assert (
            not missing
        ), f"FestimProvider is missing abstract method implementations: {sorted(missing)}"


# ---------------------------------------------------------------------------
# 2. Material validation
# ---------------------------------------------------------------------------


class TestFestimMaterialValidation:
    """validate_material should enforce D_0/E_D validity."""

    def test_valid_material(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        mat_def = MaterialDef(id=1, extra={"D_0": 1e-7, "E_D": 0.2})
        unit_cfg = UnitConfig(
            type="SolverUnit", material=0, sim_type="hydrogen_transport_tds"
        )
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert errors == []

    def test_build_material_without_solubility_law(self, fake_festim):
        from processforge.providers.festim_provider import FestimBuildHelpers

        mat = FestimBuildHelpers.build_material(
            fake_festim,
            "tungsten",
            MaterialDef(id=1, extra={"D_0": 4.1e-7, "E_D": 0.39}),
        )
        assert mat.D_0 == 4.1e-7
        assert mat.E_D == 0.39
        # No solubility_law key was passed, so FESTIM applies its own NONE default
        # (the fake defaults to "none" and rejects None).
        assert mat.solubility_law == "none"

    def test_build_material_passes_solubility_law(self, fake_festim):
        from processforge.providers.festim_provider import FestimBuildHelpers

        mat = FestimBuildHelpers.build_material(
            fake_festim,
            "tungsten",
            MaterialDef(
                id=1, extra={"D_0": 4.1e-7, "E_D": 0.39, "solubility_law": "sievert"}
            ),
        )
        assert mat.solubility_law == "sievert"

    def test_build_material_warns_on_missing_optional_fields(self, fake_festim):
        from loguru import logger

        from processforge.providers.festim_provider import FestimBuildHelpers

        captured = []
        hid = logger.add(lambda m: captured.append(m), level="WARNING")
        try:
            mat = FestimBuildHelpers.build_material(
                fake_festim,
                "tungsten",
                MaterialDef(id=1, extra={"D_0": 4.1e-7, "E_D": 0.39}),
            )
        finally:
            logger.remove(hid)

        # Defaults applied for every missing optional field.
        assert mat.solubility_law == "none"
        assert mat.K_S_0 is None
        assert mat.E_K_S is None
        assert mat.thermal_conductivity is None
        assert mat.density is None
        assert mat.heat_capacity is None

        # A WARNING was emitted for each missing optional field.
        for field in (
            "solubility_law",
            "K_S_0",
            "E_K_S",
            "thermal_conductivity",
            "density",
            "heat_capacity",
        ):
            assert any(
                f"extra.{field}" in m and "defaulting to" in m for m in captured
            ), f"expected warning for missing {field}; got: {captured}"

    def test_missing_D_0(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        mat_def = MaterialDef(id=1, extra={"E_D": 0.2})
        unit_cfg = UnitConfig(
            type="SolverUnit", material=0, sim_type="hydrogen_transport_tds"
        )
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 1
        assert "D_0" in errors[0]

    def test_missing_E_D(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        mat_def = MaterialDef(id=1, extra={"D_0": 1e-7})
        unit_cfg = UnitConfig(
            type="SolverUnit", material=0, sim_type="hydrogen_transport_tds"
        )
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 1
        assert "E_D" in errors[0]

    def test_missing_both(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        mat_def = MaterialDef(id=1, extra={})
        unit_cfg = UnitConfig(
            type="SolverUnit", material=0, sim_type="hydrogen_transport_tds"
        )
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 2

    def test_missing_extra_entirely(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        mat_def = MaterialDef(id=1)
        unit_cfg = UnitConfig(
            type="SolverUnit", material=0, sim_type="hydrogen_transport_tds"
        )
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 2

    def test_non_positive_D_0(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        mat_def = MaterialDef(id=1, extra={"D_0": 0, "E_D": 0.2})
        unit_cfg = UnitConfig(
            type="SolverUnit", material=0, sim_type="hydrogen_transport_tds"
        )
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 1
        assert "positive" in errors[0]

    def test_negative_E_D(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        mat_def = MaterialDef(id=1, extra={"D_0": 1e-7, "E_D": -1})
        unit_cfg = UnitConfig(
            type="SolverUnit", material=0, sim_type="hydrogen_transport_tds"
        )
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 1
        assert "non-negative" in errors[0]

    def test_unpaired_solubility_constants(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        mat_def = MaterialDef(id=1, extra={"D_0": 1e-7, "E_D": 0.2, "K_S_0": 1e-5})
        unit_cfg = UnitConfig(
            type="SolverUnit", material=0, sim_type="hydrogen_transport_tds"
        )
        errors = FestimProvider.validate_material("steel", mat_def, unit_cfg)
        assert len(errors) == 1
        assert "together" in errors[0]


# ---------------------------------------------------------------------------
# 3. ProviderConfig
# ---------------------------------------------------------------------------


class TestFestimProviderConfig:
    """FestimProviderConfig should parse from dict correctly."""

    def test_default_config(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig()
        assert cfg.type == "festim"
        assert cfg.url is None
        assert cfg.output_dir == "outputs/festim"

    def test_from_dict_defaults(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict({"type": "festim"})
        assert cfg.url is None
        assert cfg.output_dir == "outputs/festim"

    def test_from_dict_with_url(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict(
            {
                "type": "festim",
                "url": "http://localhost:9002",
            }
        )
        assert cfg.url == "http://localhost:9002"

    def test_from_dict_custom_output_dir(self):
        from processforge.types import FestimProviderConfig

        cfg = FestimProviderConfig.from_dict(
            {
                "type": "festim",
                "output_dir": "/tmp/festim",
            }
        )
        assert cfg.output_dir == "/tmp/festim"

    def test_in_registry(self):
        from processforge.types import _PROVIDER_CONFIG_REGISTRY

        assert "festim" in _PROVIDER_CONFIG_REGISTRY

    def test_in_union(self):
        from processforge.types import FestimProviderConfig, provider_config_from_dict

        cfg = provider_config_from_dict({"type": "festim"})
        assert isinstance(cfg, FestimProviderConfig)


# ---------------------------------------------------------------------------
# 4. Catalog metadata
# ---------------------------------------------------------------------------


class TestFestimCatalog:
    """Festim catalog entry must reflect the Docker container architecture."""

    def test_has_docker_image(self):
        from processforge.providers.registry import _PROVIDER_CATALOG

        info = _PROVIDER_CATALOG["festim"]
        assert info["docker_image"] == "ghcr.io/urjanova/processforge-festim:latest"

    def test_has_default_port(self):
        from processforge.providers.registry import _PROVIDER_CATALOG

        info = _PROVIDER_CATALOG["festim"]
        assert info["default_port"] == 9002

    def test_is_containerized(self):
        from processforge.providers.registry import is_containerized

        assert is_containerized("festim") is True

    def test_no_optional_dep(self):
        from processforge.providers.registry import _PROVIDER_CATALOG

        info = _PROVIDER_CATALOG["festim"]
        assert info["optional_dep"] is None


# ---------------------------------------------------------------------------
# 5. Sim-type strategy registry
# ---------------------------------------------------------------------------


class TestFestimSimTypeRegistry:
    """FESTIM sim_type registry must map strings to strategy classes."""

    def test_get_registered_sim_types_callable(self):
        from processforge.providers.festim_provider import get_registered_sim_types

        result = get_registered_sim_types()
        assert isinstance(result, dict)

    def test_tds_strategy_registered(self):
        from processforge.providers.festim_provider import (
            FestimSimStrategy,
            get_registered_sim_types,
        )

        registry = get_registered_sim_types()
        assert "hydrogen_transport_tds" in registry
        assert issubclass(registry["hydrogen_transport_tds"], FestimSimStrategy)

    def test_legacy_hydrogen_transport_removed(self):
        from processforge.providers.festim_provider import get_registered_sim_types

        assert "hydrogen_transport" not in get_registered_sim_types()

    def test_register_new_sim_type(self):
        from processforge.providers.festim_provider import (
            FestimSimStrategy,
            _SIM_TYPE_REGISTRY,
            get_registered_sim_types,
            register_festim_sim_type,
        )

        class _MockStrategy(FestimSimStrategy):
            def build(self, festim, festim_model, materials_map, helpers):
                return {}

        register_festim_sim_type("custom_thing", _MockStrategy)
        try:
            registry = get_registered_sim_types()
            assert "custom_thing" in registry
            assert registry["custom_thing"] is _MockStrategy
        finally:
            del _SIM_TYPE_REGISTRY["custom_thing"]

    def test_returns_copy(self):
        from processforge.providers.festim_provider import get_registered_sim_types

        r1 = get_registered_sim_types()
        r2 = get_registered_sim_types()
        assert r1 == r2
        assert r1 is not r2

    def test_unknown_sim_type_raises(self):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        provider = FestimProvider()
        provider._initialized = True
        unit_cfg = UnitConfig(type="SolverUnit", sim_type="nonexistent_type")
        with pytest.raises(ValueError, match="unknown sim_type"):
            provider.run_simulation(unit_cfg, {})


# ---------------------------------------------------------------------------
# 6. Strategy build (via fake_festim)
# ---------------------------------------------------------------------------


class TestTDSTrappingStrategy:
    """The built-in TDS strategy builds a full HydrogenTransportProblem."""

    def _build(self, fake_festim):
        from processforge.providers.festim_provider import (
            FestimBuildHelpers,
            _TDSTrappingStrategy,
        )
        from processforge.schemas.festim.festim_model import FestimModel

        sc = _unit_config_from_flowsheet(
            "flowsheets/festim/tds_tungsten.json"
        ).solver_config
        model = FestimModel.model_validate(sc)
        materials_map = {
            "tungsten": fake_festim.Material(D_0=4.1e-7, E_D=0.39, name="tungsten")
        }
        helpers = FestimBuildHelpers()
        return _TDSTrappingStrategy().build(fake_festim, model, materials_map, helpers)

    def test_builds_mesh_from_linspace_segments(self, fake_festim):
        kwargs = self._build(fake_festim)
        mesh = kwargs["mesh"]
        assert isinstance(mesh, _FakeMesh1D)
        assert mesh.vertices.shape == (700,)
        assert mesh.vertices[0] == pytest.approx(0.0)
        assert mesh.vertices[-1] == pytest.approx(20e-6)

    def test_builds_species_and_implicit(self, fake_festim):
        kwargs = self._build(fake_festim)
        species = kwargs["species"]
        assert [s.name for s in species] == ["H", "trapped_H1", "trapped_H2"]
        assert species[0].mobile is True
        assert species[1].mobile is False

    def test_builds_subdomains(self, fake_festim):
        kwargs = self._build(fake_festim)
        subdomains = kwargs["subdomains"]
        volumes = [s for s in subdomains if isinstance(s, _FakeVolumeSubdomain1D)]
        surfaces = [s for s in subdomains if isinstance(s, _FakeSurfaceSubdomain1D)]
        assert [v.id for v in volumes] == [1]
        assert volumes[0].material.name == "tungsten"
        assert [s.id for s in surfaces] == [1, 2]

    def test_builds_reactions_with_implicit_reactants(self, fake_festim):
        kwargs = self._build(fake_festim)
        reactions = kwargs["reactions"]
        assert len(reactions) == 2
        r1 = reactions[0]
        assert r1.k_0 == pytest.approx(53983071.22305305)
        assert r1.E_k == 0.39
        assert r1.p_0 == 1e13
        assert r1.E_p == 0.87
        assert any(isinstance(s, _FakeImplicitSpecies) for s in r1.reactant)

    def test_builds_gaussian_source(self, fake_festim):
        kwargs = self._build(fake_festim)
        sources = kwargs["sources"]
        assert len(sources) == 1
        assert callable(sources[0].value)
        assert sources[0].species.name == "H"

    def test_builds_boundary_conditions(self, fake_festim):
        kwargs = self._build(fake_festim)
        bcs = kwargs["boundary_conditions"]
        assert len(bcs) == 2
        assert all(isinstance(b, _FakeFixedConcentrationBC) for b in bcs)
        assert bcs[0].value == 0

    def test_builds_temperature_ramp(self, fake_festim):
        kwargs = self._build(fake_festim)
        temp = kwargs["temperature"]
        assert callable(temp)
        assert temp(0) == 300
        assert temp(449) == 300
        assert temp(451) == pytest.approx(308)
        assert temp(500) == pytest.approx(700)

    def test_builds_exports(self, fake_festim):
        kwargs = self._build(fake_festim)
        exports = kwargs["exports"]
        assert len(exports) == 8
        flux = [e for e in exports if isinstance(e, _FakeSurfaceFlux)]
        totals = [e for e in exports if isinstance(e, _FakeTotalVolume)]
        profiles = [e for e in exports if isinstance(e, _FakeProfile1DExport)]
        assert len(flux) == 2
        assert len(totals) == 3
        assert len(profiles) == 3
        assert flux[0].filename == "flux_left.csv"

    def test_builds_settings_with_adaptive_stepsize(self, fake_festim):
        kwargs = self._build(fake_festim)
        settings = kwargs["settings"]
        assert isinstance(settings, _FakeSettings)
        assert settings.atol == 1e-10
        assert settings.final_time == 500
        assert isinstance(settings.stepsize, _FakeStepsize)
        assert settings.stepsize.initial_value == 0.5
        assert settings.stepsize.milestones == [400, 450, 500]
        # after_t max_stepsize becomes a callable
        assert callable(settings.stepsize.max_stepsize)
        assert settings.stepsize.max_stepsize(100) is None
        assert settings.stepsize.max_stepsize(451) == 0.5

    def test_missing_material_raises(self, fake_festim):
        from processforge.providers.festim_provider import (
            FestimBuildHelpers,
            _TDSTrappingStrategy,
        )
        from processforge.schemas.festim.festim_model import FestimModel

        sc = _unit_config_from_flowsheet(
            "flowsheets/festim/tds_tungsten.json"
        ).solver_config
        model = FestimModel.model_validate(sc)
        helpers = FestimBuildHelpers()
        with pytest.raises(ValueError, match="not defined in the flowsheet materials"):
            _TDSTrappingStrategy().build(fake_festim, model, {}, helpers)


# ---------------------------------------------------------------------------
# 7. run_simulation behaviour
# ---------------------------------------------------------------------------


class TestFestimRunSimulation:
    """run_simulation runs the problem in-process and extracts results."""

    def test_runs_tds_flowsheet(self, tmp_path, fake_festim):
        provider = _init_provider(
            tmp_path,
            {"tungsten": _mat(1, "tungsten")},
            fake_festim=fake_festim,
        )
        unit_cfg = _unit_config_from_flowsheet("flowsheets/festim/tds_tungsten.json")

        result = provider.run_simulation(unit_cfg, {})

        assert result.status == "completed"
        assert result.sim_type == "hydrogen_transport_tds"
        assert "H flux surface 1" in result.scalars
        assert result.scalars["H flux surface 1"] == pytest.approx(12.34)
        assert "Total trapped_H1 volume 1" in result.scalars
        assert result.metadata["run_dir"] == str((tmp_path / "festim_run").resolve())
        assert "H flux surface 1" in result.metadata["series"]

    def test_exports_written_into_run_dir(self, tmp_path, fake_festim):
        provider = _init_provider(
            tmp_path,
            {"tungsten": _mat(1, "tungsten")},
            fake_festim=fake_festim,
        )
        unit_cfg = _unit_config_from_flowsheet("flowsheets/festim/tds_tungsten.json")

        provider.run_simulation(unit_cfg, {})

        assert (tmp_path / "festim_run").is_dir()
        assert (tmp_path / "festim_run" / "flux_left.csv").is_dir() is False

    def test_failure_returns_failed_status(self, tmp_path, fake_festim, monkeypatch):
        provider = _init_provider(
            tmp_path,
            {"tungsten": _mat(1, "tungsten")},
            fake_festim=fake_festim,
        )
        unit_cfg = _unit_config_from_flowsheet("flowsheets/festim/tds_tungsten.json")

        def boom(self):
            raise RuntimeError("SNES diverged")

        monkeypatch.setattr(_FakeHydrogenTransportProblem, "run", boom)

        result = provider.run_simulation(unit_cfg, {})
        assert result.status == "failed"
        assert "SNES diverged" in result.metadata["error"]
        assert result.metadata["run_dir"] == str((tmp_path / "festim_run").resolve())

    def test_run_before_init_raises(self, fake_festim):
        from processforge.providers.festim_provider import FestimProvider
        from processforge.types import UnitConfig

        provider = FestimProvider()
        unit_cfg = UnitConfig(type="SolverUnit", sim_type="hydrogen_transport_tds")
        with pytest.raises(RuntimeError, match="not been initialized"):
            provider.run_simulation(unit_cfg, {})

    def test_teardown_resets_state(self):
        from processforge.providers.festim_provider import FestimProvider

        provider = FestimProvider()
        provider._initialized = True
        provider.teardown()
        assert provider._initialized is False

    def test_compute_unit_returns_none(self):
        from processforge.providers.festim_provider import FestimProvider

        provider = FestimProvider()
        assert provider.compute_unit("SolverUnit", {}, {}) is None

    def test_get_thermo_properties_raises(self):
        from processforge.providers.festim_provider import FestimProvider

        provider = FestimProvider()
        with pytest.raises(NotImplementedError):
            provider.get_thermo_properties({})
