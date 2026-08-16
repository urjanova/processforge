"""Tests for the FESTIM solver_config Pydantic schema (festim_model.py).

Covers the happy path against the bundled flowsheets plus the cross-field
validation rules: mesh definition, species/implicit-species references,
subdomain id uniqueness, reaction rate pairing, transient consistency,
stepsize adaptivity, and ``extra="forbid"`` strictness.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from processforge.schemas.festim.festim_model import (
    FestimModel,
    GaussianProfile,
    MaterialConfig,
    RampProfile,
    AfterTProfile,
)


def _load_solver_config(flowsheet: str) -> dict:
    raw = json.loads(
        (Path(__file__).parent.parent / "flowsheets" / "festim" / flowsheet).read_text()
    )
    return copy.deepcopy(next(iter(raw["units"].values()))["solver_config"])


@pytest.fixture
def tds_config() -> dict:
    return _load_solver_config("tds_tungsten.json")


@pytest.fixture
def hydrogen_config() -> dict:
    return _load_solver_config("hydrogen_transport_1d.json")


# ---------------------------------------------------------------------------
# Valid configs
# ---------------------------------------------------------------------------


class TestValidConfigs:
    def test_tds_flowsheet_validates(self, tds_config):
        model = FestimModel.model_validate(tds_config)
        assert model.final_time == 500
        assert len(model.species) == 3
        assert len(model.implicit_species) == 2
        assert len(model.reactions) == 2
        assert len(model.exports) == 8

    def test_hydrogen_flowsheet_validates(self, hydrogen_config):
        model = FestimModel.model_validate(hydrogen_config)
        # bare float temperature is coerced to constant
        assert model.temperature == 300
        assert model.stepsize.milestones == [0.05, 0.1, 0.2, 0.5, 1]

    def test_gaussian_profile(self):
        p = GaussianProfile(
            amplitude=4.15e-5, center=4.5e-9, width=2.5e-9, active_until=400
        )
        assert p.type == "gaussian"
        assert p.active_until == 400

    def test_gaussian_profile_width_positive(self):
        with pytest.raises(ValidationError, match="greater than 0"):
            GaussianProfile(amplitude=1, center=0, width=0)

    def test_material_config_solubility_pairing(self):
        MaterialConfig(D_0=1e-7, E_D=0.2, K_S_0=1e-5, E_K_S=0.5)
        with pytest.raises(ValidationError, match="provided together"):
            MaterialConfig(D_0=1e-7, E_D=0.2, K_S_0=1e-5)

    def test_material_config_d0_positive(self):
        with pytest.raises(ValidationError, match="greater than 0"):
            MaterialConfig(D_0=0, E_D=0.2)


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------


class TestMeshValidation:
    def test_neither_vertices_nor_segments(self, tds_config):
        tds_config["mesh"] = {}
        with pytest.raises(
            ValidationError, match="exactly one of 'vertices' or 'segments'"
        ):
            FestimModel.model_validate(tds_config)

    def test_both_vertices_and_segments(self, tds_config):
        tds_config["mesh"] = {"vertices": [0, 1], "segments": []}
        with pytest.raises(
            ValidationError, match="exactly one of 'vertices' or 'segments'"
        ):
            FestimModel.model_validate(tds_config)

    def test_vertices_too_short(self, tds_config):
        tds_config["mesh"] = {"vertices": [0.0]}
        with pytest.raises(ValidationError, match="at least 2 points"):
            FestimModel.model_validate(tds_config)

    def test_segment_num_too_small(self, tds_config):
        tds_config["mesh"]["segments"][0]["num"] = 1
        with pytest.raises(ValidationError, match="greater than 1"):
            FestimModel.model_validate(tds_config)


# ---------------------------------------------------------------------------
# Species
# ---------------------------------------------------------------------------


class TestSpeciesValidation:
    def test_duplicate_species(self, tds_config):
        tds_config["species"].append({"name": "H", "mobile": True})
        with pytest.raises(ValidationError, match="names must be unique"):
            FestimModel.model_validate(tds_config)

    def test_no_mobile_species(self, tds_config):
        tds_config["species"] = [
            {"name": "trapped_H1", "mobile": False},
            {"name": "trapped_H2", "mobile": False},
        ]
        with pytest.raises(ValidationError, match="at least one mobile"):
            FestimModel.model_validate(tds_config)

    def test_empty_species(self, tds_config):
        tds_config["species"] = []
        with pytest.raises(ValidationError, match="at least one species"):
            FestimModel.model_validate(tds_config)

    def test_implicit_name_collides_with_species(self, tds_config):
        tds_config["implicit_species"][0]["name"] = "H"
        with pytest.raises(ValidationError, match="must not collide"):
            FestimModel.model_validate(tds_config)


# ---------------------------------------------------------------------------
# Subdomains
# ---------------------------------------------------------------------------


class TestSubdomainValidation:
    def test_zero_volume_id(self, tds_config):
        tds_config["subdomains"]["volume"][0]["id"] = 0
        with pytest.raises(ValidationError, match="cannot be 0"):
            FestimModel.model_validate(tds_config)

    def test_duplicate_volume_ids(self, tds_config):
        tds_config["subdomains"]["volume"].append(
            {"id": 1, "borders": [0.0, 1e-6], "material": "tungsten"}
        )
        with pytest.raises(
            ValidationError, match="volume subdomain ids must be unique"
        ):
            FestimModel.model_validate(tds_config)

    def test_duplicate_surface_ids(self, tds_config):
        tds_config["subdomains"]["surface"][1]["id"] = 1
        with pytest.raises(
            ValidationError, match="surface subdomain ids must be unique"
        ):
            FestimModel.model_validate(tds_config)


# ---------------------------------------------------------------------------
# References (sources, BCs, reactions, exports)
# ---------------------------------------------------------------------------


class TestReferenceValidation:
    def test_source_unknown_species(self, tds_config):
        tds_config["sources"][0]["species"] = "nope"
        with pytest.raises(ValidationError, match="undefined species"):
            FestimModel.model_validate(tds_config)

    def test_source_unknown_volume(self, tds_config):
        tds_config["sources"][0]["volume"] = 99
        with pytest.raises(ValidationError, match="not a defined volume"):
            FestimModel.model_validate(tds_config)

    def test_bc_unknown_species(self, tds_config):
        tds_config["boundary_conditions"][0]["species"] = "nope"
        with pytest.raises(ValidationError, match="undefined"):
            FestimModel.model_validate(tds_config)

    def test_bc_unknown_surface(self, tds_config):
        tds_config["boundary_conditions"][0]["subdomain"] = 99
        with pytest.raises(ValidationError, match="not a defined surface"):
            FestimModel.model_validate(tds_config)

    def test_reaction_unknown_species(self, tds_config):
        tds_config["reactions"][0]["reactant"] = ["H", "nope"]
        with pytest.raises(ValidationError, match="neither.*species nor an implicit"):
            FestimModel.model_validate(tds_config)

    def test_reaction_unknown_volume(self, tds_config):
        tds_config["reactions"][0]["volume"] = 99
        with pytest.raises(ValidationError, match="not a defined volume"):
            FestimModel.model_validate(tds_config)

    def test_export_unknown_species(self, tds_config):
        tds_config["exports"][0]["field"] = "nope"
        with pytest.raises(ValidationError, match="undefined"):
            FestimModel.model_validate(tds_config)

    def test_export_unknown_volume(self, tds_config):
        tds_config["exports"][5]["volume"] = 99
        with pytest.raises(ValidationError, match="not a defined volume"):
            FestimModel.model_validate(tds_config)


# ---------------------------------------------------------------------------
# Reactions
# ---------------------------------------------------------------------------


class TestReactionValidation:
    def test_product_without_rates(self, tds_config):
        r = tds_config["reactions"][0]
        r.pop("p_0")
        r.pop("E_p")
        with pytest.raises(
            ValidationError, match="required when 'product' is non-empty"
        ):
            FestimModel.model_validate(tds_config)

    def test_rates_without_product(self, tds_config):
        tds_config["reactions"][0]["product"] = []
        with pytest.raises(ValidationError, match="omitted when there is no product"):
            FestimModel.model_validate(tds_config)

    def test_irreversible_reaction_valid(self, tds_config):
        tds_config["reactions"] = [
            {
                "type": "arrhenius",
                "reactant": ["H", "empty_trap1"],
                "product": [],
                "k_0": 1e4,
                "E_k": 0.39,
                "volume": 1,
            }
        ]
        FestimModel.model_validate(tds_config)

    def test_negative_k0(self, tds_config):
        tds_config["reactions"][0]["k_0"] = -1
        with pytest.raises(ValidationError, match="greater than 0"):
            FestimModel.model_validate(tds_config)


# ---------------------------------------------------------------------------
# Transient consistency
# ---------------------------------------------------------------------------


class TestTransientConsistency:
    def test_transient_without_final_time(self, tds_config):
        del tds_config["final_time"]
        with pytest.raises(ValidationError, match="'final_time' is required"):
            FestimModel.model_validate(tds_config)

    def test_transient_without_stepsize(self, tds_config):
        del tds_config["stepsize"]
        with pytest.raises(ValidationError, match="'stepsize' is required"):
            FestimModel.model_validate(tds_config)

    def test_steady_state_ok(self, tds_config):
        tds_config["transient"] = False
        del tds_config["final_time"]
        del tds_config["stepsize"]
        FestimModel.model_validate(tds_config)


# ---------------------------------------------------------------------------
# Stepsize
# ---------------------------------------------------------------------------


class TestStepsizeValidation:
    def test_milestones_require_adaptive_factors(self, tds_config):
        tds_config["stepsize"].pop("growth_factor")
        tds_config["stepsize"].pop("cutback_factor")
        with pytest.raises(
            ValidationError, match="require 'growth_factor' and 'cutback_factor'"
        ):
            FestimModel.model_validate(tds_config)

    def test_growth_factor_below_one(self, tds_config):
        tds_config["stepsize"]["growth_factor"] = 0.9
        with pytest.raises(ValidationError, match="'growth_factor' must be >= 1"):
            FestimModel.model_validate(tds_config)

    def test_max_stepsize_below_initial(self, tds_config):
        tds_config["stepsize"]["max_stepsize"] = 0.1
        with pytest.raises(
            ValidationError, match="must not be less than 'initial_value'"
        ):
            FestimModel.model_validate(tds_config)

    def test_after_t_max_stepsize_valid(self):
        AfterTProfile(after_t=450, value=0.5)

    def test_ramp_profile(self):
        p = RampProfile(T0=300, start_t=450, ramp_rate=8)
        assert p.type == "ramp"


# ---------------------------------------------------------------------------
# Strictness
# ---------------------------------------------------------------------------


class TestStrictness:
    def test_unknown_top_level_field_rejected(self, tds_config):
        tds_config["bogus"] = True
        with pytest.raises(ValidationError):
            FestimModel.model_validate(tds_config)

    def test_unknown_bc_type_rejected(self, tds_config):
        tds_config["boundary_conditions"][0]["type"] = "magic"
        with pytest.raises(ValidationError):
            FestimModel.model_validate(tds_config)

    def test_unknown_export_type_rejected(self, tds_config):
        tds_config["exports"][0]["type"] = "vtx_species"
        with pytest.raises(ValidationError):
            FestimModel.model_validate(tds_config)

    def test_unknown_reaction_type_rejected(self, tds_config):
        tds_config["reactions"][0]["type"] = "custom"
        with pytest.raises(ValidationError):
            FestimModel.model_validate(tds_config)
