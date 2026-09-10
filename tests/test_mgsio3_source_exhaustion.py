"""Regression checks for the MgSiO3-family source-exhaustion audit."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from peritheos.eosmat import validate_eosmat_document
from scripts.reproduce_akins_2004_mgsio3_liquid import (
    candidate_hugoniot_pressure_gpa,
)
from scripts.reproduce_akins_2004_mgsio3_liquid import (
    reproduce as reproduce_akins,
)
from scripts.reproduce_mgsio3_source_exhaustion import reproduce
from scripts.validate_primary_eos_refits import validate_all

ROOT = Path(__file__).resolve().parents[1]
AUDITED_DOIS = {
    "10.1126/science.251.4992.410",
    "10.1029/96gl03027",
    "10.1098/rsta.1996.0053",
    "10.1088/0256-307x/17/3/022",
    "10.2138/am-2000-2-309",
    "10.2138/am-2000-2-310",
    "10.2465/jmps.95.236",
    "10.1029/2003gl018762",
    "10.1029/2004gl020237",
    "10.1360/cjcp2006.19(4).311.4",
    "10.1088/1674-0068/20/05/547-551",
    "10.1111/j.1365-246x.1975.tb06463.x",
    "10.1016/j.pnsc.2009.09.002",
    "10.1016/j.rgg.2015.01.011",
}


def _load(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def test_every_assigned_litcurate_row_has_a_disposition() -> None:
    candidates = _load("docs/data/litcurate-eos-candidates.json")["records"]
    selected = [
        row
        for row in candidates
        if row["publication"]["doi"].casefold() in AUDITED_DOIS
    ]
    audit = (
        ROOT / "docs/literature-reproductions/mgsio3-source-exhaustion.md"
    ).read_text(encoding="utf-8")

    assert len(selected) == 44
    assert {row["publication"]["doi"].casefold() for row in selected} == AUDITED_DOIS
    assert all(f"`{row['identifier']}`" in audit for row in selected)


@pytest.mark.parametrize(
    "material_file",
    [
        "bridgmanite.eosmat",
        "akimotoite.eosmat",
        "mgsio3_post_perovskite.eosmat",
        "mgsio3_liquid.eosmat",
    ],
)
def test_changed_material_cards_validate(material_file: str) -> None:
    validate_eosmat_document(_load(f"peritheos/data/materials/{material_file}"))


def test_source_equations_reproduce_serialized_models() -> None:
    result = reproduce()
    residuals = result["maximum_absolute_pressure_residual_gpa"]
    assert max(residuals.values()) < 1.0e-10
    assert result["hamahata_v0_a3_z4"] == pytest.approx(165.45561819988183)
    assert result["akins_v0_a3_per_formula_unit"] == pytest.approx(45.29797155879917)


def test_akins_candidate_hugoniot_reconstruction_is_similar_but_not_identified() -> (
    None
):
    result = reproduce_akins()
    published = result["published_candidate"]
    reconstructed = result["bounded_unweighted_reconstruction"]

    assert result["shots"] == [318, 322, 319]
    assert published["calculated_pressure_gpa"] == pytest.approx(
        [177.572, 202.954, 195.798], abs=0.001
    )
    assert reconstructed["parameters"]["rho0_g_cm3"] == pytest.approx(
        3.68912, rel=2.0e-5
    )
    assert reconstructed["parameters"]["K0S_gpa"] == pytest.approx(125.463, rel=2.0e-5)
    assert reconstructed["parameters"]["V0_a3_per_formula_unit"] == pytest.approx(
        45.18599, rel=2.0e-5
    )
    assert reconstructed["pressure_rmse_gpa"] == pytest.approx(6.79454, rel=2.0e-5)
    assert reconstructed["degrees_of_freedom"] == 1
    assert reconstructed["parameter_correlation"] > 0.999
    assert result["identifiability"]["K0S_prime_refittable"] is False

    weighted = result["weighted_sensitivity_check"]["parameters"]
    assert abs(weighted["rho0_g_cm3"] - reconstructed["parameters"]["rho0_g_cm3"]) > 0.2
    assert abs(weighted["K0S_gpa"] - reconstructed["parameters"]["K0S_gpa"]) > 20


def test_akins_hugoniot_is_not_a_direct_bm3_pressure_fit() -> None:
    result = reproduce_akins()
    semantics = result["data_semantics"]
    assert "shock velocity" in semantics["directly_measured"]
    assert "pressure" in semantics["impedance_matched"]
    assert semantics["calculated_model_output"] == "candidate Hugoniot pressure"

    pressure = candidate_hugoniot_pressure_gpa([3.199], [5.68], 3.68, 125.0)
    assert pressure[0] == pytest.approx(177.572, abs=0.001)


def test_akins_primary_refit_ledger_uses_bounded_reconstruction() -> None:
    outcome = next(
        item
        for item in validate_all()["records"]
        if item["record_identifier"] == "mgsio3_liquid_akins_2004_adiabatic_bm3"
    )
    assert outcome["status"] == "similar"
    assert outcome["fit_kind"] == "mie_gruneisen_hugoniot_bounded_reconstruction"
    assert outcome["observations"] == 3
    assert outcome["free_parameters"] == ["V0", "K0"]
    assert outcome["fixed_parameters"] == []
    assert outcome["reconstruction_fixed_parameters"] == [
        "K0_prime",
        "gamma0",
        "q",
        "transition_energy",
    ]
