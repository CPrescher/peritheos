import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from scripts.audit_fortes_2019_fcc_pb import (
    KUZNETSOV_FCC,
    KUZNETSOV_HCP,
    PUBLISHED,
    audit,
    fortes_pressure,
    kuznetsov_hcp_pressure,
    kuznetsov_hcp_temperature_parameters,
    kuznetsov_pressure,
    kuznetsov_temperature_parameters,
    temperature_parameters,
)

ROOT = Path(__file__).parents[1]
RESOURCE = "lead-fcc-kuznetsov-2002-table1-transition-pvt.csv"


def _record():
    document = get_material_document("lead_fcc")
    return next(
        row
        for row in document["eos_records"]
        if row["identifier"] == "lead_fcc_fortes_2019_bm4_1"
    )


def _kuznetsov_record():
    document = get_material_document("lead_fcc")
    return next(
        row
        for row in document["eos_records"]
        if row["identifier"] == "lead_fcc_kuznetsov_2002_bm3_2"
    )


def test_fortes_300k_slice_matches_public_bm4():
    document = get_material_document("lead_fcc")
    record = _record()
    executable = Material.from_eosmat(document).get_eos_record(record["identifier"])
    volumes = np.asarray([118.0, 110.0, 102.0])
    assert executable.pressure(volumes) == pytest.approx(
        fortes_pressure(volumes, 300.0), abs=1.0e-12
    )
    assert record["eos"]["parameters"] == {
        "V0": PUBLISHED["V0_300"],
        "K0": PUBLISHED["K0_300"],
        "K0_prime": PUBLISHED["K0_prime_300"],
        "K0_double_prime": PUBLISHED["K0_double_prime"],
    }


def test_fortes_temperature_polynomials_have_300k_constant_terms():
    v0, k0, k0_prime, k0_double_prime = temperature_parameters(300.0)
    assert float(v0) == PUBLISHED["V0_300"]
    assert float(k0) == PUBLISHED["K0_300"]
    assert float(k0_prime) == PUBLISHED["K0_prime_300"]
    assert float(k0_double_prime) == PUBLISHED["K0_double_prime"]


def test_kuznetsov_table1_partial_pvt_fixture_is_exact_and_linked():
    path = ROOT / "peritheos" / "data" / "datasets" / RESOURCE
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert [(row["temperature_k"], row["pressure_gpa"]) for row in rows] == [
        ("296", "13.1"),
        ("402", "13.9"),
        ("469", "12.6"),
    ]
    assert [float(row["volume_a3_conventional_cell"]) for row in rows] == [
        100.48,
        100.56,
        101.04,
    ]
    document = get_material_document("lead_fcc")
    dataset = next(
        row
        for row in document["datasets"]
        if row["identifier"] == "lead_fcc_kuznetsov_2002_table1_transition_pvt"
    )
    assert (
        dataset["resource"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    )
    assert dataset["used_by_eos_records"] == [
        "lead_fcc_fortes_2019_bm4_1",
        "lead_fcc_kuznetsov_2002_bm3_2",
    ]


def test_kuznetsov_published_predecessor_model_is_separately_executable():
    v0, k0, k0_prime = kuznetsov_temperature_parameters(296.0)
    assert float(v0) == KUZNETSOV_FCC["V0_296_atomic"]
    assert float(k0) == KUZNETSOV_FCC["K0_296"]
    assert float(k0_prime) == KUZNETSOV_FCC["K0_prime_296"]
    result = audit()["kuznetsov_published_predecessor_model"]
    assert "not a reconstruction" in result["role"]
    assert result["transition_anchor_pressure_rmse_gpa"] == pytest.approx(
        0.456, abs=0.001
    )
    assert kuznetsov_pressure(30.307, 296.0) == pytest.approx(0.0, abs=1.0e-12)


def test_kuznetsov_296k_slice_is_a_published_fcc_eos_record():
    document = get_material_document("lead_fcc")
    record = _kuznetsov_record()
    executable = Material.from_eosmat(document).get_eos_record(record["identifier"])
    volumes = np.asarray([118.0, 110.0, 102.0])
    assert executable.pressure(volumes) == pytest.approx(
        kuznetsov_pressure(volumes / 4.0, 296.0), abs=1.0e-12
    )
    assert record["record_kind"] == "published"
    assert record["equation_kind"] == "isothermal"
    assert record["temperature_ref"] == 296.0
    assert record["eos"]["parameters"] == {
        "V0": 4.0 * KUZNETSOV_FCC["V0_296_atomic"],
        "K0": KUZNETSOV_FCC["K0_296"],
        "K0_prime": KUZNETSOV_FCC["K0_prime_296"],
    }
    assert "ambient-pressure phase" in document["phase"]
    assert (
        "transition onset near 12 GPa"
        in record["scientific_validation"]["primary_source_check"]["finding"]
    )


def test_kuznetsov_hcp_296k_slice_is_verified_before_cataloging():
    v0, k0, k0_prime = kuznetsov_hcp_temperature_parameters(296.0)
    assert float(v0) == KUZNETSOV_HCP["V0_296_atomic"]
    assert float(k0) == KUZNETSOV_HCP["K0_296"]
    assert float(k0_prime) == KUZNETSOV_HCP["K0_prime_296"]
    assert kuznetsov_hcp_pressure(29.908, 296.0) == pytest.approx(0.0, abs=1.0e-12)

    result = audit()["kuznetsov_hcp_published_model"]
    anchors = result["table1_transition_anchors"]
    markers = result["figure3_reduced_room_temperature_markers"]
    assert anchors["observations"] == 3
    assert anchors["pressure_rmse_gpa"] == pytest.approx(0.7538, abs=0.001)
    assert markers["observations"] == 21
    assert markers["pressure_rmse_gpa"] == pytest.approx(0.6418, abs=0.001)
    assert markers["solid_room_temperature_observations"] == 7
    assert markers["solid_room_temperature_pressure_rmse_gpa"] == pytest.approx(
        0.2662, abs=0.001
    )


def test_kuznetsov_hcp_296k_slice_is_a_published_eos_record():
    document = get_material_document("lead_hcp")
    record = next(
        row
        for row in document["eos_records"]
        if row["identifier"] == "lead_hcp_kuznetsov_2002_bm3_3"
    )
    executable = Material.from_eosmat(document).get_eos_record(record["identifier"])
    cell_volumes = np.asarray([56.0, 52.0, 48.0])
    assert executable.pressure(cell_volumes) == pytest.approx(
        kuznetsov_hcp_pressure(cell_volumes / 2.0, 296.0), abs=1.0e-12
    )
    assert record["record_kind"] == "published"
    assert record["equation_kind"] == "isothermal"
    assert record["temperature_ref"] == 296.0
    assert record["eos"]["parameters"] == {
        "V0": 2.0 * KUZNETSOV_HCP["V0_296_atomic"],
        "K0": KUZNETSOV_HCP["K0_296"],
        "K0_prime": KUZNETSOV_HCP["K0_prime_296"],
    }
    dataset = next(
        item
        for item in document["datasets"]
        if record["identifier"] in item.get("used_by_eos_records", [])
    )
    assert dataset["reference"]["doi"] == record["reference"]["doi"]


def test_reported_fit_limitations_are_machine_checkable():
    result = audit()
    assert result["outcome"] == "not_independently_refittable"
    assert result["published_surface_executable"] is True
    fit = result["reported_fit_attempt"]
    assert fit["is_reported_fortes_fit"] is False
    assert fit["observations"] == 3
    assert fit["published_parameters_pressure_rmse_gpa"] == pytest.approx(
        0.8146, abs=0.001
    )
    assert fit["diagnostic_fit_pressure_rmse_gpa"] == pytest.approx(0.3844, abs=0.001)
    assert fit["jacobian_condition_number"] > 1.0e8
    assert result["room_temperature_identifiability"]["d_pressure_d_e_gpa_k"] == [
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    assert len(result["missing_for_independent_refit"]) == 5


def test_material_provenance_does_not_overclaim_the_partial_fit():
    record = _record()
    assert record["reference"]["doi"] == "10.5286/raltr.2019002"
    assert record["fixed_parameters"] == ["V0", "K0"]
    provenance = record["published_fit_provenance"]
    assert provenance["refined_parameters"] == [
        "K0_prime_300",
        "e",
        "K0_double_prime",
    ]
    assert provenance["weights"] == "not reported"
    check = record["scientific_validation"]["primary_data_check"]
    assert check["status"] == "combined_primary_data_partially_bundled"
    assert check["dataset_identifiers"] == [
        "lead_fcc_kuznetsov_2002_table1_transition_pvt"
    ]
    assert "cannot be independently reproduced" in check["finding"]
