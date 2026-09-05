import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document, validate_eosmat_document
from scripts.reproduce_lee_wan_2008_mgo import bm3_pressure, reproduce

ROOT = Path(__file__).parents[1]
DOI = "10.1103/physrevb.78.224103"


def _records() -> list[dict]:
    document = get_material_document("mgo")
    return [
        record
        for record in document["eos_records"]
        if record["reference"].get("doi", "").lower() == DOI
    ]


def test_lee_wan_preserves_both_source_computed_branches():
    records = _records()
    assert [record["identifier"] for record in records] == [
        "mgo_lee_wan_2008_lda_static_bm3",
        "mgo_lee_wan_2008_gga_static_bm3",
    ]
    assert [tuple(record["eos"]["parameters"].values()) for record in records] == [
        (pytest.approx(70.059879275145, abs=1e-12), 177.486, 4.026),
        (pytest.approx(74.952382755288, abs=1e-12), 149.32, 4.08),
    ]
    assert all("default_for" not in record for record in records)
    assert all(record["temperature_ref"] == 0.0 for record in records)


def test_lee_wan_document_validates_and_curves_are_executable():
    document = get_material_document("mgo")
    validate_eosmat_document(document)
    material = Material.from_eosmat(document)
    for record in _records():
        parameters = record["eos"]["parameters"]
        executable = material.get_eos_record(record["identifier"])
        volumes = parameters["V0"] * np.array([1.0, 0.9, 0.8])
        assert np.asarray(executable.pressure(volumes)) == pytest.approx(
            bm3_pressure(
                volumes,
                parameters["V0"],
                parameters["K0"],
                parameters["K0_prime"],
            )
        )


def test_lee_wan_reproduction_checks_conversion_and_pressure_difference():
    result = reproduce()
    assert result["converted_parameters_z4"]["LDA"]["V0_a3"] == pytest.approx(
        70.059879275145, abs=1e-12
    )
    assert result["converted_parameters_z4"]["GGA"]["V0_a3"] == pytest.approx(
        74.952382755288, abs=1e-12
    )
    check = result["figure_curve_check"]
    assert len(check["volumes_bohr3"]) == 5
    assert check["maximum_absolute_difference_from_polynomial_gpa"] < 0.15
    assert all(value > 0 for value in check["gga_minus_lda_gpa"])


def test_lee_wan_audit_disposes_both_same_doi_candidates_once():
    audit = (
        ROOT
        / "docs"
        / "literature-reproductions"
        / "lee-wan-2008-mgo-pressure-correction.md"
    ).read_text(encoding="utf-8")
    candidates = json.loads(
        (ROOT / "docs" / "data" / "litcurate-eos-candidates.json").read_text(
            encoding="utf-8"
        )
    )["records"]
    same_doi = [
        row for row in candidates if row["publication"].get("doi", "").lower() == DOI
    ]
    assert len(same_doi) == 2
    assert all(audit.count(row["identifier"]) == 1 for row in same_doi)
