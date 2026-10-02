"""Prevent the Chen supporting study from becoming a false validated EOS."""

import csv
import hashlib
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import get_material, get_material_document
from scripts.reproduce_chen_2010_argon import (
    MASS_PER_CELL,
    POWERS,
    ROOT,
    density_pressure_derivative,
    heat_capacity_ratio,
    local_bm3_coefficients,
    read_points,
    reproduce,
)


def test_supporting_study_datasets_load_without_an_invented_eos():
    doc = get_material_document("argon_fcc")
    schema = json.loads((ROOT / "peritheos/data/eosmat-v3.schema.json").read_text())
    assert not list(Draft202012Validator(schema).iter_errors(doc))
    assert all(
        r.get("record_kind") in {"refit", "diagnostic"}
        for r in doc["eos_records"]
        if "chen_2010" in r["identifier"]
    )
    material = get_material("argon_fcc")
    datasets = [d for d in doc["datasets"] if "chen_2010" in d["identifier"]]
    assert len(datasets) == 3
    for dataset in datasets:
        assert dataset["used_by_eos_records"] == (
            ["argon_fcc_chen_2010_bm3_digitized_refit"]
            if dataset["identifier"].endswith("figure5_brillouin")
            else []
        )
        path = ROOT / "peritheos/data" / dataset["resource"]["path"]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )
        with path.open() as stream:
            assert next(csv.reader(stream)) == [c["name"] for c in dataset["columns"]]
        assert material.get_dataset(dataset["identifier"]).checksum_verified
    obs = material.get_dataset("argon_fcc_chen_2010_figure5_brillouin")
    assert len(obs) == 80
    np.testing.assert_allclose(
        obs["volume_a3"] * obs["density_g_cm3"], MASS_PER_CELL, rtol=3e-9
    )
    assert np.all(obs["pressure_plot_halfwidth_gpa"] > 0)
    assert np.all(obs["density_plot_halfwidth_g_cm3"] > 0)
    assert len(read_points("figure5-curve")) == 262
    assert len(read_points("reported-values")) == 14


def test_published_nonzero_reference_is_not_relabelled_zero_pressure():
    coefficients = local_bm3_coefficients()
    np.testing.assert_allclose(coefficients, [32.5325, -47.415, 16.8825], atol=1e-12)
    assert sum(coefficients) == pytest.approx(2)
    assert sum(POWERS * coefficients) == pytest.approx(15.1)
    assert sum(POWERS**2 * coefficients) / 15.1 == pytest.approx(5.4)
    result = reproduce()
    assert result["reproduction_status"] == "not_reproduced"
    audit = result["standard_bm3_local_constraint_diagnostic"]
    assert audit["native_pressure_max_difference_gpa"] < 2e-12
    assert audit["predicted_rho0_g_cm3"] == pytest.approx(1.6746473822)
    assert audit["density_marker_pressure_rms_gpa"] > 6
    shifted = result["shifted_bm3_at_2gpa_alternative_not_published_fit"]
    assert shifted["predicted_K_Ksecond_at_2gpa"] == pytest.approx(-7.3, abs=0.06)
    assert shifted["predicted_rho0_g_cm3"] == pytest.approx(1.6337365859)
    assert shifted["density_marker_pressure_rms_gpa"] > 6
    curve = result["figure5_curve_diagnostic"]
    assert curve["rho_at_2_gpa_g_cm3"] == pytest.approx(2.18, abs=0.01)
    assert 10.8 < curve["local_cubic_1_to_3_gpa_bulk_modulus_at_2_gpa"] < 11.5


def test_printed_acoustic_equation_units_and_invalid_denominator():
    assert heat_capacity_ratio(0) == 1.25
    assert heat_capacity_ratio(2) == pytest.approx(1.168588514060111)
    # Eq. 4 in km/s, g/cm^3 and GPa needs no extra factor of 1000.
    assert density_pressure_derivative(2, 3.228, 1.153) == pytest.approx(
        0.1351369531610182
    )
    with pytest.raises(ValueError, match="Positive acoustic bulk modulus"):
        density_pressure_derivative(2, 1, 2)


def test_packaged_source_manifest_preserves_status_and_pdf_identity():
    study = json.loads(
        (ROOT / "peritheos/data/studies/argon-fcc-chen-2010.json").read_text()
    )
    assert study["eos_record_identifiers"] == []
    assert study["reproduction_status"] == "not_reproduced"
    assert study["reference_temperature_k"] == 290
    assert study["reference_pressure_gpa"] == 2
    assert (
        study["pdf"]["sha256"]
        == "17a0bcd30604bd81f125894a18a45a91255b235fb80643cddd4bcbc6e5bd01f1"
    )
