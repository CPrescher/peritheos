"""Source normalization, partial observations and exact thermal-law regression."""

import copy
import csv
import hashlib
import json
import subprocess
import sys
from collections import Counter

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    eosmat_schema,
    get_eos_record,
    get_material,
    get_material_document,
)
from scripts.reproduce_crichton_2016_vanadium import (
    CURVES,
    DATA,
    POINTS,
    RECORD,
    ROOT,
    reproduce,
    source_pressure,
)


def test_structure_reference_volume_and_uncertainty_provenance():
    document = get_material_document("vanadium_bcc")
    Draft202012Validator(eosmat_schema()).validate(document)
    assert document["formula"] == "V"
    assert document["phase"] == "bcc"
    assert document["space_group_number"] == 229
    assert document["formula_units_per_cell"] == 2
    site = document["atom_sites"][0]
    assert (site["x"], site["y"], site["z"], site["occupancy"]) == (0, 0, 0, 1)
    assert site["multiplicity"] * site["occupancy"] == 2
    assert document["lattice"]["a"] ** 3 == pytest.approx(27.655941287521)
    structure = document["source"]["structure"]
    assert (
        hashlib.sha256(
            (DATA / "vanadium-james-1960-9012770.cif").read_bytes()
        ).hexdigest()
        == structure["sha256"]
    )
    record = document["eos_records"][0]
    assert record["eos"]["parameters"] == {"V0": 27.995, "K0": 150.4, "K0_prime": 5.5}
    assert record["parameter_errors"] == {"V0": None, "K0": 6.2, "K0_prime": 1.0}
    assert record["parameter_error_confidence"] is None
    normalization = record["source_volume_normalization"]
    assert normalization["plot_denominator_a3"] * 1.005 == pytest.approx(27.995)
    assert normalization["reported_uncertainty"] == 0.001
    thermal = record["thermal"]
    assert thermal["thermal_expansion_law"] == "linear_temperature"
    assert thermal["reference_volume_law"] == "integrated_expansivity"
    assert thermal["parameter_errors"] == {
        "alpha0": 6e-6,
        "alpha1": 9e-9,
        "dK_dT": 0.0007,
    }
    assert thermal["fixed_parameters"] == ["Tr"]
    assert "fit_datasets" not in record  # Incomplete plot is not the original table.
    restored = Material.from_eosmat(Material.from_eosmat(document).to_eosmat())
    assert restored.to_eosmat()["datasets"] == document["datasets"]
    assert restored.eos_records[0].pressure(27.0, 1000) == pytest.approx(
        source_pressure(27.0, 1000)
    )


def test_digitized_observations_are_separate_from_curve_and_cited_data():
    document = get_material_document("vanadium_bcc")
    material = get_material("vanadium_bcc")
    for dataset in document["datasets"]:
        resource = dataset["resource"]
        payload = (ROOT / "peritheos/data" / resource["path"]).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == resource["sha256"]
        loaded = material.get_dataset(dataset["identifier"])
        assert loaded.checksum_verified
    with POINTS.open() as stream:
        rows = list(csv.DictReader(stream))
    assert Counter(r["series"] for r in rows) == {
        "filled": 29,
        "cited_ming_1978": 10,
        "recovered_anchor": 1,
    }
    filled = [r for r in rows if r["series"] == "filled"]
    assert Counter(r["temperature_nominal_k"] for r in filled) == {
        "300": 6,
        "450": 4,
        "700": 5,
        "900": 6,
        "1000": 8,
    }
    assert all(float(r["temperature_selection_half_width_k"]) == 20 for r in filled)
    assert float(filled[0]["pressure_gpa"]) == pytest.approx(80 * 12 / 690)
    assert float(filled[0]["relative_volume"]) == pytest.approx(0.94 + 489 * 0.08 / 693)
    with CURVES.open() as stream:
        curves = list(csv.DictReader(stream))
    assert len(curves) == 24
    assert max(float(r["pressure_gpa"]) for r in curves) == pytest.approx(11.35652174)
    # Test the retained pixel calibration, not model-generated expectations.
    for row in rows + curves:
        assert float(row["relative_volume"]) == pytest.approx(
            0.94 + (894 - float(row["pixel_y"])) * 0.08 / 693, abs=5e-9
        )


def test_independent_curve_reproduction_and_proxy_does_not_claim_parity():
    report = reproduce()
    assert report["curve_max_abs_relative_volume_difference"] < 0.001
    assert report["native_independent_pressure_max_difference_gpa"] < 2e-12
    assert report["alpha_300k_per_k"] == pytest.approx(4.08e-5)
    assert report["alpha_k_300k_gpa_per_k"] == pytest.approx(0.00616, abs=3e-5)
    assert report["k0_1000k_gpa"] == pytest.approx(119.18)
    assert report["expansion_300_to_1000k_percent"]["0"] == pytest.approx(2.2939147)
    assert abs(report["expansion_300_to_1000k_percent"]["0"] - 2.7) > 0.4
    proxy = report["nominal_temperature_proxy_fit"]
    assert proxy["status"] == "diagnostic_only"
    assert not proxy["original_fit_reproduced"]
    assert proxy["observations"] == 29
    assert proxy["published_rmse_gpa"] == pytest.approx(0.24263848)
    assert proxy["refit_rmse_gpa"] == pytest.approx(0.18247682)
    assert proxy["jacobian_rank"] == 6
    subprocess.run(
        [sys.executable, "scripts/reproduce_crichton_2016_vanadium.py", "--check"],
        cwd=ROOT,
        check=True,
    )


def test_native_fallback_derivatives_roundtrip_and_validity():
    record = get_eos_record(RECORD)
    model = record.eos
    fallback = copy.copy(model)
    del fallback._native
    volume = np.array([26.5, 27.0, 27.995])[:, None]
    temperature = np.array([300.0, 700.0, 1000.0])[None, :]
    assert model.pressure(volume, temperature) == pytest.approx(
        source_pressure(volume, temperature), abs=2e-12
    )
    assert fallback.pressure(volume, temperature) == pytest.approx(
        model.pressure(volume, temperature), abs=2e-12
    )
    h = 1e-5
    derivative = (
        -volume
        * (
            source_pressure(volume + h, temperature)
            - source_pressure(volume - h, temperature)
        )
        / (2 * h)
    )
    assert model.bulk_modulus(volume, temperature) == pytest.approx(
        derivative, rel=2e-8
    )
    assert model.pressure(27.995, 300) == pytest.approx(0.0, abs=1e-12)
    assert model.bulk_modulus(27.995, 300) == pytest.approx(150.4)
    assert model.thermal_pressure(volume, 300) == pytest.approx(
        np.zeros((3, 1)), abs=1e-12
    )
    for p, t in [(0, 300), (5, 700), (11.4, 1000)]:
        v = record.volume(p, t, check_validity=True)
        assert record.pressure(v, t, check_validity=True) == pytest.approx(p, abs=1e-10)
    for p, t in [(12, 300), (5, 1001), (5, 299), (60, 300)]:
        with pytest.raises(ValueError):
            record.volume(p, t, check_validity=True)


def test_investigation_and_refit_ledgers_keep_incomplete_alternatives():
    audit = json.loads((ROOT / "peritheos/data/primary-source-audit.json").read_text())
    entry = next(r for r in audit["records"] if r["record"] == RECORD)
    assert len(entry["withheld_candidates"]) == 2
    assert {r["parameters"]["K0_prime"] for r in entry["withheld_candidates"]} == {
        3.5,
        4,
    }
    ledger = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    row = next(r for r in ledger["records"] if r["record_identifier"] == RECORD)
    assert row["status"] == "not_refittable"
    assert row["partial_validation"]["original_fit_reproduced"] is False
    assert row["curve_validation"]["observations"] == 24
    paper_ledger = (ROOT / "docs/paper-investigation-ledger.md").read_text()
    assert "BM3 K0_prime fixed at 3.5" in paper_ledger
    assert "BM3 K0_prime fixed at 4" in paper_ledger
