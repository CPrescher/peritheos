"""Primary-source checks for Hirose et al. (2005) MgGeO3."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_hirose_2005_mggeo3 import reproduce

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/mggeo3-post-perovskite-hirose-2005-table3-pv.csv"
IDENTIFIERS = {
    "mggeo3_post_perovskite_hirose_2005_bm3_fixed_1",
    "mggeo3_post_perovskite_hirose_2005_bm3_free_2",
}


def test_dataset_hash_and_complete_table3():
    assert hashlib.sha256(DATA.read_bytes()).hexdigest() == (
        "280487791ea68da488f3c321f59ee84fe05f3bf4d512a5a5d6f5ebb9e9b1160d"
    )
    assert len(DATA.read_text(encoding="utf-8").splitlines()) == 10


def test_both_published_fits_are_executable():
    document = get_material_document("mggeo3_post_perovskite")
    loaded = Material.from_eosmat(document, record_identifiers=sorted(IDENTIFIERS))
    assert {record.identifier for record in loaded.eos_records} == IDENTIFIERS
    for record in loaded.eos_records:
        assert record.eos.pressure(record.eos.V0) == pytest.approx(0.0)
        assert record.eos.bulk_modulus(record.eos.V0) == pytest.approx(record.eos.K0)


def test_diagnostic_refits_recover_published_uncertainty_regions():
    metrics = reproduce()
    assert metrics["fixed"]["rows"] == 9
    assert metrics["fixed"]["refit_parameters"] == pytest.approx(
        (183.01609525, 192.6917299, 4.0), rel=2.0e-8
    )
    assert metrics["free"]["refit_parameters"] == pytest.approx(
        (182.17460987, 212.62137272, 3.46315314), rel=1.0e-6
    )
    assert metrics["fixed"]["published_rmse_gpa"] < 0.88
    assert metrics["free"]["published_rmse_gpa"] < 0.90


def test_source_preference_and_pressure_calibration_are_explicit():
    document = get_material_document("mggeo3_post_perovskite")
    records = {
        record["identifier"]: record
        for record in document["eos_records"]
        if record["identifier"] in IDENTIFIERS
    }
    assert records["mggeo3_post_perovskite_hirose_2005_bm3_fixed_1"][
        "fixed_parameters"
    ] == ["K0_prime"]
    assert (
        records["mggeo3_post_perovskite_hirose_2005_bm3_free_2"]["fixed_parameters"]
        == []
    )
    assert all(
        record["pressure_calibration"]["methods"][0]["reference"]["doi"]
        == "10.1063/1.344177"
        for record in records.values()
    )
