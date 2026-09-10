import csv
import hashlib
import json
from pathlib import Path

import pytest

from peritheos import get_material_document
from scripts.reproduce_dorogokupets_2015_mgsio3 import (
    GUIGNOT,
    KOMABAYASHI,
    RECONSTRUCTION,
    REPORT,
    WANG,
    ZHOU,
    build_observations,
    reproduce,
    write_reconstruction,
)
from scripts.validate_primary_eos_refits import validate_all

ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = (
    ROOT / "peritheos/data/datasets/dorogokupets-2015-298k-reconstruction-source.json"
)

TARGETS = {
    "akimotoite": (
        "akimotoite_dorogokupets_2015_298k_rydberg_stacey",
        "dorogokupets_2015_akimotoite_298k_reconstruction",
    ),
    "bridgmanite": (
        "bridgmanite_dorogokupets_2015_298k_rydberg_stacey",
        "dorogokupets_2015_bridgmanite_298k_reconstruction",
    ),
    "mgsio3_post_perovskite": (
        "mgsio3_post_perovskite_dorogokupets_2015_298k_rydberg_stacey",
        "dorogokupets_2015_post_perovskite_298k_reconstruction",
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_reconstruction_writer_preserves_git_normalized_bytes(tmp_path):
    output = tmp_path / "reconstruction.csv"
    write_reconstruction(output)
    assert output.read_bytes() == RECONSTRUCTION.read_bytes()
    assert b"\r" not in output.read_bytes()


def test_reconstructed_source_rows_and_pressure_coordinates():
    observations = build_observations()
    grouped = {
        phase: [item for item in observations if item.phase == phase]
        for phase in {
            "bridgmanite",
            "akimotoite_ruby",
            "akimotoite_ice_vii",
            "akimotoite_wang",
            "post_perovskite",
        }
    }

    assert {phase: len(rows) for phase, rows in grouped.items()} == {
        "bridgmanite": 19,
        "akimotoite_ruby": 16,
        "akimotoite_ice_vii": 14,
        "akimotoite_wang": 23,
        "post_perovskite": 27,
    }
    assert len({item.source_order for item in grouped["akimotoite_ruby"]}) == 16
    assert grouped["post_perovskite"][0].pressure_gpa == pytest.approx(112.0738968631)
    assert all(
        "Sokolova-2013 MgO" in item.pressure_coordinate
        for item in grouped["post_perovskite"]
    )

    with RECONSTRUCTION.open(newline="", encoding="utf-8") as handle:
        generated_rows = list(csv.DictReader(handle))
    assert len(generated_rows) == len(observations) == 99
    with GUIGNOT.open(newline="", encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle))) == 24
    with WANG.open(newline="", encoding="utf-8") as handle:
        wang = list(csv.DictReader(handle))
    assert len(wang) == 134
    assert {row["run"] for row in wang} == {"T0133", "T0150"}
    assert sum(row["run"] == "T0133" for row in wang) == 78
    with KOMABAYASHI.open(newline="", encoding="utf-8") as handle:
        komabayashi = list(csv.DictReader(handle))
    assert len(komabayashi) == 22
    assert sum(row["ppv_fit_included"] == "true" for row in komabayashi) == 21
    with ZHOU.open(newline="", encoding="utf-8") as handle:
        zhou = list(csv.DictReader(handle))
    assert len(zhou) == 58
    assert sum(bool(row["ks_gpa"]) for row in zhou) == 55


def test_diagnostic_fits_are_stable_and_do_not_claim_published_parity(
    assert_audit_close,
):
    fits = reproduce()["fits"]

    assert fits["bridgmanite"]["unweighted"]["K0_gpa"] == pytest.approx(256.4434912813)
    assert fits["bridgmanite"]["unweighted"]["K0_prime"] == pytest.approx(3.8501759713)
    assert fits["akimotoite_ruby"]["unweighted"]["K0_prime"] == pytest.approx(
        6.5448216288
    )
    assert fits["akimotoite_ice_vii"]["unweighted"]["K0_prime"] == pytest.approx(
        3.1801339544
    )
    assert fits["akimotoite_wang"]["unweighted"]["K0_prime"] == pytest.approx(1.0)
    assert fits["akimotoite_wang_plus_reynard_ruby"]["unweighted"][
        "K0_gpa"
    ] == pytest.approx(254.9276897403)
    assert fits["post_perovskite"]["unweighted"]["K0_gpa"] == pytest.approx(
        246.7833658099
    )
    for result in fits.values():
        assert result["unweighted"]["success"]
        assert result["reported_uncertainty_diagnostic"]["success"]
        assert (
            "not the unpublished Dorogokupets objective weights"
            in result["reported_uncertainty_diagnostic"]["weighting"]
        )

    elasticity = reproduce()["zhou_2014_elasticity"]
    assert elasticity["observations"] == 55
    assert elasticity["unweighted"]["K0S_gpa"] == pytest.approx(219.8672626016)
    assert elasticity["unweighted"]["dKS_dP"] == pytest.approx(4.6070287481)
    assert elasticity["reported_KS_uncertainty_sensitivity"]["dKS_dP"] == pytest.approx(
        4.6197293203
    )

    checked_report = json.loads(REPORT.read_text(encoding="utf-8"))
    assert_audit_close(fits, checked_report["fits"])
    assert checked_report["reconstruction_sha256"] == _sha256(RECONSTRUCTION)


def test_material_provenance_and_common_ledger_expose_partial_reconstruction():
    reconstruction_hash = _sha256(RECONSTRUCTION)
    provenance_hash = _sha256(PROVENANCE)
    for material_id, (record_id, dataset_id) in TARGETS.items():
        document = get_material_document(material_id)
        record = next(
            item for item in document["eos_records"] if item["identifier"] == record_id
        )
        dataset = next(
            item for item in document["datasets"] if item["identifier"] == dataset_id
        )
        assert record["diagnostic_datasets"] == [dataset_id]
        assert record["scientific_validation"]["primary_data_check"]["status"] == (
            "partial_reconstruction"
        )
        assert dataset["resource"]["sha256"] == reconstruction_hash
        assert dataset["provenance_resource"]["sha256"] == provenance_hash
        assert dataset["uncertainty"]["covariance"] == "not_reported"

    post_perovskite = get_material_document("mgsio3_post_perovskite")
    guignot = next(
        item
        for item in post_perovskite["datasets"]
        if item["identifier"]
        == "mgsio3_post_perovskite_guignot_2007_table1_300k_compression"
    )
    assert guignot["resource"]["sha256"] == _sha256(GUIGNOT)

    akimotoite = get_material_document("akimotoite")
    for dataset_id, path in (
        ("akimotoite_wang_2004_tables1_2_pvt", WANG),
        ("akimotoite_zhou_2014_table1_elasticity", ZHOU),
    ):
        dataset = next(
            item for item in akimotoite["datasets"] if item["identifier"] == dataset_id
        )
        assert dataset["resource"]["sha256"] == _sha256(path)

    komabayashi = next(
        item
        for item in post_perovskite["datasets"]
        if item["identifier"] == "mgsio3_post_perovskite_komabayashi_2008_table1_pvt"
    )
    assert komabayashi["resource"]["sha256"] == _sha256(KOMABAYASHI)

    ledger_records = {
        item["record_identifier"]: item for item in validate_all()["records"]
    }
    for record_id, _dataset_id in TARGETS.values():
        outcome = ledger_records[record_id]
        assert outcome["status"] == "parity_not_achieved"
        assert "partial" in outcome["reason"].lower()
        assert outcome["solver_success"]
