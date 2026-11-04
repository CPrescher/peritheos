"""Preserve the distinct author fit and published structural observations."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import get_material, get_material_document
from peritheos.eosmat import convert_eosmat, load_eosmat, validate_eosmat_document
from scripts.reproduce_prescher_2015_fe7c3 import (
    DATA,
    DATASET_ID,
    RECORD_ID,
    ledger_outcome,
    observations,
)

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "peritheos/data/datasets"


def test_author_input_preserves_every_raw_token_and_table_s4_is_unchanged():
    raw = (DATA_DIR / "iron-carbide-fe7c3-prescher-2015-author-input.txt").read_text()
    with DATA.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 28
    for row, line in zip(rows, raw.splitlines()[1:]):
        assert [
            row[name]
            for name in (
                "pressure_gpa",
                "pressure_uncertainty_gpa",
                "unit_cell_volume_a3",
                "unit_cell_volume_uncertainty_a3",
            )
        ] == line.split()
    assert rows[1]["run_id"] == "a_p1_rot"
    assert rows[14]["unit_cell_volume_a3"] == "571.8"
    assert rows[25]["pressure_gpa"] == "153.43"
    table = DATA_DIR / "iron-carbide-fe7c3-prescher-2015-table-s4-compression.csv"
    assert hashlib.sha256(table.read_bytes()).hexdigest() == (
        "744fc747827c776d347dcc7f221782a23dbd17534357cafd817952fab686617e"
    )


def test_fit_association_and_table_prohibition_survive_native_legacy_roundtrip():
    path = ROOT / "peritheos/data/materials/iron_carbide_fe7c3.eosmat"
    native = load_eosmat(path)
    validate_eosmat_document(native)
    roundtrip = convert_eosmat(convert_eosmat(native, version=3), version=4)
    assert (
        roundtrip["eos_records"][0]["parameters"]
        == native["eos_records"][0]["parameters"]
    )
    for document in (native, roundtrip):
        record = document["eos_records"][0]
        assert record["fit_datasets"] == [DATASET_ID]
        assert record["parameters"]["K0"]["value"] == 168.0
        table, author = document["datasets"]
        assert table["usage_recommendation"] == "do_not_use_for_scientific_calculations"
        assert "should not be used" in table["notes"]
        assert table["used_by_eos_records"] == []
        assert author["used_by_eos_records"] == [RECORD_ID]
        historical = record["evidence"]["reproduction"][
            "historical_table_s4_diagnostic"
        ]
        assert historical["assessment"]["fit"]["observations"] == 27
        assert historical["assessment"]["parameters"][1]["refit"] == pytest.approx(
            188.4191063
        )


def test_reconstruction_recovers_saved_coefficients_and_errors():
    record = get_material_document("iron_carbide_fe7c3")["eos_records"][0]
    result = ledger_outcome(record)
    assert result["status"] == "parity"
    assert result["dataset_identifiers"] == [DATASET_ID]
    assert result["observations"] == 28
    assert result["degrees_of_freedom"] == 25
    assert result["author_workflow_reproduced"] is False
    assert result["author_fit_export_reproduced_at_printed_precision"] is True
    assert [p["refit"] for p in result["parameters"]] == pytest.approx(
        [748.09, 168.15, 6.06], abs=0.005
    )
    assert [p["refit_error"] for p in result["parameters"]] == pytest.approx(
        [1.05, 4.11, 0.13], abs=0.005
    )
    assert result["published_rmse_gpa"] == pytest.approx(0.77625928, abs=1e-7)
    assert result["absolute_sigma"] is False
    assert np.linalg.eigvalsh(result["parameter_covariance"]).min() > 0


def test_installed_dataset_access_and_generic_reaudit_keep_the_two_sources_separate():
    from scripts.apply_primary_source_audit import audit_record

    material = get_material("iron_carbide_fe7c3")
    author = material.get_dataset(DATASET_ID)
    table = material.get_dataset(
        "iron_carbide_fe7c3_prescher_2015_table_s4_compression"
    )
    assert len(author) == 28 and len(table) == 27
    assert author.checksum_verified and table.checksum_verified
    np.testing.assert_array_equal(
        author.as_pressure_volume().volume,
        observations()["unit_cell_volume_a3"],
    )
    record = get_material_document("iron_carbide_fe7c3")["eos_records"][0]
    reaudited = audit_record(record, "iron_carbide_fe7c3.eosmat")
    assert reaudited["fit_datasets"] == [DATASET_ID]
    assert (
        reaudited["scientific_validation"]["primary_data_check"]
        == record["scientific_validation"]["primary_data_check"]
    )


def test_plot_replays_author_observations_and_keeps_printed_zero_widths():
    from scripts.reconstruction_plot_data import observations as plot_observations

    document = get_material_document("iron_carbide_fe7c3")
    record = document["eos_records"][0]
    fit = record["reproduction"]["assessment"]["fit"]
    plot = plot_observations(document, record, fit)
    rows = observations()
    np.testing.assert_array_equal(plot.volume, rows["unit_cell_volume_a3"])
    np.testing.assert_array_equal(plot.pressure_sigma, rows["pressure_uncertainty_gpa"])
    assert np.count_nonzero(plot.pressure_sigma == 0) == 3


def test_published_curve_is_unchanged_and_matches_independent_bm3_calculation():
    record = get_material("iron_carbide_fe7c3").get_eos_record(RECORD_ID)
    v = np.array([732.3, 571.8, 532.1])
    eta = (748.0 / v) ** (1 / 3)
    expected = 1.5 * 168 * (eta**7 - eta**5) * (1 + 0.75 * (6.1 - 4) * (eta**2 - 1))
    p = record.pressure(volume=v, temperature=300)
    np.testing.assert_allclose(p, expected, rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(
        record.volume(pressure=p, temperature=300), v, rtol=1e-10
    )


def test_packaged_provenance_preserves_log_discrepancy_and_raw_resource_checksums():
    source = json.loads(
        (DATA_DIR / "iron-carbide-fe7c3-prescher-2015-source.json").read_text()
    )
    for key in ("author_input", "author_fit_export", "normalized_csv"):
        resource = source[key]
        path = ROOT / "peritheos/data" / resource["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == resource["sha256"]
    assert len(source["observations"]) == 28
    assert len(source["lineage"]["logs"]) == 4
    p9 = source["lineage"]["logs"][0]
    assert "571.8(3)" in p9["orthorhombic_refined_cell"]["text"]
    assert "564.726811" in p9["reduction_cell"]["text"]
