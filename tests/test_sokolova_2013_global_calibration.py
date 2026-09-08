import json
from pathlib import Path

import pytest

from peritheos import get_material_document
from scripts.reproduce_sokolova_2013_global_calibration import reconstruct

ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "docs" / "data" / "sokolova-2013-global-calibration.json"
REFIT_LEDGER_PATH = ROOT / "docs" / "data" / "primary-eos-refits.json"
MANIFEST_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "sokolova-2013-global-calibration-manifest.json"
)


def test_sokolova_global_reconstruction_uses_all_markers_in_one_objective():
    committed = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    replay = reconstruct()

    assert replay == committed
    assert committed["classification"] == (
        "source_constrained_cross_calibration_reconstruction"
    )
    assert committed["independent_eos_refit"] is False
    assert committed["markers"] == 11
    assert committed["observations"] == 392
    assert committed["published_parameters"] == {"A_gpa": 1870.0, "m": 6.0}
    assert committed["sensitivity_fit"]["parameters"]["A_gpa"] == pytest.approx(
        1869.3837841
    )
    assert committed["sensitivity_fit"]["parameters"]["m"] == pytest.approx(6.10334763)
    assert committed["published_calibration_marker_equal_rmse_gpa"] == pytest.approx(
        1.2774517231
    )
    assert committed["fitted_calibration_marker_equal_rmse_gpa"] == pytest.approx(
        1.2639397186
    )
    assert committed["table4_closure_marker_equal_rmse_gpa"] == pytest.approx(
        1.1474132284
    )


def test_sokolova_manifest_preserves_sources_constraints_and_workbook_custody():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    assert manifest["references"] == {
        "method": "10.1103/PhysRevB.75.024115",
        "calibration": "10.1016/j.rgg.2013.01.005",
        "workbooks": "10.1016/j.cageo.2016.06.002",
    }
    assert len(manifest["markers"]) == 11
    assert sum(marker["n"] for marker in manifest["markers"]) == 12
    assert manifest["reconstructable_cross_calibration"]["constraints"] == {
        "A_gpa": [1500.0, 2200.0],
        "m": [0.0, 20.0],
    }
    assert len(manifest["official_2016_workbooks"]["files"]) == 11
    assert (
        len({item["sha256"] for item in manifest["official_2016_workbooks"]["files"]})
        == 11
    )
    assert "excluded" in manifest["official_2016_workbooks"]["graphs_sheet"]
    assert any(
        "residual weights" in item
        for item in manifest["published_objective"]["not_published"]
    )


def test_sokolova_records_classify_reconstruction_separately_from_refitting():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    refit_ledger = json.loads(REFIT_LEDGER_PATH.read_text(encoding="utf-8"))
    outcomes = {item["record_identifier"]: item for item in refit_ledger["records"]}
    for marker in manifest["markers"]:
        record = next(
            item
            for item in get_material_document(marker["material"])["eos_records"]
            if item["identifier"] == marker["record_identifier"]
        )
        check = record["scientific_validation"]["primary_data_check"]
        reconstruction = record["scientific_validation"]["source_reconstruction"]
        assert check["status"] == "parameterization_only"
        assert check["comparison_dataset_identifiers"] == [marker["dataset_identifier"]]
        assert reconstruction["status"] == "coupled_cross_calibration_reconstructed"
        assert reconstruction["observations"] == 392
        assert reconstruction["independent_eos_refit"] is False
        assert "not an independent refit" in record["notes"]
        outcome = outcomes[marker["record_identifier"]]
        assert outcome["status"] == "source_reconstruction"
        assert outcome["dataset_identifiers"] == [marker["dataset_identifier"]]
        assert outcome["independent_eos_refit"] is False
