import csv
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import get_material_document, recalculate_ruby_pressure
from scripts.reproduce_datchi_2007_diamond import reproduce

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "diamond-occelli-2003-figure2-digitized.csv"
)
REPORT = ROOT / "docs" / "data" / "datchi-2007-diamond-refit.json"


def test_occelli_figure2_bundle_and_record_provenance():
    document = get_material_document("diamond")
    record = next(
        row
        for row in document["eos_records"]
        if row["identifier"] == "diamond_datchi_2007_vinet_1"
    )
    dataset = next(
        row
        for row in document["datasets"]
        if row["identifier"] == "diamond_occelli_2003_figure2_digitized"
    )
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 24
    assert sum(row["run"] == "2" for row in rows) == 14
    assert sum(row["run"] == "3" for row in rows) == 10
    assert record["fixed_parameters"] == ["V0"]
    assert record["eos"]["parameters"]["V0"] == pytest.approx(5.6733 * 8)
    assert record["fit_datasets"] == [dataset["identifier"]]
    assert record["pressure_calibration"]["status"] == "resolved"
    assert record["pressure_calibration"]["recalculation"]["status"] == "ready"
    assert record["scientific_validation"]["primary_data_check"]["status"] == (
        "plot_only"
    )
    assert dataset["provenance"]["quality_control"]["covariance_inferred"] is False


def test_exact_mxb1986_to_h2005_columns():
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    source = np.array([float(row["pressure_mxb1986_gpa"]) for row in rows])
    stored = np.array([float(row["pressure_h2005_gpa"]) for row in rows])
    recalculated = recalculate_ruby_pressure(
        source,
        source_calibration="ruby_mao_1986",
        target_calibration="ruby_holzapfel_2005",
    )

    assert recalculated == pytest.approx(stored, abs=5.0e-9)
    assert stored[-1] == pytest.approx(154.975338747)


def test_dedicated_table_ii_audit_is_current_and_covariance_free():
    checked_in = json.loads(REPORT.read_text(encoding="utf-8"))
    regenerated = reproduce()
    assert regenerated["dataset_sha256"] == checked_in["dataset_sha256"]
    assert regenerated["pressure_recalibration"] == checked_in["pressure_recalibration"]

    mxb = checked_in["fits"]["MXB1986"]["best_supported_run2_selection"][
        "equal_volume_weight_fit"
    ]
    h05 = checked_in["fits"]["H2005"]["best_supported_run2_selection"][
        "equal_volume_weight_fit"
    ]
    assert [mxb["parameters"]["K0"], mxb["parameters"]["K0_prime"]] == (
        pytest.approx([446.7119197742, 2.9901052362])
    )
    assert [h05["parameters"]["K0"], h05["parameters"]["K0_prime"]] == (
        pytest.approx([444.3462890945, 3.9607004387])
    )
    assert h05["volume_chi_square_using_reported_constant_sigma"] == pytest.approx(
        1.8176627891
    )
    assert all(value is None for value in h05["parameter_errors"].values())
    regenerated_h05 = regenerated["fits"]["H2005"]["best_supported_run2_selection"][
        "equal_volume_weight_fit"
    ]
    assert regenerated_h05["parameters"] == pytest.approx(h05["parameters"])
    assert regenerated_h05["volume_chi_square_using_reported_constant_sigma"] == (
        pytest.approx(h05["volume_chi_square_using_reported_constant_sigma"])
    )
    assert checked_in["fit_protocol"]["published_exclusions"] is None
    assert checked_in["fit_protocol"]["published_weights"] is None
    assert "not estimated" in checked_in["fit_protocol"]["covariance"]
