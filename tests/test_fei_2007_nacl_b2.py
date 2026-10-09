"""NaCl-specific checks of the primary Fei (2007) equation and figure recovery."""

import csv
import hashlib
import json

import numpy as np
import pytest

from peritheos import Material
from peritheos.eosmat import validate_eosmat_document
from peritheos.materials import NACL_B2_FEI_2007
from scripts.reproduce_fei_2007_nacl_b2 import (
    DATA,
    OUTPUT,
    ROOT,
    SOURCE,
    pressure,
    subset,
)


def document():
    return json.loads((ROOT / "peritheos/data/materials/nacl_b2.eosmat").read_text())


def test_canonical_record_preserves_published_equation_and_formula_unit_basis():
    raw = document()
    validate_eosmat_document(raw)
    record = Material.from_eosmat(
        raw, record_identifiers=["nacl_b2_fei_2007"]
    ).eos_records[0]
    v = np.array([20.5, 23.0, 27.0, 35.0, 41.35])[:, None]
    t = np.array([300.0, 1000.0, 2000.0])[None, :]
    expected = pressure(v, t)
    assert record.pressure(v, t) == pytest.approx(expected, abs=1e-9)
    assert record.pressure(v, t) == pytest.approx(
        NACL_B2_FEI_2007.pressure(v, t), abs=1e-9
    )
    assert record.eos.n == 2
    assert record.eos.configuration_values() == {
        "debye_temperature_law": "variable_exponent"
    }
    assert record.eos.thermal_pressure(v * record.volume_scale, 300) == pytest.approx(
        np.zeros((5, 1)), abs=1e-12
    )
    source = next(r for r in raw["eos_records"] if r["identifier"] == record.identifier)
    assert source["volume_basis"]["formula_units"] == 1
    assert source["parameter_errors"] == {"V0": None, "K0": 2.9, "K0_prime": 0.26}
    assert source["thermal"]["parameter_errors"]["q"] == 0.3
    assert source["parameter_error_confidence"] is None
    assert source["parameter_covariance"] is None
    assert source.get("default") is not True


def test_recovered_markers_include_touching_hot_pair_and_exclude_comparisons():
    source = json.loads(SOURCE.read_text())
    assert len(source["observations"]) == 36
    assert len(subset(source, 3, 300)[0]) == 12
    assert len(subset(source, 4, 300)[0]) == 12
    assert len(subset(source, 4, 1000)[0]) == 12
    hot = [r for r in source["observations"] if r["temperature_k"] == 1000]
    assert [r["pdf_curve_index"] for r in hot] == list(range(126, 138))
    assert hot[7]["x_pt"] != hot[8]["x_pt"]
    assert {
        r["pdf_curve_index"] for r in source["observations"] if r["figure"] == 3
    } == set(range(84, 96))
    assert all(
        c["pdf_curve_index"]
        not in {r["pdf_curve_index"] for r in source["observations"]}
        for c in source["calculated_curves"]
    )


def test_csv_provenance_hashes_and_missing_errors_are_explicit():
    for dataset in document()["datasets"]:
        if not dataset["identifier"].startswith("nacl_b2_fei_2007_"):
            continue
        path = DATA / dataset["resource"]["path"].split("/")[-1]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == (12 if "figure3" in path.name else 24)
        assert all(
            row[column] == ""
            for row in rows
            for column in (
                "pressure_error_gpa",
                "volume_error_a3",
                "temperature_error_k",
                "platinum_volume_a3",
            )
        )


def test_report_keeps_equation_curve_and_author_fit_claims_separate():
    report = json.loads(OUTPUT.read_text())
    assert report["source_sha256"] == hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert report["native_max_difference_gpa"] < 1e-9
    assert report["author_fit_reproduced"] is False
    assert report["source_error_weighted_fit"]["status"] == "not_possible"
    assert report["cold_fit_figure3"]["observations"] == 12
    assert report["hot_q_only_published_cold"]["observations"] == 12
    assert report["hot_excluding_touching_points_sensitivity"]["observations"] == 10
    assert report["hot_q_only_published_cold"]["parameters"]["q"] > 0.9
    assert report["hot_excluding_touching_points_sensitivity"]["parameters"]["q"] < 0.7
    assert report["figure3_vs_figure4_same_cold_rows"][
        "duplicate_observations_not_pooled"
    ]
    assert all(
        c["role"] == "calculated curve, excluded from fits"
        for c in report["curve_checks"]
    )
    assert all(
        c["max_horizontal_difference_pt"] < c["printed_linewidth_pt"]
        for c in report["curve_checks"]
    )


def test_refit_ledger_adapter_does_not_pool_duplicate_cold_markers():
    from scripts.reproduce_fei_2007_nacl_b2 import ledger_outcome

    record = next(
        r for r in document()["eos_records"] if r["identifier"] == "nacl_b2_fei_2007"
    )
    outcome = ledger_outcome(record)
    assert outcome["status"] == "parity_not_achieved"
    assert outcome["observations"] == 24
    assert outcome["solver_success"]
    assert outcome["parameters"][-1]["parameter"] == "q"
    assert outcome["parameters"][-1]["similar"] is False
    assert outcome["parameters"][-1]["within_combined_2sigma"] is None
