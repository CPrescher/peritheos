import copy
import json
from pathlib import Path

import pytest

from peritheos import get_material_document
from scripts.reproduce_tange_2009_mgo import build_report

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs" / "data" / "tange-2009-mgo-partial-validation.json"
APPROXIMATE_REPORT = ROOT / "docs" / "data" / "tange-2009-mgo-approximate-refit.json"
RECORD_ID = "mgo_b1_tange_2009_vinet"
DATASET_IDS = [
    "mgo_dubrovinsky_1997_cod_thermal_expansion",
    "mgo_fiquet_1999_cod_thermal_expansion",
    "mgo_li_2006_table1_elasticity",
    "mgo_marsh_1980_lasl_single_crystal_hugoniot",
]


@pytest.fixture(scope="module")
def report():
    return build_report()


def test_tange_partial_validation_report_is_current(report, assert_audit_close):
    saved = json.loads(REPORT.read_text(encoding="utf-8"))
    replay = copy.deepcopy(report)
    # This incomplete global-fit diagnostic has weakly constrained parameters.
    # Keep the existing 0.1% relative / 0.0005 absolute comparison for its
    # diagnostics, while checking convergence and objective quality separately.
    partial = replay["bundled_partial_validation"].pop("partial_refit")
    expected_partial = saved["bundled_partial_validation"].pop("partial_refit")
    # Different platforms can meet different successful stopping criteria.
    # Reject failures/unknown messages instead of pinning one SciPy message.
    for fit in (partial, expected_partial):
        assert fit["solver_success"] is True
        assert fit.pop("solver_message") in {
            "`gtol` termination condition is satisfied.",
            "`ftol` termination condition is satisfied.",
            "`xtol` termination condition is satisfied.",
            "Both `ftol` and `xtol` termination conditions are satisfied.",
        }
    # Windows shifts the weakly constrained b by 5.08e-4. Bound this one
    # coefficient to 0.001 (~0.2%); do not relax other coefficients or metrics.
    assert partial["coefficients"].pop("b") == pytest.approx(
        expected_partial["coefficients"].pop("b"), rel=0, abs=1e-3
    )
    # Coefficient drift must still preserve the quality of the fitted objective
    # to 0.001%, substantially tighter than the individual group diagnostics.
    assert partial.pop("weighted_sum_of_squares") == pytest.approx(
        expected_partial.pop("weighted_sum_of_squares"), rel=1e-5, abs=0
    )
    assert_audit_close(partial, expected_partial, rel=1e-3, abs=5e-4)
    assert_audit_close(replay, saved)
    assert report["scope"] == "partial_validation_not_global_refit"
    assert report["global_refit_reproduced"] is False
    assert "Zha" in report["global_refit_blocker"]


def test_tange_source_inventory_keeps_missing_rows_explicit(report):
    inventory = {row["source"]: row for row in report["source_inventory"]}

    assert inventory["Dubrovinsky and Saxena (1997)"]["rows_recovered"] == 25
    assert inventory["Fiquet et al. (1999)"]["rows_recovered"] == 36
    assert inventory["Fiquet et al. (1999)"]["source_rows_stated_by_tange"] == 37
    assert inventory["Zha et al. (2000)"]["rows_recovered"] == 0
    assert report["bundled_partial_validation"]["observations"] == 102


def test_tange_local_reconstruction_is_qualified_similarity():
    approximate = json.loads(APPROXIMATE_REPORT.read_text(encoding="utf-8"))

    assert approximate["global_refit_reproduced"] is False
    assert approximate["result"] == "similar"
    assert approximate["observations"] == {
        "source_rows_stated_by_tange": 165,
        "rows_reconstructed": 164,
        "unresolved_rows": 1,
        "unresolved_source": "Fiquet et al. (1999)",
    }
    comparisons = approximate["approximate_refit"]["coefficient_comparisons"]
    assert all(item["similar"] for item in comparisons)
    assert all(
        item["all_parameters_similar"] for item in approximate["zha_weight_sensitivity"]
    )


def test_tange_published_coefficients_recover_source_group_residuals(report):
    metrics = report["bundled_partial_validation"]["published_coefficient_metrics"]

    assert metrics["lasl_hugoniot_pressure_gpa"]["rmse"] == pytest.approx(1.4, abs=0.05)
    assert metrics["li_2006_adiabatic_bulk_modulus_gpa"]["rmse"] == pytest.approx(
        1.5, abs=0.02
    )
    assert metrics["thermal_expansion_pressure_gpa"]["rmse"] < 0.3


def test_tange_partial_refit_is_not_misreported_as_global_parity(report):
    fit = report["bundled_partial_validation"]["partial_refit"]

    assert fit["solver_success"] is True
    assert fit["coefficients"]["a"] == pytest.approx(1.0, abs=1e-4)
    assert fit["coefficients"]["a"] != pytest.approx(0.138, abs=0.05)


def test_tange_record_links_only_the_redistributable_subset():
    document = get_material_document("mgo")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == RECORD_ID
    )
    datasets = {row["identifier"]: row for row in document["datasets"]}

    assert record["fit_datasets"] == DATASET_IDS
    assert record["scientific_validation"]["primary_data_check"]["status"] == "bundled"
    assert set(DATASET_IDS).issubset(datasets)
    assert (
        "partial validation"
        in record["scientific_validation"]["primary_data_check"]["finding"]
    )
