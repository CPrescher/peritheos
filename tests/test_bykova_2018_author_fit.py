"""Scientific provenance and independent source-curve checks for the author files."""

import json
from pathlib import Path

import pytest

from scripts.reproduce_bykova_2018_author_fit import replay

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def source():
    return json.loads(
        (ROOT / "docs/data/bykova-2018-author-fit-inputs.json").read_text()
    )


@pytest.fixture(scope="module")
def result(source):
    return replay(source)


def test_recovered_selection_preserves_derived_ambient_and_curve_separately(source):
    assert len(source["points"]) == 22
    assert len(source["calculated_curve"]) == 100
    anchor = source["points"][14]
    assert anchor["origin"] == "derived_literature_ambient_eos_output"
    assert anchor["line"] == 23
    assert [float(t) for t in anchor["tokens"]] == [0.0001, 34.204, 0.0001, 0.041248]
    assert source["workbook_cells"]["LIT!Q2"]["value"] == 547.26
    assert source["workbook_cells"]["LIT!S2"]["formula"] == "=Q2/P2"
    assert source["workbook_cells"]["LIT!X2"]["number_format"] == "0.0"
    assert source["workbook_cells"]["LIT!V7"]["value"] == pytest.approx(35.9246)
    assert all(
        p["workbook_cells_p_sigma_p_v_sigma_v"][0] != "LIT!V7" for p in source["points"]
    )
    assert all(
        p["workbook_cells_p_sigma_p_v_sigma_v"][0] not in ("LIT!AA10", "LIT!AA11")
        for p in source["points"]
    )


def test_author_reported_pressure_weighted_replay_matches_printed_coefficients(result):
    fit = result["author_protocol_pressure_weighted_replay"]
    assert fit["parameters_v0_per_formula_unit_k0_gpa_kprime"] == pytest.approx(
        [34.2040, 102.73, 3.024], abs=0.005
    )
    assert fit["observations"] == 22
    assert fit["measured_rows"] == 21
    assert fit["volume_errors_used_in_objective"] is False
    assert result["parameter_agreement"]["central_coefficient_parity_achieved"] is True
    assert (
        result["parameter_agreement"][
            "author_fit_values_reproduced_at_reported_precision"
        ]
        is True
    )
    assert (
        result["parameter_agreement"]["each_parameter_within_printing_interval"]
        == [True] * 3
    )


def test_independent_saved_curve_checks_resolve_export_precision(result):
    curve = result["saved_calculated_curve"]
    assert curve["author_communicated_curve_compatible_rows"] == 100
    assert curve["workbook_precision_reconstruction_compatible_rows"] == 100
    # Shortened DAT observations reproduce published precision; the workbook
    # precision reconstruction also agrees with every original curve row at
    # that row's own printing resolution.
    assert curve["weighted_replay_compatible_rows"] == 90
    fit = result["workbook_precision_reconstruction"]
    assert fit["parameters_v0_per_formula_unit_k0_gpa_kprime"] == pytest.approx(
        [34.2038, 102.75, 3.023], abs=0.005
    )


def test_derived_anchor_is_a_material_constraint_not_an_independent_observation(result):
    fit = result["sensitivity_controls_not_replacement_eos"][
        "same_pressure_weights_without_derived_anchor_21"
    ]
    assert fit["parameters_v0_per_formula_unit_k0_gpa_kprime"][2] == pytest.approx(
        2.622779, abs=1e-5
    )
    assert fit["observations"] == fit["measured_rows"] == 21
    assert result["author_workflow"]["status"] == "not_assessed"


def test_generic_source_reaudit_preserves_author_file_evidence():
    from peritheos import get_material_document
    from scripts.apply_primary_source_audit import audit_record
    from scripts.reproduce_bykova_2018_author_fit import RECORD_ID

    record = next(
        r
        for r in get_material_document("coesite")["eos_records"]
        if r["identifier"] == RECORD_ID
    )
    previous = record["scientific_validation"]
    restored = audit_record(record, "coesite.eosmat")["scientific_validation"]
    for field in ("note", "audit_date", "primary_source_check", "primary_data_check"):
        assert restored[field] == previous[field]


def test_catalog_keeps_published_coefficients_and_original_author_input_types(source):
    import hashlib

    import numpy as np

    from peritheos import (
        convert_eosmat,
        get_eos_record,
        get_material,
        get_material_document,
    )
    from scripts.reproduce_bykova_2018_author_fit import DATASET_ID, DIGESTS, RECORD_ID

    document = convert_eosmat(get_material_document("coesite"), version=4, inline=True)
    record = next(r for r in document["eos_records"] if r["identifier"] == RECORD_ID)
    assert record["parameters"] == {
        "V0": {"value": 547.2, "error": 0.16},
        "K0": {"value": 103.0, "error": 2.0},
        "K0_prime": {"value": 3.02, "error": 0.15},
    }
    assert record["record_kind"] == "published"
    assessment = record["evidence"]["reproduction"]["assessment"]
    assert assessment["status"] == "parity"
    assert assessment["scope"] == "parameter_agreement"
    assert assessment["basis_code"] == "published_decimal_precision"
    assert "Central-coefficient parity is achieved" in assessment["qualifications"][0]
    assert "author's fit values are reproduced" in assessment["qualifications"][0]
    assert all(
        p["refit_error"] is None and p["within_combined_2sigma"] is None
        for p in assessment["parameters"]
    )
    assert (
        record["evidence"]["reproduction"]["historical_assessments"][0]["assessment"][
            "status"
        ]
        == "not_reproduced"
    )
    restored = convert_eosmat(
        convert_eosmat(document, version=3), version=4, inline=True
    )
    restored_record = next(
        r for r in restored["eos_records"] if r["identifier"] == RECORD_ID
    )
    assert (
        restored_record["evidence"]["reproduction"]
        == record["evidence"]["reproduction"]
    )
    data = get_material("coesite").get_dataset(DATASET_ID)
    assert len(data) == 22 and data.checksum_verified
    assert data["observation_origin"][14] == "derived_literature_ambient_eos_output"
    for i, point in enumerate(source["points"]):
        p, v, sp, sv = map(float, point["tokens"])
        np.testing.assert_allclose(
            [
                data["pressure_gpa"][i],
                data["volume_a3_z16_equivalent_cell"][i],
                data["pressure_sigma_gpa"][i],
                data["volume_sigma_a3_z16_equivalent_cell"][i],
            ],
            [p, v * 16, sp, sv * 16],
            rtol=0,
            atol=1e-12,
        )
    raw = (
        ROOT
        / "peritheos/data/datasets/coesite-i-iii-bykova-2018-author-eosfit-plot.dat"
    )
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == DIGESTS["coeII-plot.dat"]
    eos = get_eos_record(RECORD_ID)
    pressures = np.array([1.0, 15.0, 35.0])
    np.testing.assert_allclose(
        eos.pressure(eos.volume(pressures)), pressures, rtol=1e-10, atol=1e-9
    )


def test_current_ledger_and_generated_curves_replay_the_same_22_inputs():
    from peritheos import convert_eosmat, get_material_document
    from scripts.reconstruction_plot_data import curve_comparison
    from scripts.reproduce_bykova_2018_author_fit import RECORD_ID, ledger_outcome
    from scripts.standardize_reconstruction_metadata import assessment_for_record
    from scripts.validate_primary_eos_refits import check_ledger_values

    document = convert_eosmat(get_material_document("coesite"), version=4, inline=True)
    record = next(r for r in document["eos_records"] if r["identifier"] == RECORD_ID)
    rows = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())[
        "records"
    ]
    row = next(r for r in rows if r["record_identifier"] == RECORD_ID)
    outcome = ledger_outcome(record)
    assert outcome["central_coefficient_parity_achieved"] is True
    assert outcome["author_fit_values_reproduced_at_reported_precision"] is True
    check_ledger_values(outcome, {k: row[k] for k in outcome})
    assert (
        assessment_for_record(record, row)
        == record["evidence"]["reproduction"]["assessment"]
    )
    comparison, missing = curve_comparison(document, record)
    assert not missing and comparison is not None
    assert len(comparison.series.pressure) == 22
    assert comparison.source_rmse == pytest.approx(row["published_rmse_gpa"])
    assert comparison.refit_rmse == pytest.approx(row["rmse_gpa"])
    assert "derived ambient" in comparison.note
