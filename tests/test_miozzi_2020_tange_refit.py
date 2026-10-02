"""Selectable Tange-calibrated Fe refit with distinct source provenance."""

import csv
import hashlib
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    eosmat_schema,
    get_eos_record,
    get_eos_record_document,
    get_material,
    get_material_document,
)
from scripts.reconstruct_miozzi_2020_iron import fe_pressure, read_data
from scripts.refit_miozzi_2020_tange import (
    CALIBRATION_ID,
    CSV_PATH,
    DATASET_ID,
    NAMES,
    RECORD_ID,
    REPORT,
    TANGE,
    derived_rows,
    ledger_outcome,
    reproduce,
    tange_pressure,
    target_pressures,
)
from scripts.reproduce_miozzi_2020_iron import RECORDS


@pytest.mark.parametrize("model", ["vinet", "bm3"])
def test_primary_tange_table5_checkpoint(model):
    # Independent Table 5 values at x=.8, T=300/2000 K, rounded to 0.01 GPa.
    expected = [57.57, 68.91] if model == "vinet" else [57.55, 68.07]
    np.testing.assert_allclose(
        tange_pressure(0.8 * TANGE[model]["V0"], [300, 2000], model),
        expected,
        atol=0.01,
        rtol=0,
    )


def test_derived_pressures_preserve_sources_and_native_calibration():
    original, data, hashes = read_data()
    rows, text = derived_rows(original, data)
    assert CSV_PATH.read_text() == text
    assert len(rows) == 131
    assert {row["source_medium"] for row in rows} == {"he", "mgo"}
    target = target_pressures(data, "vinet")
    mask = data["mgo"]
    np.testing.assert_array_equal(target[:15], data["pressure_gpa"][:15])
    native = get_eos_record(CALIBRATION_ID).pressure(
        data["mgo_volume_a3"][mask], data["temperature_k"][mask]
    )
    np.testing.assert_allclose(target[mask], native, rtol=1e-11, atol=1e-9)
    for source, derived in zip(original, rows):
        assert derived["printed_pressure_gpa"] == source["pressure_gpa"]
        assert derived["printed_pressure_error_gpa"] == source["pressure_error_gpa"]
        for name in (
            "temperature_k",
            "temperature_error_k",
            "volume_a3",
            "volume_error_a3",
        ):
            assert derived[name] == source[name]
        if derived["source_medium"] == "he":
            assert derived["pressure_temperature_covariance_gpa_k"] == ""
    assert rows[15 + 26]["volume_error_a3"] == ""
    assert rows[15 + 75]["conditional_pressure_error_gpa"] == ""
    report = reproduce()
    assert report["original_csv_sha256"] == hashes
    assert report["derived_csv_sha256"] == hashlib.sha256(text.encode()).hexdigest()
    ds = next(
        d
        for d in get_material_document("iron")["datasets"]
        if d["identifier"] == DATASET_ID
    )
    assert ds["resource"]["sha256"] == report["derived_csv_sha256"]
    with CSV_PATH.open() as stream:
        assert list(csv.DictReader(stream)) == rows


def test_refit_is_selectable_with_errors_covariance_and_lossless_export():
    document = get_material_document("iron")
    Draft202012Validator(eosmat_schema()).validate(document)
    raw = get_eos_record_document(RECORD_ID)
    assert raw["record_kind"] == "refit"
    assert raw["default"] is False
    assert "Peritheos refit" in raw["label"]
    assert "Tange" in raw["label"]
    assert (
        raw["scientific_validation"]["original_publication_reproduction_status"]
        == "not_reproduced"
    )
    assert (
        get_material("iron").default_record().identifier
        == "iron_sakai_2025_rydberg_stacey_1"
    )
    # Strict loading is permitted for this numerically reproduced independent
    # refit; unreproduced author coefficients keep their own restriction.
    material = Material.from_eosmat(document, record_identifiers=[RECORD_ID])
    exported = material.to_eosmat()
    for key in (
        "parameter_errors",
        "parameter_covariance",
        "fit_provenance",
        "scientific_validation",
        "reproduction",
    ):
        assert exported["eos_records"][0][key] == raw[key]
    assert exported["datasets"][0]["identifier"] == DATASET_ID
    assert (
        exported["datasets"][0]["pressure_reconstruction"]
        == document["datasets"][-1]["pressure_reconstruction"]
    )
    covariance = np.array(raw["parameter_covariance"]["matrix"])
    errors = {**raw["parameter_errors"], **raw["thermal"]["parameter_errors"]}
    assert covariance.shape == (5, 5)
    assert np.all(np.linalg.eigvalsh(covariance) > 0)
    np.testing.assert_allclose(
        np.diag(covariance), [errors[name] ** 2 for name in NAMES]
    )
    assert raw["thermal"]["parameter_errors"]["theta0"] is None
    assert raw["thermal"]["fixed_parameters"] == ["Tr", "theta0", "n"]
    original, data, _ = read_data()
    parameters = [
        raw["eos"]["parameters"][name]
        if name in NAMES[:3]
        else raw["thermal"]["parameters"][name]
        for name in NAMES
    ]
    native = get_eos_record(RECORD_ID)
    independent = fe_pressure(data["volume_a3"], data["temperature_k"], parameters)
    np.testing.assert_allclose(
        native.pressure(data["volume_a3"], data["temperature_k"]),
        independent,
        rtol=2e-11,
        atol=1e-9,
    )
    np.testing.assert_allclose(
        native.volume(independent, data["temperature_k"]), data["volume_a3"], rtol=1e-10
    )
    for identifier, published in RECORDS.items():
        r = get_eos_record_document(identifier)
        assert list(r["eos"]["parameters"].values()) == published[:3]
        assert r["scientific_validation"]["status"] == "not_reproduced"


def test_reproduction_covariance_multistart_and_ledger_meaning():
    report = reproduce()
    stored = json.loads(REPORT.read_text())
    primary = report["primary"]
    assert primary["solver_success"]
    assert primary["degrees_of_freedom"] == 126
    assert primary["fit_row_metrics"]["row_count"] == 131
    assert primary["active_bounds"] == [0] * 5
    assert primary["excluded_combined_row_numbers"] == []
    np.testing.assert_allclose(
        [primary["parameters"][name] for name in NAMES],
        [22.578319, 150.484868, 5.801501, 1.989345, 0.643085],
        rtol=2e-6,
    )
    assert primary["fit_row_metrics"]["rmse_gpa"] == pytest.approx(0.9844675646)
    constrained = report["sensitivity"]["K0_160_gpa_fixed"]
    assert constrained["parameters"]["K0"] == 160
    assert "K0" not in constrained["free_parameter_names"]
    assert "K0" not in constrained["standard_errors"]
    assert constrained["degrees_of_freedom"] == 127
    assert constrained["fit_row_metrics"]["rmse_gpa"] == pytest.approx(0.9894090061)
    np.testing.assert_allclose(
        [constrained["parameters"][name] for name in NAMES],
        [
            stored["sensitivity"]["K0_160_gpa_fixed"]["parameters"][name]
            for name in NAMES
        ],
        rtol=2e-6,
    )
    for candidate in [primary, *report["multistart"], *report["sensitivity"].values()]:
        assert candidate["solver_success"]
        assert np.all(np.array(candidate["active_bounds"]) == 0)
        assert np.all(np.linalg.eigvalsh(candidate["covariance"]) > 0)
        assert set(candidate["standard_errors"]) == set(
            candidate["free_parameter_names"]
        )
    np.testing.assert_allclose(
        primary["covariance"], stored["primary"]["covariance"], rtol=1e-4, atol=1e-6
    )
    outcome = ledger_outcome(get_eos_record_document(RECORD_ID))
    assert outcome["status"] == "parity"
    assert outcome["fit_kind"] == "stored_refit_reproduction"
    assert outcome["original_publication_reproduction_status"] == "not_reproduced"
    paper = (REPORT.parents[1] / "paper-investigation-ledger.md").read_text()
    assert "Miozzi et al. (2020)" in paper
    assert "Original unreproduced; independent refit available" in paper
