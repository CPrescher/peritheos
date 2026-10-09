"""Author setting and the separately attributed, normalized Speziale Fe refit."""

import hashlib
import json

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import (
    Material,
    get_eos_record,
    get_eos_record_document,
    get_material,
    get_material_document,
    validate_eosmat_document,
)
from scripts.reconstruct_miozzi_2020_iron import (
    NA,
    fe_pressure,
    pressure_column,
    read_data,
)
from scripts.refit_miozzi_2020_speziale import (
    AUTHOR_ID,
    CSV_PATH,
    DATASET_ID,
    MODEL,
    NAMES,
    RECORD_ID,
    REPORT,
    author_record,
    derived_rows,
    ledger_outcome,
    reproduce,
)


def test_author_n2_record_preserves_printed_coefficients_and_communication():
    raw = get_eos_record_document(AUTHOR_ID)
    assert raw["eos"]["parameters"] == {"V0": 22.81, "K0": 129, "K0_prime": 6.24}
    assert raw["thermal"]["parameters"] == {
        "Tr": 300,
        "theta0": 420,
        "gamma0": 1.11,
        "q": 0.3,
        "n": 2,
    }
    assert raw["scientific_validation"]["status"] == "not_reproduced"
    assert (
        raw["scientific_validation"]["personal_communication"]["attribution"]
        == "Miozzi et al."
    )
    assert raw["parameter_covariance"] is None
    native = get_eos_record(AUTHOR_ID)
    assert native.eos.n == 2
    assert native.eos.rt_eos.V0 == pytest.approx(22.81 * NA / 2e25)
    assert author_record(raw) == raw


def test_speziale_derived_pressures_preserve_source_columns_and_holes():
    original, data, hashes = read_data()
    rows, table = derived_rows(original, data)
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    assert CSV_PATH.read_text(encoding="utf-8") == table
    assert report["original_csv_sha256"] == hashes
    assert report["derived_csv_sha256"] == hashlib.sha256(table.encode()).hexdigest()
    for path, fingerprint in report["source_code_sha256"].items():
        from scripts.reconstruct_miozzi_2020_iron import ROOT

        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == fingerprint
    pressure = pressure_column(data, MODEL)[0]
    assert len(rows) == 131
    for i, (source, derived) in enumerate(zip(original, rows)):
        assert derived["printed_pressure_gpa"] == source["pressure_gpa"]
        assert derived["printed_pressure_error_gpa"] == source["pressure_error_gpa"]
        for key in (
            "temperature_k",
            "temperature_error_k",
            "volume_a3",
            "volume_error_a3",
        ):
            assert derived[key] == source[key]
        assert float(derived["pressure_gpa"]) == pytest.approx(pressure[i], abs=5e-10)
        if source["medium"] == "he":
            assert derived["pressure_gpa"] == source["pressure_gpa"]
            assert derived["pressure_temperature_covariance_gpa_k"] == ""
    assert rows[15 + 26]["volume_error_a3"] == ""
    assert rows[15 + 75]["conditional_pressure_error_gpa"] == ""
    ds = next(
        d
        for d in get_material_document("iron")["datasets"]
        if d["identifier"] == DATASET_ID
    )
    assert ds["resource"]["sha256"] == report["derived_csv_sha256"]
    assert ds["pressure_reconstruction"]["model"] == MODEL
    assert "calibration_record" not in ds["pressure_reconstruction"]


def test_speziale_refit_selection_covariance_export_and_native_pressure():
    doc = get_material_document("iron")
    validate_eosmat_document(doc)
    raw = get_eos_record_document(RECORD_ID)
    assert raw["record_kind"] == "refit"
    assert raw["thermal"]["parameters"]["n"] == 1
    assert (
        raw["scientific_validation"]["reproduction_status"]
        == "independent_refit_reproduced"
    )
    assert (
        raw["scientific_validation"]["original_publication_reproduction_status"]
        == "not_reproduced"
    )
    assert (
        get_material("iron").default_record().identifier
        == "iron_sakai_2025_rydberg_stacey_1"
    )
    material = Material.from_eosmat(doc, record_identifiers=[RECORD_ID])
    exported = material.to_eosmat()
    validate_eosmat_document(exported)
    restored = next(r for r in exported["eos_records"] if r["identifier"] == RECORD_ID)
    for key in (
        "parameter_errors",
        "parameter_covariance",
        "fit_provenance",
        "scientific_validation",
    ):
        assert restored[key] == raw[key]
    fit = reproduce()["primary"]
    errors = [fit["conditional_standard_errors"][name] for name in NAMES]
    covariance = np.array(raw["parameter_covariance"]["matrix"])
    np.testing.assert_allclose(np.diag(covariance), np.array(errors) ** 2, rtol=1e-12)
    assert np.all(np.linalg.eigvalsh(covariance) > 0)
    _, data, _ = read_data()
    independent = fe_pressure(
        data["volume_a3"],
        data["temperature_k"],
        [fit["parameters"][name] for name in NAMES],
    )
    native = get_eos_record(RECORD_ID)
    np.testing.assert_allclose(
        native.pressure(data["volume_a3"], data["temperature_k"]),
        independent,
        atol=1e-8,
        rtol=0,
    )
    np.testing.assert_allclose(
        native.volume(independent, data["temperature_k"]), data["volume_a3"], rtol=1e-10
    )
    assert ledger_outcome(raw)["status"] == "parity"


def test_speziale_native_staged_fit_matches_independent_least_squares():
    _, data, _ = read_data()
    pressure = pressure_column(data, MODEL)[0]
    check = least_squares(
        lambda p: fe_pressure(data["volume_a3"], data["temperature_k"], p) - pressure,
        [22.81, 129, 6.24, 1.11, 0.3],
        x_scale=[23, 100, 5, 1, 1],
        bounds=([20, 20, 1, 0.1, -5], [25, 400, 10, 5, 10]),
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    fit = reproduce()["primary"]
    assert check.success
    np.testing.assert_allclose(
        check.x, [fit["parameters"][name] for name in NAMES], rtol=2e-6, atol=2e-6
    )
    assert fit["rmse_gpa"] == pytest.approx(np.sqrt(np.mean(check.fun**2)), abs=1e-9)
    for stage in fit["stages"]:
        assert stage["solver"]["success"]
