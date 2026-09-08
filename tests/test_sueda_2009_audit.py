"""Published Sueda thermal fits, complete transcription, and independent checks."""

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.errors import MaterialError
from scripts.reproduce_sueda_2009_mgal2o4 import (
    DATA_PATH,
    STATIC,
    bm3,
    htbm_pressure,
    load_data,
    mgd_pressure,
    reproduce,
)

ROOT = Path(__file__).resolve().parents[1]
DOI = "10.1016/j.pepi.2008.07.046"
PREFIX = "mgal2o4_cafe2o4_sueda_2009_"
CASES = [
    ("htbm_2", htbm_pressure, 0.307385763448),
    ("mgd_3", mgd_pressure, 0.288157430963),
]


def test_table1_transcription_preserves_all_observations_and_missing_errors():
    rows = load_data()
    assert len(rows) == 46
    assert Counter(row["run"] for row in rows) == {
        "M114": 21,
        "M140": 4,
        "M194": 5,
        "CF-1": 5,
        "CF-2": 11,
    }
    assert sum(row["static_refit_included"] == "true" for row in rows) == 16
    assert all(row["thermal_fit_included"] == "true" for row in rows)
    assert sum(bool(row["au_volume_a3"]) for row in rows) == 33
    assert sum(bool(row["ruby_pressure_gpa"]) for row in rows) == 4
    assert rows[28]["volume_a3_uncertainty"] == "2.0"  # Printed 221.8(20).
    assert rows[26]["a_angstrom"] == "8.492"  # Retain printed axes.
    assert rows[30]["pressure_gpa"] == "1.5"
    assert rows[30]["ruby_pressure_gpa"] == "1.6"
    assert rows[35]["temperature_k"] == "836.0"
    assert rows[-1]["temperature_k"] == rows[-2]["temperature_k"] == "298.0"
    assert rows[-1]["volume_a3"] == "239.6"
    assert rows[-1]["pressure_gpa"] == "0.0001"
    assert rows[-1]["pressure_gpa_uncertainty"] == ""
    assert rows[-1]["au_volume_a3"] == ""
    dataset = get_material_document("mgal2o4_cafe2o4")["datasets"][-1]
    assert (
        hashlib.sha256(DATA_PATH.read_bytes()).hexdigest()
        == dataset["resource"]["sha256"]
    )
    assert dataset["used_by_eos_records"] == [PREFIX + case[0] for case in CASES]


@pytest.mark.parametrize("suffix,independent,rmse", CASES)
def test_native_model_matches_independent_cell_normalization_and_table(
    suffix, independent, rmse
):
    document = get_material_document("mgal2o4_cafe2o4")
    material = Material.from_eosmat(document)
    stored = next(
        r for r in document["eos_records"] if r["identifier"] == PREFIX + suffix
    )
    record = material.get_eos_record(PREFIX + suffix)
    assert stored["reference"]["doi"] == DOI
    assert stored["eos"]["parameters"] == STATIC
    assert stored["fixed_parameters"] == ["V0", "K0", "K0_prime"]
    assert stored["parameter_error_confidence"] is None
    assert stored["parameter_covariance"] is None
    assert (
        stored["pressure_calibration"]["methods"][0]["reference_eos_record"]
        == "gold_anderson_1989_bm3_1"
    )
    assert stored["scientific_validation"]["status"] == "primary_source_validated"
    rows = load_data()
    v, t, p = (
        np.array([float(row[k]) for row in rows])
        for k in ("volume_a3", "temperature_k", "pressure_gpa")
    )
    calculated = record.pressure(v, t)
    np.testing.assert_allclose(calculated, independent(v, t), atol=2e-10, rtol=0)
    assert np.sqrt(np.mean((calculated - p) ** 2)) == pytest.approx(rmse, abs=1e-10)
    np.testing.assert_allclose(record.pressure(v, 300.0), bm3(v), atol=2e-10, rtol=0)
    assert record.pressure(240.1, 300.0) == pytest.approx(0.0, abs=1e-10)
    assert record.eos.bulk_modulus(240.1 * record.volume_scale, 300.0) == pytest.approx(
        205.0
    )
    # Independent measured high-P/high-T states; tolerance includes propagated
    # printed P/V errors (roughly 0.3-0.6 GPa), with a conservative two-error margin.
    assert record.pressure(218.9, 1800.0) == pytest.approx(32.4, abs=0.6)
    assert record.pressure(223.9, 2400.0) == pytest.approx(30.1, abs=1.0)
    assert record.volume(0.0001, 836.0) == pytest.approx(243.8, abs=0.3)
    np.testing.assert_allclose(record.volume(calculated, t), v, rtol=2e-9)
    roundtrip = Material.from_eosmat(material.to_eosmat())
    np.testing.assert_allclose(
        roundtrip.get_eos_record(PREFIX + suffix).pressure(v, t), calculated
    )
    with pytest.raises(MaterialError):
        record.pressure(218.9, 2500.0, check_validity=True)
    if suffix == "mgd_3":
        assert stored["thermal"]["parameters"]["n"] == 7
        assert stored["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    else:
        assert stored["thermal"]["thermal_expansion_law"] == "linear_temperature"
        assert stored["thermal"]["parameters"]["alpha0"] == 1.96e-5


def test_independent_staged_refits_recover_published_coefficients():
    result = reproduce()
    assert result["static"]["observations"] == 16
    assert result["static"]["solver_success"]
    for name, error in (("V0", 0.2), ("K0", 6.0), ("K0_prime", 0.3)):
        assert result["static"]["parameters"][name] == pytest.approx(
            STATIC[name], abs=error
        )
    for name in ("htbm", "mgd"):
        assert result[name]["solver_success"]
        assert (
            max(
                abs(x)
                for x in result[name][
                    "parameter_differences_in_reported_errors"
                ].values()
            )
            < 1.0
        )
        assert result[name]["refit_rmse_gpa"] < result[name]["published_rmse_gpa"]


def test_sueda_promoted_once_without_fabricating_imported_candidates():
    def load(path):
        return json.loads((ROOT / path).read_text(encoding="utf-8"))

    assert not any(
        r["doi"].lower() == DOI
        for r in load("docs/data/nonproduction-paper-investigations.json")["papers"]
    )
    assert not any(
        (r["publication"].get("doi") or "").lower() == DOI
        for r in load("docs/data/litcurate-eos-candidates.json")["records"]
    )
    records = [
        r
        for r in load("peritheos/data/primary-source-audit.json")["records"]
        if (r.get("doi") or "").lower() == DOI
    ]
    assert {r["record"] for r in records} == {PREFIX + case[0] for case in CASES}
    document = get_material_document("mgal2o4_cafe2o4")
    assert [r["identifier"] for r in document["eos_records"] if r.get("default")] == [
        "mgal2o4_cafe2o4_funamori_1998_bm2_1"
    ]
