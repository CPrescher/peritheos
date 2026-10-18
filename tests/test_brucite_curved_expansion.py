"""The purple-law variant is a separate ambient refit with adopted moduli."""

import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, eosmat_schema, get_material
from peritheos.errors import MaterialError
from scripts.brucite_curved_expansion import PREFIX, REPORT, ledger_outcome, reproduce
from scripts.brucite_replay_checks import check_report
from scripts.brucite_reproduce import bm3, observations, source
from scripts.brucite_stage import audit_document


def audit_record(identifier):
    return Material.from_eosmat(
        audit_document(), require_primary_validation=False
    ).get_eos_record(identifier)


@pytest.fixture(scope="module")
def report():
    # Runtime checks use the registered coefficients; replay agreement is separate.
    saved = json.loads(REPORT.read_text(encoding="utf-8"))
    check_report(reproduce(), saved)
    return saved


@pytest.mark.parametrize("index,temperature", [(0, 300), (1, 473), (2, 673), (3, 873)])
def test_curved_isotherms_preserve_published_moduli_and_source_ranges(
    report, index, temperature
):
    original = audit_record(f"brucite_fukui_2003_{temperature}k_bm3")
    refit = audit_record(f"{PREFIX}_{temperature}k_bm3_refit")
    meta = refit.eosmat_metadata
    assert not refit.is_default
    assert meta["record_kind"] == "refit"
    for key in ("K0", "K0_prime"):
        assert (
            meta["eos"]["parameters"][key]
            == original.eosmat_metadata["eos"]["parameters"][key]
        )
        assert (
            meta["parameter_errors"][key]
            == original.eosmat_metadata["parameter_errors"][key]
        )
    assert (
        meta["validity"]["pressure_gpa"]
        == original.eosmat_metadata["validity"]["pressure_gpa"]
    )
    rows = [
        r
        for r in observations("fukui_2003", "table4")
        if float(r["temperature_k"]) == temperature
    ]
    residual = [
        refit.pressure(float(r["cell_volume_a3"]), temperature)
        - float(r["pressure_gpa"])
        for r in rows
    ]
    assert np.sqrt(np.mean(np.array(residual) ** 2)) == pytest.approx(
        report["isotherms"][index]["pressure_rmse_gpa"], abs=1e-12
    )
    lower, upper = meta["validity"]["pressure_gpa"]
    for pressure in [lower, 10, upper]:
        v = refit.volume(pressure, temperature, check_validity=True)
        assert refit.pressure(v, temperature) == pytest.approx(pressure)
    with pytest.raises(MaterialError):
        refit.volume(10, temperature + 1, check_validity=True)
    # One propagated V0 variance is conditional on the adopted moduli.
    v = refit.volume(10, temperature)
    expected_sigma = (
        refit.eos.bulk_modulus(v)
        / refit.reference_volume
        * meta["parameter_errors"]["V0"]
    )
    assert refit.pressure_with_uncertainty(
        v, temperature
    ).standard_error == pytest.approx(expected_sigma, rel=1e-5)


def test_curved_thermal_uses_exact_purple_law_and_published_slopes(report):
    r = audit_record(f"{PREFIX}_thermal_refit")
    meta = r.eosmat_metadata
    fit = report["fit"]
    assert not r.is_default
    assert meta["thermal"]["thermal_expansion_law"] == "linear_reference_temperature"
    assert meta["fit_datasets"] == ["brucite_fukui_2003_table1"]
    assert meta["fit_provenance"]["selection"]["source_rows"] == list(range(1, 12))
    assert meta["fit_provenance"]["statistics"]["observations"] == 11
    v, t = np.meshgrid([34.0, 36.0, 40.0], [300.0, 573.0, 873.0])
    delta = t - 300
    v0 = fit["V300_a3"] * np.exp(
        fit["alpha_300_per_k"] * delta + 0.5 * fit["dalpha_dT_per_k2"] * delta**2
    )
    expected = bm3(v, v0, 41.8 - 0.0319 * delta, 6.6 + 0.00571 * delta)
    np.testing.assert_allclose(r.pressure(v, t), expected, rtol=1e-12, atol=1e-12)
    native = r.eos._native
    del r.eos._native
    try:
        np.testing.assert_allclose(r.pressure(v, t), expected, rtol=1e-12, atol=1e-12)
    finally:
        r.eos._native = native
    np.testing.assert_allclose(
        native.evaluate_array("pressure", v.ravel(), t.ravel()),
        expected.ravel(),
        rtol=1e-12,
        atol=1e-12,
    )
    for p, temp in [(5.95, 300), (10, 573), (15.48, 873)]:
        volume = r.volume(p, temp, check_validity=True)
        assert r.pressure(volume, temp) == pytest.approx(p)
    for p, temp in [(5.94, 873), (15.49, 873), (10, 874), (0, 873)]:
        with pytest.raises(MaterialError):
            r.volume(p, temp, check_validity=True)
    # Slopes preserve the printed continuous law; Table 6 slices are distinct.
    assert report["isotherms"][1]["continuous_variant"]["K0_gpa"] != 32.6
    assert report["isotherms"][1]["published_moduli"]["K0_gpa"] == 32.6


def test_curved_expansion_covariance_roundtrip_and_registered_fit(report):
    check_report(report, json.loads(REPORT.read_text(encoding="utf-8")))
    doc = audit_document()
    Draft202012Validator(eosmat_schema()).validate(doc)
    m = Material.from_eosmat(doc, require_primary_validation=False)
    restored = Material.from_eosmat(m.to_eosmat(), require_primary_validation=False)
    r = restored.get_eos_record(f"{PREFIX}_thermal_refit")
    assert r.eosmat_metadata["parameter_covariance"]["parameter_order"] == [
        "rt_eos.V0",
        "alpha0",
        "alpha1",
    ]
    parameters = np.array(
        [
            report["fit"][key]
            for key in ["V300_a3", "alpha_300_per_k", "dalpha_dT_per_k2"]
        ]
    )

    def pressure(x):
        v0 = x[0] * np.exp(x[1] * 573 + 0.5 * x[2] * 573**2)
        return bm3(35.5, v0, 41.8 - 0.0319 * 573, 6.6 + 0.00571 * 573)

    gradient = []
    for index, h in enumerate([1e-4, 1e-9, 1e-12]):
        shift = np.zeros(3)
        shift[index] = h
        gradient.append(
            (pressure(parameters + shift) - pressure(parameters - shift)) / (2 * h)
        )
    covariance = np.array(r.eosmat_metadata["parameter_covariance"]["matrix"])
    expected = np.sqrt(np.array(gradient) @ covariance @ np.array(gradient))
    assert r.pressure_with_uncertainty(35.5, 873).standard_error == pytest.approx(
        expected, rel=1e-5
    )
    assert (
        get_material("brucite").default_record().identifier
        == "brucite_fukui_2003_300k_bm3"
    )
    for record in doc["eos_records"]:
        if "curved_expansion" in record["identifier"]:
            outcome = ledger_outcome(record)
            assert outcome["parity_basis"] == "registered_peritheos_refit_reproduced"
            assert not outcome["original_author_workflow_parity"]
            assert outcome["observations"] == 11
    # Published coefficients remain traceable directly to the same source.
    assert (
        source()["sources"]["fukui_2003"]["parameters"]["isotherms"][-1]["K0_gpa"]
        == 23.9
    )
