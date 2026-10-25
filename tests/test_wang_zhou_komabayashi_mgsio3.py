"""Primary-scale EOS reproduction, not a reconstruction of the 2015 synthesis."""

import copy
import hashlib
import json

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.eosmat import validate_eosmat_document
from scripts.reproduce_wang_zhou_komabayashi_mgsio3 import (
    DATA,
    IDS,
    REPORT,
    observations,
    reproduce,
    static_pressure_uncertainty,
    thermal_model,
)


def test_explicit_source_masks_and_original_pressure_coordinates():
    assert (
        hashlib.sha256(
            (DATA / "akimotoite-wang-2004-tables1-2-pvt.csv").read_bytes()
        ).hexdigest()
        == "e242d075f32537b336528ef2fb9280847a418c931d0808a3117e51f305cdd04b"
    )
    data = observations()
    assert [len(data[k]) for k in IDS] == [13, 10, 17, 21]
    assert min(r[2] for r in data[IDS[1]]) == -0.05
    assert max(r[2] for r in data[IDS[1]]) == 17.044
    assert sum(r[2] == 0 for r in data[IDS[2]]) == 3
    assert sum(r[3] == 300 for r in data[IDS[3]]) == 3
    assert 18 not in [r[0] for r in data[IDS[3]]]  # Pv-only observation
    assert data[IDS[2]][0][1] == pytest.approx(
        6 * 100.387 * 1e24 / 6.02214076e23 / 3.8099
    )


def test_refits_recover_published_coefficients_without_claiming_exact_parity():
    report = reproduce()
    stored = json.loads(REPORT.read_text())
    for filename, digest in report["source_sha256"].items():
        assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == digest
    expected = [
        (264.3240246, 4.4577030),
        (263.9298820, 5.5689043),
        (262.4091655, 206.7705206),
        (223.1705607, 1.0969957e-9, -0.0084155966),
    ]
    for identifier, coefficients in zip(IDS, expected):
        fit = report["fits"][identifier]
        assert fit["status"] == ("parity" if identifier == IDS[1] else "similar")
        assert fit["solver_success"]
        assert [p["refit"] for p in fit["parameters"]] == pytest.approx(
            coefficients, rel=2e-5, abs=1e-14
        )
        assert all(p["within_reported_error"] for p in fit["parameters"])
        assert all(p["within_combined_2sigma"] is None for p in fit["parameters"])
        assert fit["source_orders"] == stored["fits"][identifier]["source_orders"]
        assert [p["refit"] for p in fit["parameters"]] == pytest.approx(
            [p["refit"] for p in stored["fits"][identifier]["parameters"]],
            rel=2e-5,
            abs=1e-14,
        )
    ppv = report["fits"][IDS[3]]
    assert [s["observations"] for s in ppv["stages"]] == [3, 21]
    assert ppv["derived_alpha0"] + 300 * ppv["parameters"][1]["refit"] == pytest.approx(
        1.7e-5
    )
    assert ppv["rounded_K0_sensitivity"]["dK_dT"] == pytest.approx(
        -0.00851051, rel=2e-5
    )


@pytest.mark.parametrize("identifier", IDS)
def test_published_records_validate_and_preserve_coefficients(identifier):
    material_id = (
        "akimotoite"
        if identifier.startswith("akimotoite")
        else "mgsio3_post_perovskite"
    )
    document = get_material_document(material_id)
    validate_eosmat_document(document)
    record = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    fit = json.loads(REPORT.read_text())["fits"][identifier]
    params = {
        **record["eos"]["parameters"],
        **record.get("thermal", {}).get("parameters", {}),
    }
    errors = {
        **record["parameter_errors"],
        **record.get("thermal", {}).get("parameter_errors", {}),
    }
    for comparison in fit["parameters"]:
        name = comparison["parameter"]
        assert params[name] == comparison["published"]
        assert errors[name] == comparison["published_error"]
    assert not record.get("default", False)
    assert (
        record["experimental_pressure_range_gpa"] == fit["observed_pressure_range_gpa"]
    )
    dataset = next(
        d for d in document["datasets"] if d["identifier"] == record["fit_datasets"][0]
    )
    assert identifier in dataset["used_by_eos_records"]
    loaded = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).eos_records[0]
    source = np.array(observations()[identifier])
    computed = loaded.pressure(source[:, 1], source[:, 3], check_validity=False)
    assert np.sqrt(np.mean((computed - source[:, 2]) ** 2)) == pytest.approx(
        fit["published_rmse_gpa"]
    )


def test_thermal_constraint_and_equation_are_not_constant_alpha_or_berman():
    k0, alpha1, dkdt = 223.2, 1.13e-9, -0.0085
    temperature, volume = np.array([300, 1000, 2330]), np.array([121, 123, 125])
    alpha0 = 1.7e-5 - 300 * alpha1
    vt = 163.813 * np.exp(
        alpha0 * (temperature - 300) + 0.5 * alpha1 * (temperature**2 - 300**2)
    )
    kt = k0 + dkdt * (temperature - 300)
    ratio = vt / volume
    # K'=4 removes the third-order correction term, leaving this BM expression.
    expected = 1.5 * kt * (ratio ** (7 / 3) - ratio ** (5 / 3))
    assert thermal_model(k0, alpha1, dkdt).pressure(
        volume, temperature
    ) == pytest.approx(expected)


def test_common_ledger_has_four_independent_records_with_source_constraints():
    ledger = json.loads((REPORT.parent / "primary-eos-refits.json").read_text())
    by_id = {r["record_identifier"]: r for r in ledger["records"]}
    assert all(by_id[k]["status"] == "parity" for k in IDS)
    assert by_id[IDS[3]]["stages"][1]["dependent_parameters"] == {
        "alpha0": "1.7e-5 - 300*alpha1"
    }


def test_t0150_covariance_matches_independent_bm3_finite_differences():
    """Derive the covariance from a direct expression, outside the library EOS."""
    fit = reproduce()["fits"][IDS[1]]
    volume, pressure = np.array(observations()[IDS[1]])[:, 1:3].T
    values = np.array([p["refit"] for p in fit["parameters"]])

    def pressure_direct(params):
        ratio = params[0] / volume
        return (
            1.5
            * 210
            * (ratio ** (7 / 3) - ratio ** (5 / 3))
            * (1 + 0.75 * (params[1] - 4) * (ratio ** (2 / 3) - 1))
        )

    jacobian = []
    for index, step in enumerate((1e-3, 1e-3)):
        delta = np.zeros(2)
        delta[index] = step
        jacobian.append(
            (pressure_direct(values + delta) - pressure_direct(values - delta))
            / (2 * step)
        )
    jacobian = np.array(jacobian).T
    residual = pressure_direct(values) - pressure
    expected = np.linalg.inv(jacobian.T @ jacobian) * (residual @ residual) / 8
    covariance = np.array(fit["parameter_covariance"])
    np.testing.assert_allclose(covariance, expected, rtol=2e-7)
    errors = np.sqrt(np.diag(covariance))
    assert errors == pytest.approx([0.2280527, 0.75821594], rel=2e-6)
    assert [p["refit_error"] for p in fit["parameters"]] == pytest.approx(errors)
    assert np.linalg.eigvalsh(covariance).min() > 0
    assert covariance[0, 1] / np.prod(errors) == pytest.approx(-0.83905425)
    assert fit["degrees_of_freedom"] == 8
    assert fit["parity_basis"] == "within_reported_parameter_errors"
    assert all(p["within_combined_2sigma"] is None for p in fit["parameters"])
    # The variance estimate must follow residual scatter, not an assumed zero error.
    doubled = static_pressure_uncertainty(volume, values, 2 * residual)
    np.testing.assert_allclose(doubled["matrix"], 4 * covariance)
    with pytest.raises(ValueError, match="n > 2"):
        static_pressure_uncertainty(volume[:2], values, residual[:2])
    with pytest.raises(ValueError, match="rank-deficient"):
        static_pressure_uncertainty(np.repeat(250.0, 10), values, residual)


def test_t0150_uncertainty_assessment_round_trips_without_replacing_source(
    assert_audit_close,
):
    from peritheos import convert_eosmat, load_eosmat
    from scripts.generate_reconstruction_reports import ROOT, explicit_reports
    from scripts.standardize_reconstruction_metadata import assessment_for_record
    from scripts.validate_primary_eos_refits import _mgsio3_primary_outcome

    fit = reproduce()["fits"][IDS[1]]
    ledger = json.loads((REPORT.parent / "primary-eos-refits.json").read_text())
    row = next(r for r in ledger["records"] if r["record_identifier"] == IDS[1])
    native = load_eosmat(DATA.parent / "materials/akimotoite.eosmat")
    record = next(r for r in native["eos_records"] if r["identifier"] == IDS[1])
    translated = _mgsio3_primary_outcome(record, fit)
    saved = copy.deepcopy({key: row[key] for key in translated})
    replay = copy.deepcopy(translated)
    # Cancellation magnifies relative solver roundoff in small differences;
    # compare the coefficients and uncertainties, then check each stored identity.
    for outcome in (saved, replay):
        for parameter in outcome["parameters"]:
            difference = parameter.pop("difference")
            assert difference == pytest.approx(
                parameter["refit"] - parameter["published"], abs=1e-12
            )
            relative = parameter.pop("relative_difference")
            assert relative == pytest.approx(
                abs(difference / parameter["published"]), abs=1e-12
            )
    assert_audit_close(saved, replay)
    expected = assessment_for_record(record, row)
    assert record["evidence"]["reproduction"]["assessment"] == expected
    assert expected["scope"] == "parameter_agreement"
    assert expected["status"] == "parity"
    assert expected["fit"]["covariance"]["degrees_of_freedom"] == 8
    assert (
        "literature-reproductions/wang-zhou-komabayashi-mgsio3.md#t0150-uncertainty-reconstruction"
        in explicit_reports(record, ROOT)
    )
    assert (
        "not" in expected["basis"]
    )  # Printed widths do not imply a confidence convention.
    assert record["parameters"] == {
        "V0": {"value": 263.9, "error": 0.2},
        "K0": {"value": 210, "error": None, "fixed": True},
        "K0_prime": {"value": 5.6, "error": 0.8},
    }
    for version in (3, 4):
        converted = convert_eosmat(native, version=version)
        validate_eosmat_document(converted)
        restored = convert_eosmat(converted, version=4)
        restored_record = next(
            r for r in restored["eos_records"] if r["identifier"] == IDS[1]
        )
        assert restored_record["evidence"]["reproduction"]["assessment"] == expected
    diagnostic = fit["covariance"]["volume_weighting_sensitivity"]
    assert diagnostic["parameters"]["K0_prime"] == pytest.approx(6.14407, rel=2e-5)
    assert diagnostic["weighted_sum_squares_per_dof"] > 40
    invalid = copy.deepcopy(fit)
    invalid["parameters"][0]["refit_error"] = None
    with pytest.raises(AssertionError, match="uncertainty audit"):
        _mgsio3_primary_outcome(record, invalid)
