"""Primary-scale EOS reproduction, not a reconstruction of the 2015 synthesis."""

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
    thermal_model,
)


def test_explicit_source_masks_and_original_pressure_coordinates():
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
        assert fit["status"] == "similar"
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
    assert all(by_id[k]["status"] == "similar" for k in IDS)
    assert by_id[IDS[3]]["stages"][1]["dependent_parameters"] == {
        "alpha0": "1.7e-5 - 300*alpha1"
    }
