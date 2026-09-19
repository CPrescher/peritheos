"""Verify independent B2 refits, conditional covariance and source separation."""

import json

import numpy as np
import pytest
from scipy.constants import Avogadro

from peritheos import Material, get_material_document
from scripts.refit_fischer_2014_b2 import REPORT_PATH, reproduce
from scripts.reproduce_fischer_2014_fesi import load_data, parameters, pressure
from scripts.validate_primary_eos_refits import _fit_record


@pytest.fixture(scope="module")
def refit_evidence():
    return reproduce()


@pytest.mark.parametrize("model", ("bm3", "vinet"))
def test_independent_refit_covariance_and_sensitivities_are_reproducible(
    model, refit_evidence, assert_audit_close
):
    saved = json.loads(REPORT_PATH.read_text())
    assert_audit_close(refit_evidence, saved, rel=2e-4, abs=2e-6)
    evidence = refit_evidence["models"][model]
    fit = evidence["fit"]
    assert fit["observations"] == 114
    assert fit["degrees_of_freedom"] == 111
    assert fit["jacobian_rank"] == 3
    assert fit["active_bounds"] == [0, 0, 0]
    assert 2.3 < fit["pressure_rmse_gpa"] < 2.5
    assert np.linalg.eigvalsh(fit["covariance"]).min() > 0
    assert fit["correlation"][1][2] > 0.85
    for result in evidence["multistart"]:
        np.testing.assert_allclose(
            list(result["parameters"].values()),
            list(fit["parameters"].values()),
            rtol=2e-5,
        )
    # Real source-scope sensitivities must remain visible, not hidden in a tiny
    # local standard error. These are not probabilities or independent runs.
    weighted = evidence["pressure_sigma_only_sensitivity"]
    assert abs(weighted["parameters"]["gamma0"] - fit["parameters"]["gamma0"]) > 0.25
    blocks = evidence["pressure_block_holdouts"]
    assert sorted(
        row for block in blocks for row in block["held_out_source_rows"]
    ) == sorted(saved["source_rows"])
    assert all(block["observations"] == 76 for block in blocks)
    assert all(block["held_out_pressure_rmse_gpa"] > 3 for block in blocks)


@pytest.mark.parametrize("index,model", ((0, "bm3"), (1, "vinet")))
def test_b2_refit_native_curve_inverse_and_full_covariance(index, model):
    document = get_material_document("fesi_b2")
    stored = document["eos_records"][index]
    material = Material.from_eosmat(document)
    record = material.eos_records[index]
    rows = load_data("fesi_b2")
    cell_v = np.array([float(row["b2_volume_a3"]) for row in rows])
    atomic_v = cell_v * Avogadro / 1e24 / 2
    t = np.array([float(row["temperature_k"]) for row in rows])
    coeff = parameters("fesi_b2", model)
    coeff[[1, 4, 5]] = [
        stored["eos"]["parameters"]["K0"],
        stored["thermal"]["parameters"]["gamma0"],
        stored["thermal"]["parameters"]["q"],
    ]
    expected = pressure("fesi_b2", atomic_v, t, model, coeff)
    np.testing.assert_allclose(record.pressure(cell_v, t), expected, atol=1e-9, rtol=0)
    np.testing.assert_allclose(record.volume(expected, t), cell_v, atol=2e-8, rtol=0)
    assert record.pressure(record.reference_volume, 300) == pytest.approx(0, abs=1e-10)
    assert stored["thermal"]["parameters"]["n"] == 2
    # Independently differentiate the physical-atom implementation to verify
    # public cell-volume conversion and the cross-reference K0 covariance key.
    volume, temperature = 17.0, 2000.0
    derivative = []
    for i in (1, 4, 5):
        step = 1e-4 * max(1.0, abs(coeff[i]))
        plus, minus = coeff.copy(), coeff.copy()
        plus[i] += step
        minus[i] -= step
        derivative.append(
            (
                pressure(
                    "fesi_b2", volume * Avogadro / 1e24 / 2, temperature, model, plus
                )
                - pressure(
                    "fesi_b2", volume * Avogadro / 1e24 / 2, temperature, model, minus
                )
            )
            / (2 * step)
        )
    derivative = np.asarray(derivative)
    covariance = np.asarray(stored["parameter_covariance"]["matrix"])
    expected_error = np.sqrt(derivative @ covariance @ derivative)
    prediction = record.pressure_with_uncertainty(volume, temperature)
    assert prediction.standard_error == pytest.approx(expected_error, rel=1e-6)
    assert 0.2 < prediction.standard_error < 0.5
    with pytest.raises(ValueError, match="validity|calibration"):
        record.volume(330, 2000, check_validity=True)
    result = _fit_record(document, stored, document["datasets"][0])
    assert result["status"] == "parity"
    assert result["observations"] == 114
    assert "not Fischer's published B2 parameters" in result["qualification"]
    assert (
        Material.from_eosmat(material.to_eosmat())
        .eos_records[index]
        .parameter_covariance
        == record.parameter_covariance
    )


def test_refit_only_b2_has_real_structure_and_no_fictitious_published_parent():
    document = get_material_document("fesi_b2")
    assert document["space_group_number"] == 221
    assert document["formula_units_per_cell"] == 1
    assert [
        (s["element"], s["wyckoff"], s["occupancy"]) for s in document["atom_sites"]
    ] == [("Fe", "1a", 1.0), ("Si", "1b", 1.0)]
    assert document["lattice"]["a"] == pytest.approx(2.666351276708074)
    assert "42.7476826557 GPa" in document["source"]["structure_location"]
    for record in document["eos_records"]:
        assert record["record_kind"] == "refit"
        assert record["identifier"].endswith("_refit")
        assert record["label"].startswith("Peritheos")
        assert "derived_from_record" not in record
        assert record["parameter_covariance"]["parameter_order"] == [
            "rt_eos.K0",
            "gamma0",
            "q",
        ]
        assert record["fixed_parameters"] == ["V0", "K0_prime"]
        assert (
            record["pressure_calibration"]["recalculation"]["status"]
            == "reference_eos_not_bundled"
        )
