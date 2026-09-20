"""Scientific scope and numerical checks for the conditional protocol audit."""

import json

import numpy as np
import pytest

from peritheos import get_eos_record
from scripts.audit_campbell_2009_protocol import (
    OUTPUT,
    calibrant_pressure,
    dataset,
    effective_pressure_variance,
)
from scripts.reconstruct_seagle_2008_fcc import (
    CELL_TO_MOLAR,
    V0_1400,
    calibration_checks,
    completed_rows,
    fcc_pressure,
)
from scripts.reproduce_campbell_2009_buffers import SOURCE, pressure


def test_protocol_pressure_scales_and_missing_errors_remain_distinct():
    data, ids, covariance_pt = dataset("feo", "seagle_2006")
    other, other_ids, other_covariance = dataset("feo", "dewaele_2006_hypothesis")
    assert data.shape == other.shape == (6, 90)
    assert ids == other_ids
    assert np.array_equal(data[:2], other[:2])
    assert np.array_equal(data[2, :25], other[2, :25])
    assert covariance_pt[:25].tolist() == [0] * 25
    assert np.all(covariance_pt[25:] > 0)
    assert np.all(other_covariance[25:] > 0)
    assert [key for key, missing in zip(ids, np.isnan(data[5])) if missing] == [
        "seagle:33",
        "seagle:34",
        "seagle:35",
        "seagle:36",
    ]
    assert np.mean(other[2, 25:] - data[2, 25:]) == pytest.approx(-2.21755, abs=1e-5)
    # Cross-check the audit's independent Dewaele expression against the native EOS.
    native = get_eos_record("iron_dewaele_2006_vinet_thermal")
    for v, t in [(19.2, 1800), (17.8, 2400)]:
        expected = native.pressure(v, t, check_validity=False)
        assert calibrant_pressure(v, t, "dewaele_2006_hypothesis") == pytest.approx(
            expected, abs=1e-7
        )


def test_protocol_matrix_keeps_matched_controls_and_distinct_residual_units():
    report = json.loads(OUTPUT.read_text())
    assert len(report["groups"]) == 12
    assert "hypothetical" in report["qualification"]
    for group in report["groups"].values():
        assert set(group["weighted_row_ids"]) | set(
            group["missing_error_row_ids"]
        ) == set(group["row_ids"])
        assert set(group["weighted_row_ids"]).isdisjoint(group["missing_error_row_ids"])
        assert len(group["fits"]) == 7
        for key, fit in group["fits"].items():
            assert fit["solver_success"]
            ids = (
                group["row_ids"]
                if key.startswith("all_rows")
                else group["weighted_row_ids"]
            )
            assert fit["observations"] == len(ids)
            if fit["weighting"] == "unweighted":
                metric = (
                    "rmse_pressure_gpa"
                    if fit["objective"] == "pressure"
                    else "rmse_volume_cm3_mol"
                )
                assert fit[metric] <= fit["published_" + metric]
                # Ensures the optimizer actually minimizes its declared units,
                # rather than relabelling a pressure fit as a volume fit.
                assert fit["sum_squared_scaled_residuals"] == pytest.approx(
                    fit["observations"] * fit[metric] ** 2, rel=1e-10
                )


def test_fcc_reconstruction_reference_units_and_independent_checkpoints():
    # Basinski's lattice spacings are kX, not angstroms; interpolation precedes
    # cubing. A four-atom cell must be converted to a one-mole Fe volume.
    assert V0_1400 == pytest.approx(7.408637774549367, abs=1e-12)
    assert float(fcc_pressure(V0_1400 / CELL_TO_MOLAR, 1400)) == pytest.approx(
        0, abs=1e-10
    )
    checks = calibration_checks()
    assert len(checks["seagle_checkpoints"]) == 8
    assert all(r["within_reported_uncertainty"] for r in checks["seagle_checkpoints"])
    assert checks["seagle_rmse_gpa"] < 2.85
    # The positive discrepancy is retained rather than silently bias-corrected.
    assert checks["seagle_mean_difference_gpa"] > 2.7
    assert len(checks["funamori_table1"]) == 18
    isotherm = [
        r for r in checks["funamori_table1"] if float(r["temperature_k"]) == 1400
    ]
    assert len(isotherm) == 2
    assert max(abs(r["difference_gpa"]) for r in isotherm) < 0.61


def test_completed_pressures_preserve_missing_data_and_explicit_nominal_choices():
    rows = completed_rows()
    assert len(rows) == 81
    assert [r["source_row"] for r in rows if r["reconstructed_pressure_gpa"] == ""] == [
        "61",
        "62",
    ]
    assert sum(r["status"].startswith("conditional_fcc") for r in rows) == 14
    assert rows[5]["measurement_only_pressure_error_gpa"] == ""
    data, ids, _ = dataset("feo", "seagle_2006", reconstruct_fcc=True)
    nominal, nominal_ids, _ = dataset("feo", "seagle_2006", True, True)
    assert data.shape == (6, 104)
    assert nominal.shape == (6, 106)
    assert set(nominal_ids) - set(ids) == {"seagle:61", "seagle:62"}
    assert [key for key, missing in zip(ids, np.isnan(data[5])) if missing] == [
        "seagle:6",
        "seagle:33",
        "seagle:34",
        "seagle:35",
        "seagle:36",
    ]


def test_shared_fcc_marker_covariance_matches_direct_residual_propagation():
    data, _, cov_pt = dataset("fe_fcc", "funamori_boehler", reconstruct_fcc=True)
    assert data.shape == (6, 29)
    data, cov_pt = data[:, 15:16], cov_pt[15:16]
    v, t, _, sv, st, _ = data
    source = np.array(SOURCE["fe_fcc"])

    def residual(volume, temperature):
        return pressure(volume, temperature, source) - fcc_pressure(
            volume / CELL_TO_MOLAR, temperature
        )

    drdv = (residual(v + 1e-4, t) - residual(v - 1e-4, t)) / 2e-4
    drdt = (residual(v, t + 0.01) - residual(v, t - 0.01)) / 0.02
    direct_variance = (drdv * sv) ** 2 + (drdt * st) ** 2
    cell = v / CELL_TO_MOLAR
    marker_dpdv = (fcc_pressure(cell + 1e-4, t) - fcc_pressure(cell - 1e-4, t)) / (
        2e-4 * CELL_TO_MOLAR
    )
    variance, _ = effective_pressure_variance(
        "fe_fcc", data, cov_pt, marker_dpdv * sv**2
    )
    assert variance == pytest.approx(direct_variance, rel=1e-6)
