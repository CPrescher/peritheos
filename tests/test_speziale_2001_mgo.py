"""Source-equation, normalization, and calibration-convention regression evidence."""

import hashlib
import json

import numpy as np
import pytest

from scripts.reproduce_speziale_2001_mgo import (
    DATA,
    OUTPUT,
    RESIDUALS,
    bm3,
    gamma_theta,
    pressure,
    reproduce,
)


@pytest.fixture(scope="module")
def report():
    return reproduce()[0]


def test_printed_differential_identities_and_constant_q_limit():
    u = np.array([-0.5, -0.3, -0.1, 0.0, 0.02])
    v = 74.71 * np.exp(u)
    h = 1e-5
    gm, thm = gamma_theta(v * np.exp(-h))
    gp, thp = gamma_theta(v * np.exp(h))
    gamma, _ = gamma_theta(v)
    assert (np.log(gp) - np.log(gm)) / (2 * h) == pytest.approx(
        1.65 * np.exp(11.8 * u), rel=1e-8
    )
    assert -(np.log(thp) - np.log(thm)) / (2 * h) == pytest.approx(gamma, rel=1e-8)
    g, th = gamma_theta(v, q1=0)
    assert g == pytest.approx(1.524 * np.exp(1.65 * u))
    assert th == pytest.approx(773 * np.exp((1.524 - g) / 1.65), rel=1e-12)
    g, th = gamma_theta(v, q0=0)
    assert g == pytest.approx(np.full(5, 1.524))
    assert th == pytest.approx(773 * np.exp(-1.524 * u))


def test_reference_temperature_and_equivalent_molar_cell_normalizations(report):
    v = np.array([48.0, 60.0, 74.71])
    t = np.array([3663.0, 1500.0, 300.0])
    assert pressure(v, 300) == pytest.approx(bm3(v), abs=1e-12)
    assert pressure(74.71, 300) == pytest.approx(0, abs=1e-12)
    assert pressure(v, t, atoms=8, formula_units=1) == pytest.approx(
        pressure(v, t), abs=1e-12
    )
    assert report["independent_quadrature"]["max_difference_gpa"] < 1e-9
    assert report["existing_Au_diagnostic_max_difference_gpa"] < 1e-9
    for v, t in [(0, 300), (74.71, 0), (np.nan, 300), (74.71, np.inf)]:
        with pytest.raises(ValueError):
            pressure(v, t)


@pytest.mark.parametrize(
    "law", ["integrated", "eq4_local_q", "eq4_q0", "power", "fixed"]
)
def test_diagnostic_theta_laws_support_scalar_and_array_inputs(law):
    assert pressure(60.0, 2000.0, theta_law=law) == pytest.approx(
        pressure([60.0], [2000.0], theta_law=law)[0]
    )


def test_original_constraints_are_partial_not_a_global_reproduction(report):
    assert report["cold_table1"]["rows"] == 32
    assert report["cold_table1"]["V0"] == pytest.approx(74.71, abs=0.01)
    assert report["cold_table1"]["K0_prime"] == pytest.approx(3.99, abs=0.01)
    assert report["Dewaele2000_hot41"]["q"] == pytest.approx(1.3, abs=0.5)
    assert report["Dewaele2000_hot41"][
        "eq11_mean_abs_relative_volume_difference_percent"
    ] == pytest.approx(0.382298944, abs=1e-7)
    assert report["Fiquet1999_COD_hot34"]["q"] == pytest.approx(0.60775879, abs=1e-6)
    assert report["exact_author_thermal_fit_reproduced"] is False
    assert report["exact_later_pressure_reducer_recovered"] is False
    assert report["published_coefficients_changed"] is False


def test_later_pressure_convention_is_diagnostic_not_speziale_eq11(report):
    c = report["calibration_conventions"]
    assert c["eq11_integrated_300K"]["Fei2004"]["rmse_gpa"] == pytest.approx(
        0.2546472496
    )
    assert c["local_power_integrated_300K"]["Fei2004"]["rmse_gpa"] < 0.01
    assert c["local_power_integrated_300K"]["Hirose2006"]["max_abs_residual_gpa"] < 0.02
    assert c["local_power_integrated_298K"]["Fei2004"]["max_abs_residual_gpa"] < 0.005
    # Distinguish the inferred reducer from the printed monotonic gamma law.
    assert gamma_theta(74.71 * 0.64)[0] == pytest.approx(1.32608359, abs=1e-8)
    assert gamma_theta(74.71 * 0.64, gamma_law="local_power")[0] > 1.51
    assert report["rounding"]["original_eq11_outside_bounds"] == 28


def test_recovered_shock_temperatures_are_not_model_temperature_targets():
    source = json.loads((DATA / "mgo-svendsen-1987-source.json").read_text())
    assert source["model_columns_are_observations"] is False
    assert source["uncertainties"]["covariance"] is None
    assert (
        source["csv_sha256"]
        == hashlib.sha256(
            (DATA / "mgo-svendsen-1987-table6.csv").read_bytes()
        ).hexdigest()
    )


def test_archived_report_and_row_predictions_reproduce(report):
    from scripts.reproduce_fei_2007_gold import check_saved

    check_saved(json.loads(OUTPUT.read_text()), report)
    assert RESIDUALS.read_text() == reproduce()[1]
