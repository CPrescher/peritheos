"""Independent scientific checks of coordinate-error fit conventions."""

import numpy as np
import pytest

from scripts.check_prescher_2012_coordinate_errors import coordinate_fit
from scripts.recover_prescher_2012_fe3c import volume_at_pressure


def test_coordinate_fit_recovers_known_bm3_with_ambient_observation():
    truth = [155.2, 160.0, 6.0]
    p = np.array([1.0, 3.0, 7.0, 12.0, 20.0, 30.0, 45.0])
    v = volume_at_pressure(p, *truth)
    result = coordinate_fit(p, v, np.full(7, 0.1), np.full(7, 0.07), (truth[0], 0.19))
    assert list(result["parameters"].values()) == pytest.approx(truth, rel=1e-7)
    assert result["weighted_sum_squares"] < 1e-16
    assert result["degrees_of_freedom"] == 5
    assert np.all(np.linalg.eigvalsh(result["absolute_covariance"]) > 0)


def test_scaling_all_coordinate_errors_changes_only_absolute_error_scale():
    p = np.array([1.0, 3.0, 7.0, 12.0, 20.0, 30.0, 45.0])
    v = volume_at_pressure(p, 155.2, 160, 6) + np.array(
        [0.03, -0.02, 0.04, -0.04, 0.02, -0.01, 0.03]
    )
    sp, sv = np.full(7, 0.1), np.full(7, 0.07)
    first = coordinate_fit(p, v, sp, sv, (155.2, 0.19))
    twice = coordinate_fit(p, v, 2 * sp, 2 * sv, (155.2, 0.38))
    assert list(twice["parameters"].values()) == pytest.approx(
        list(first["parameters"].values()), rel=1e-6
    )
    np.testing.assert_allclose(
        twice["absolute_covariance"],
        4 * np.array(first["absolute_covariance"]),
        rtol=1e-5,
    )
    np.testing.assert_allclose(
        twice["residual_scaled_covariance"],
        first["residual_scaled_covariance"],
        rtol=1e-5,
    )
    assert twice["weighted_sum_squares"] == pytest.approx(
        first["weighted_sum_squares"] / 4
    )


def test_exact_ambient_pressure_does_not_fix_reference_volume():
    p = np.array([1.0, 3.0, 7.0, 12.0, 20.0, 30.0, 45.0])
    v = volume_at_pressure(p, 155.2, 160, 6)
    result = coordinate_fit(p, v, np.full(7, 0.1), np.full(7, 0.07), (155.7, 0.19))
    assert 155.2 < result["parameters"]["V0"] < 155.69
    assert result["absolute_errors_V0_K0_K0prime"][0] > 0
