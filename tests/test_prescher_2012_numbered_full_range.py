"""Scientific controls for the fixed-V0 coordinate-error comparison."""

import numpy as np
import pytest

from scripts.check_prescher_2012_numbered_full_range import (
    fixed_coordinate_fit,
    fixed_odr_check,
)
from scripts.recover_prescher_2012_fe3c import volume_at_pressure


def test_fixed_reference_recovers_known_bm3_and_matches_inverse_odr():
    truth = [154.7, 162.4, 5.7]
    p = np.array([1.3, 5.1, 9.7, 17.2, 25.6, 37.8, 47.5])
    v = volume_at_pressure(p, *truth)
    sp, sv = np.full(7, 0.1), np.full(7, 0.07)
    result = fixed_coordinate_fit(p, v, sp, sv, truth[0])
    check = fixed_odr_check(p, v, sp, sv, truth[0])
    assert list(result["parameters"].values()) == pytest.approx(truth, rel=1e-8)
    assert result["fixed_parameters"] == {"V0": truth[0]}
    assert result["degrees_of_freedom"] == 5
    assert result["weighted_sum_squares"] < 1e-16
    assert np.all(np.linalg.eigvalsh(result["absolute_covariance"]) > 0)
    np.testing.assert_allclose(
        result["absolute_covariance"],
        check["absolute_covariance_K0_K0prime"],
        rtol=5e-4,
    )


def test_fixed_reference_error_scale_keeps_residual_covariance_distinct():
    p = np.array([1.3, 5.1, 9.7, 17.2, 25.6, 37.8, 47.5])
    v = volume_at_pressure(p, 154.7, 162.4, 5.7) + np.array(
        [0.03, -0.02, 0.04, -0.04, 0.02, -0.01, 0.03]
    )
    sp, sv = np.full(7, 0.1), np.full(7, 0.07)
    first = fixed_coordinate_fit(p, v, sp, sv, 154.7)
    twice = fixed_coordinate_fit(p, v, 2 * sp, 2 * sv, 154.7)
    assert list(twice["parameters"].values()) == pytest.approx(
        list(first["parameters"].values()), rel=1e-7
    )
    np.testing.assert_allclose(
        twice["absolute_covariance"],
        4 * np.array(first["absolute_covariance"]),
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        twice["residual_scaled_covariance"],
        first["residual_scaled_covariance"],
        rtol=1e-6,
    )
