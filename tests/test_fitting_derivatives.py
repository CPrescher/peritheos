"""Matched public fits using analytic and finite-difference native Jacobians."""

import numpy as np
import pytest

from peritheos import _rust
from peritheos.eos.rt import BM3, Vinet
from peritheos.eos.thermal import QCompromiseDebye
from peritheos.fitting import fit_joint_eos, fit_rt_eos, fit_thermal_eos


def compare_fits(analytic, numerical):
    assert analytic.success and numerical.success
    assert analytic.free_parameters == numerical.free_parameters
    assert analytic.degrees_of_freedom == numerical.degrees_of_freedom
    assert analytic.loss == numerical.loss
    assert analytic.parameters == pytest.approx(
        numerical.parameters, rel=2e-5, abs=2e-7
    )
    assert analytic.standard_errors == pytest.approx(
        numerical.standard_errors, rel=2e-5, abs=1e-9
    )
    np.testing.assert_allclose(
        analytic.covariance, numerical.covariance, rtol=5e-5, atol=1e-9
    )
    np.testing.assert_allclose(analytic.correlation, numerical.correlation, atol=2e-5)
    np.testing.assert_allclose(
        analytic.adjusted_volume, numerical.adjusted_volume, atol=2e-7
    )
    if analytic.adjusted_temperature is not None:
        np.testing.assert_allclose(
            analytic.adjusted_temperature, numerical.adjusted_temperature, atol=2e-4
        )
    np.testing.assert_allclose(
        analytic.weighted_residuals, numerical.weighted_residuals, atol=2e-4, rtol=2e-5
    )
    assert analytic.chi_square == pytest.approx(
        numerical.chi_square, rel=1e-8, abs=1e-9
    )
    assert analytic.cost == pytest.approx(numerical.cost, rel=1e-8, abs=1e-9)


def numerical_fit(monkeypatch, function, *args, **kwargs):
    """Force only the private derivative switch; public fit options stay identical."""
    calls = []
    with monkeypatch.context() as patch:
        for name in ("fit_rt_eos_native", "fit_thermal_eos_native"):
            original = getattr(_rust, name)

            def finite_difference(*args, _original=original, **kwargs):
                calls.append(_original)
                return _original(*args, **kwargs, analytic_derivatives=False)

            patch.setattr(_rust, name, finite_difference)
        result = function(*args, **kwargs)
        assert len(calls) == 1, "comparison must exercise the native EOS fitter"
        return result


@pytest.mark.parametrize(
    "errors", ["pressure", "independent", "correlated", "unspecified"]
)
@pytest.mark.parametrize("absolute_sigma", [False, True])
@pytest.mark.parametrize("loss", ["linear", "soft_l1", "huber", "cauchy", "arctan"])
def test_bm3_public_fit_preserves_weighting_and_uncertainty(
    monkeypatch, errors, absolute_sigma, loss
):
    rng = np.random.default_rng(171)
    volume = np.linspace(8.5, 11.3, 24)
    pressure = BM3(11.3, 171.0, 4.7).pressure(volume) + rng.normal(
        0, 0.025, volume.size
    )
    # An outlier exercises robust-loss weighting and its covariance transform.
    pressure[4] += 0.15
    kwargs = dict(
        initial={"K0_prime": 4.1, "V0": 11.25, "K0": 163.0},
        bounds={"K0": (150, 190), "V0": (11, 12), "K0_prime": (3, 6)},
        absolute_sigma=absolute_sigma,
        loss=loss,
        max_nfev=600,
    )
    if errors == "correlated":
        kwargs["observation_covariance"] = [[0.04**2, 0.0001], [0.0001, 0.008**2]]
    elif errors != "unspecified":
        kwargs["pressure_sigma"] = np.linspace(0.03, 0.06, volume.size)
        if errors == "independent":
            kwargs["volume_sigma"] = 0.008
    analytic = fit_rt_eos(BM3, volume, pressure, **kwargs)
    numerical = numerical_fit(monkeypatch, fit_rt_eos, BM3, volume, pressure, **kwargs)
    compare_fits(analytic, numerical)


@pytest.mark.parametrize(
    "errors", ["pressure", "volume", "temperature", "both", "correlated"]
)
@pytest.mark.parametrize("absolute_sigma", [False, True])
@pytest.mark.parametrize("loss", ["linear", "soft_l1"])
def test_q_compromise_joint_fit_preserves_latents_and_profiled_covariance(
    monkeypatch, errors, absolute_sigma, loss
):
    rng = np.random.default_rng(822)
    volume = np.tile(np.linspace(8.8, 11.3, 6), 8)
    temperature = np.repeat(np.linspace(250, 2100, 8), 6)
    expected = QCompromiseDebye(BM3(11.3, 171.0, 4.7), 298.0, 820.0, 1.1, 20.0)
    pressure = expected.pressure(volume, temperature) + rng.normal(0, 0.02, volume.size)
    kwargs = dict(
        initial={
            "theta0": 750.0,
            "rt_eos.K0_prime": 4.3,
            "gamma0": 1.0,
            "rt_eos.V0": 11.25,
            "rt_eos.K0": 163.0,
        },
        fixed={"Tr": 298.0, "n": 20.0},
        bounds={"theta0": (100, 1800), "gamma0": (0.2, 2.0), "rt_eos.V0": (11, 12)},
        absolute_sigma=absolute_sigma,
        loss=loss,
        max_nfev=600,
    )
    if errors == "correlated":
        sigma = np.array([0.04, 0.008, 4.0])
        correlation = np.array([[1.0, 0.25, -0.1], [0.25, 1.0, 0.2], [-0.1, 0.2, 1.0]])
        # Per-point factors exercise nonconstant covariance whitening.
        kwargs["observation_covariance"] = (correlation * np.outer(sigma, sigma))[
            None, :, :
        ] * np.linspace(0.8, 1.2, volume.size)[:, None, None]
    else:
        kwargs["pressure_sigma"] = 0.04
        if errors in ("volume", "both"):
            kwargs["volume_sigma"] = 0.008
        if errors in ("temperature", "both"):
            kwargs["temperature_sigma"] = 4.0
    analytic = fit_joint_eos(
        QCompromiseDebye, BM3, volume, temperature, pressure, **kwargs
    )
    numerical = numerical_fit(
        monkeypatch,
        fit_joint_eos,
        QCompromiseDebye,
        BM3,
        volume,
        temperature,
        pressure,
        **kwargs,
    )
    compare_fits(analytic, numerical)


def test_q_compromise_thermal_only_and_other_reference_fallback(monkeypatch):
    volume = np.linspace(9.0, 11.0, 16)
    temperature = np.linspace(350.0, 1900.0, 16)
    for reference in (BM3(11.3, 171.0, 4.7), Vinet(11.3, 171.0, 4.7)):
        expected = QCompromiseDebye(reference, 298.0, 820.0, 0.0, 20.0)
        pressure = expected.pressure(volume, temperature) + 0.001 * np.sin(temperature)
        kwargs = dict(
            initial={"gamma0": 0.2},
            fixed={"Tr": 298.0, "theta0": 820.0, "n": 20.0},
            pressure_sigma=0.04,
            volume_sigma=0.008,
            temperature_sigma=4.0,
            absolute_sigma=True,
        )
        analytic = fit_thermal_eos(
            QCompromiseDebye, reference, volume, temperature, pressure, **kwargs
        )
        numerical = numerical_fit(
            monkeypatch,
            fit_thermal_eos,
            QCompromiseDebye,
            reference,
            volume,
            temperature,
            pressure,
            **kwargs,
        )
        compare_fits(analytic, numerical)


@pytest.mark.parametrize("correlated", [False, True])
@pytest.mark.parametrize("loss", ["linear", "soft_l1"])
def test_q_compromise_matches_previous_callback_objective(
    monkeypatch, correlated, loss
):
    import peritheos.fitting as fitting

    volume = np.tile(np.linspace(9.0, 11.0, 6), 4)
    temperature = np.repeat([350.0, 700.0, 1200.0, 1800.0], 6)
    reference = BM3(11.3, 171.0, 4.7)
    truth = QCompromiseDebye(reference, 298.0, 820.0, 1.1, 20.0)
    pressure = truth.pressure(volume, temperature) + 0.02 * np.sin(volume * temperature)
    kwargs = dict(
        initial={"theta0": 750.0, "gamma0": 1.0},
        fixed={"Tr": 298.0, "n": 20.0},
        loss=loss,
        absolute_sigma=True,
    )
    if correlated:
        sigma = np.array([0.04, 0.008, 4.0])
        correlation = np.array([[1, 0.25, -0.1], [0.25, 1, 0.2], [-0.1, 0.2, 1]])
        kwargs["observation_covariance"] = correlation * np.outer(sigma, sigma)
    else:
        kwargs.update(pressure_sigma=0.04, volume_sigma=0.008, temperature_sigma=4.0)
    analytic = fit_thermal_eos(
        QCompromiseDebye, reference, volume, temperature, pressure, **kwargs
    )
    monkeypatch.setattr(fitting, "_native_fitting_model", lambda model: None)
    callback = fit_thermal_eos(
        QCompromiseDebye, reference, volume, temperature, pressure, **kwargs
    )
    compare_fits(analytic, callback)


def test_q_compromise_subclasses_keep_custom_pressure_evaluation(monkeypatch):
    class ShiftedQCompromise(QCompromiseDebye):
        def pressure(self, volume, temperature):
            return super().pressure(volume, temperature) + 0.5

    def forbidden(*args, **kwargs):
        raise AssertionError("custom scientific behavior must use callbacks")

    monkeypatch.setattr(_rust, "fit_thermal_eos_native", forbidden)
    volume = np.linspace(9, 11, 12)
    temperature = np.linspace(350, 1800, 12)
    reference = BM3(11.3, 171.0, 4.7)
    model = ShiftedQCompromise(reference, 298.0, 820.0, 1.1, 20.0)
    result = fit_thermal_eos(
        ShiftedQCompromise,
        reference,
        volume,
        temperature,
        model.pressure(volume, temperature),
        initial={"gamma0": 1.0},
        fixed={"Tr": 298.0, "theta0": 820.0, "n": 20.0},
    )
    assert result.success
    assert result.parameters["gamma0"] == pytest.approx(1.1)
