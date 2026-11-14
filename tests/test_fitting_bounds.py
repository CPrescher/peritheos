"""Bounds must constrain free parameters across the public fitting APIs."""

import numpy as np
import pytest

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import LinearThermalPressure
from peritheos.errors import FitValidationError
from peritheos.fitting import (
    fit_acoustic_finite_strain,
    fit_joint_eos,
    fit_linear_us_up,
    fit_rt_eos,
    fit_thermal_eos,
)


@pytest.mark.parametrize("kind", ["rt", "thermal", "joint"])
@pytest.mark.parametrize("bounded", ["typo", "fixed"])
def test_eos_fit_rejects_ignored_bounds(kind, bounded):
    reference = BM3(10, 100, 4)
    volume = np.linspace(8, 10, 10)
    temperature = np.linspace(300, 1000, 10)
    if kind == "rt":

        def call(**kw):
            return fit_rt_eos(
                BM3,
                volume,
                reference.pressure(volume),
                initial={"K0": 100},
                fixed={"V0": 10, "K0_prime": 4},
                **kw,
            )

        fixed_name = "V0"
    elif kind == "thermal":
        eos = LinearThermalPressure(reference, 300, 0.01)

        def call(**kw):
            return fit_thermal_eos(
                LinearThermalPressure,
                reference,
                volume,
                temperature,
                eos.pressure(volume, temperature),
                initial={"alpha_KT": 0.01},
                fixed={"Tr": 300},
                **kw,
            )

        fixed_name = "Tr"
    else:
        eos = LinearThermalPressure(reference, 300, 0.01)

        def call(**kw):
            return fit_joint_eos(
                LinearThermalPressure,
                BM3,
                volume,
                temperature,
                eos.pressure(volume, temperature),
                initial={"alpha_KT": 0.01},
                fixed={
                    "rt_eos.V0": 10,
                    "rt_eos.K0": 100,
                    "rt_eos.K0_prime": 4,
                    "Tr": 300,
                },
                **kw,
            )

        fixed_name = "rt_eos.V0"
    with pytest.raises(FitValidationError, match="Unknown bounded|fixed parameters"):
        call(bounds={fixed_name if bounded == "fixed" else "typo": (1, 20)})


def test_acoustic_and_hugoniot_fits_reject_fixed_bounds_without_rejecting_defaults():
    from peritheos.acoustics import EulerianFiniteStrainAcoustic

    model = EulerianFiniteStrainAcoustic(4.11, 247, 4.5, 176, 1.6)
    density = np.linspace(4.13, 4.5, 10)
    vp, vs = model.velocities(density)
    with pytest.raises(FitValidationError, match="fixed parameters"):
        fit_acoustic_finite_strain(
            density,
            vp,
            vs,
            rho0=4.11,
            initial={"K_S0": 240, "K_S0_prime": 4, "G0_prime": 1.5},
            fixed={"G0": 176},
            bounds={"G0": (100, 200)},
        )
    particle = np.linspace(0, 2, 10)
    shock = 4 + 1.5 * particle
    options = {"initial": {"c0": 3.9}, "fixed": {"s": 1.5}}
    with pytest.raises(FitValidationError, match="fixed parameters"):
        fit_linear_us_up(
            particle, shock, V0=10, rho0=5, bounds={"s": (1, 2)}, **options
        )
    result = fit_linear_us_up(particle, shock, V0=10, rho0=5, **options)
    assert result.success
    assert result.parameters["c0"] == pytest.approx(4)
