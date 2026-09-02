"""Tests for the generic Vinet/double-Debye Helmholtz EOS."""

import numpy as np
import pytest
from scipy.constants import Avogadro, R, electron_volt

from peritheos.eos.rt import BM3, Vinet
from peritheos.eos.thermal import (
    DoubleDebyeHelmholtz,
    DoubleDebyeLogMomentHelmholtz,
)
from peritheos.fitting import fit_thermal_eos

# One A^3/atom expressed as Peritheos' J bar^-1 mol^-1 thermal volume.
ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR = Avogadro * 1.0e-25
EV_PER_ATOM_TO_J_PER_MOL = electron_volt * Avogadro


@pytest.fixture
def diamond_eos():
    """Benedict et al. (2014) Table I diamond parameter set."""
    volume_factor = ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    cold_curve = Vinet(
        V0=5.7034 * volume_factor,
        K0=432.4,
        K0_prime=3.793,
    )
    return DoubleDebyeHelmholtz(
        rt_eos=cold_curve,
        Vp=5.571 * volume_factor,
        theta_a0=1887.8,
        a_a=-0.316 / volume_factor,
        b_a=0.913,
        theta_b0=1887.8,
        a_b=0.168 / volume_factor,
        b_b=0.429,
        theta_1_0=1887.8,
        a_1=0.0846 / volume_factor,
        b_1=0.499,
        n=1.0,
        alpha0=3.79e-5,
        Ve=5.785 * volume_factor,
        kappa=0.0,
        phi0=-9.066 * EV_PER_ATOM_TO_J_PER_MOL,
    )


@pytest.fixture
def correa_diamond_eos():
    """Correa et al. (2008) Table I diamond parameter set."""
    volume_factor = ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    return DoubleDebyeLogMomentHelmholtz(
        rt_eos=Vinet(
            V0=5.785 * volume_factor,
            K0=368.2,
            K0_prime=4.038,
        ),
        Vp=5.571 * volume_factor,
        theta_a0=1887.8,
        a_a=-0.316 / volume_factor,
        b_a=0.913,
        theta_b0=1887.8,
        a_b=0.168 / volume_factor,
        b_b=0.429,
        theta_0_0=1887.8,
        a_0=0.131 / volume_factor,
        b_0=0.202,
        n=1.0,
        anharmonic_a=3.8e-5,
        phi0=-155.059 * EV_PER_ATOM_TO_J_PER_MOL,
    )


def test_table_i_volume_and_energy_conversions(diamond_eos):
    conversion = ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR

    assert conversion == pytest.approx(0.0602214076)
    assert diamond_eos.rt_eos.V0 == pytest.approx(5.7034 * conversion)
    assert diamond_eos.phi0 == pytest.approx(-9.066 * EV_PER_ATOM_TO_J_PER_MOL)
    assert diamond_eos.cold_energy(diamond_eos.rt_eos.V0) == pytest.approx(
        diamond_eos.phi0
    )


def test_temperature_law_and_weights_follow_published_equations(diamond_eos):
    volume = 5.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    theta_a, theta_b, theta_1 = diamond_eos.debye_temperatures(volume)
    weight_a, weight_b = diamond_eos.double_debye_weights(volume)

    expected_theta_a = 1887.8 * (5.0 / 5.571) ** -0.913 * np.exp(-0.316 * (5.571 - 5.0))
    assert theta_a == pytest.approx(expected_theta_a)
    assert weight_a == pytest.approx((theta_b - theta_1) / (theta_b - theta_a))
    assert weight_a + weight_b == pytest.approx(1.0)
    assert weight_a * theta_a + weight_b * theta_b == pytest.approx(theta_1)


def test_degenerate_reference_temperatures_have_finite_limiting_weights(diamond_eos):
    weight_a, weight_b = diamond_eos.double_debye_weights(diamond_eos.Vp)
    pressure = diamond_eos.pressure(diamond_eos.Vp, 300.0)

    assert weight_a == pytest.approx(0.1783709190711836)
    assert weight_a + weight_b == pytest.approx(1.0)
    assert np.isfinite(pressure)


def test_zero_point_energy_is_included(diamond_eos):
    volume = diamond_eos.rt_eos.V0
    theta_a, theta_b, _ = diamond_eos.debye_temperatures(volume)
    weight_a, weight_b = diamond_eos.double_debye_weights(volume)
    expected = 9.0 * R / 8.0 * (weight_a * theta_a + weight_b * theta_b)

    assert diamond_eos.zero_point_energy(volume) == pytest.approx(expected)
    assert diamond_eos.ion_helmholtz_free_energy(volume, 0.0) == pytest.approx(expected)
    assert diamond_eos.molar_heat_capacity_v(volume, 0.0) == 0.0


@pytest.mark.parametrize("temperature", [0.0, 1.0e-8, 1.0, 300.0, 5000.0])
def test_near_zero_and_wide_temperature_range_is_finite(diamond_eos, temperature):
    volume = 0.9 * diamond_eos.rt_eos.V0

    assert np.isfinite(diamond_eos.helmholtz_free_energy(volume, temperature))
    assert np.isfinite(diamond_eos.pressure(volume, temperature))
    assert np.isfinite(diamond_eos.molar_heat_capacity_v(volume, temperature))


def test_free_energy_is_sum_of_individual_contributions(diamond_eos):
    volume = 5.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    temperature = 2000.0
    expected = (
        diamond_eos.cold_energy(volume)
        + diamond_eos.ion_helmholtz_free_energy(volume, temperature)
        + diamond_eos.anharmonic_helmholtz_free_energy(volume, temperature)
    )

    assert diamond_eos.helmholtz_free_energy(volume, temperature) == pytest.approx(
        expected
    )
    assert diamond_eos.pressure(volume, temperature) == pytest.approx(
        diamond_eos.rt_eos.pressure(volume)
        + diamond_eos.ion_pressure(volume, temperature)
        + diamond_eos.anharmonic_pressure(volume, temperature)
    )


@pytest.mark.parametrize(
    "atomic_volume,temperature",
    [(5.7034, 0.0), (5.571, 300.0), (5.0, 2000.0), (4.5, 5000.0)],
)
def test_analytic_pressure_equals_minus_volume_derivative_of_free_energy(
    diamond_eos, atomic_volume, temperature
):
    volume = atomic_volume * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    step = 1.0e-6 * volume
    numerical_pressure = -(
        diamond_eos.helmholtz_free_energy(volume + step, temperature)
        - diamond_eos.helmholtz_free_energy(volume - step, temperature)
    ) / (2.0 * step * 1.0e4)

    assert diamond_eos.pressure(volume, temperature) == pytest.approx(
        numerical_pressure, rel=2.0e-9, abs=2.0e-8
    )


def test_volume_dependent_weight_derivative_is_part_of_ion_pressure(diamond_eos):
    volume = 5.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    temperature = 2000.0
    step = 1.0e-6 * volume
    numerical_pressure = -(
        diamond_eos.ion_helmholtz_free_energy(volume + step, temperature)
        - diamond_eos.ion_helmholtz_free_energy(volume - step, temperature)
    ) / (2.0 * step * 1.0e4)

    assert diamond_eos.ion_pressure(volume, temperature) == pytest.approx(
        numerical_pressure, rel=2.0e-9
    )


def test_vinet_cold_pressure_is_negative_cold_energy_derivative(diamond_eos):
    volume = 5.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    step = 1.0e-6 * volume
    numerical_pressure = -(
        diamond_eos.cold_energy(volume + step) - diamond_eos.cold_energy(volume - step)
    ) / (2.0 * step * 1.0e4)

    assert diamond_eos.rt_eos.pressure(volume) == pytest.approx(
        numerical_pressure, rel=1.0e-9
    )


def test_anharmonic_pressure_is_thermodynamically_consistent(diamond_eos):
    eos = diamond_eos.with_parameters(kappa=1.25)
    volume = 5.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    temperature = 2000.0
    step = 1.0e-6 * volume
    numerical_pressure = -(
        eos.anharmonic_helmholtz_free_energy(volume + step, temperature)
        - eos.anharmonic_helmholtz_free_energy(volume - step, temperature)
    ) / (2.0 * step * 1.0e4)

    assert eos.anharmonic_pressure(volume, temperature) == pytest.approx(
        numerical_pressure, rel=1.0e-9
    )
    assert diamond_eos.anharmonic_pressure(volume, temperature) == 0.0


@pytest.mark.parametrize(
    "atomic_volume,temperature,expected_pressure",
    [
        (5.7034, 300.0, 5.097381463860107),
        (5.4, 1000.0, 34.72410455615707),
        (5.0, 2000.0, 87.83106753791182),
        (4.654270411587497, 3000.0, 149.99999999999955),
    ],
)
def test_benedict_table_i_diamond_pressure_regression(
    diamond_eos, atomic_volume, temperature, expected_pressure
):
    """Pin states calculated directly from the published analytic model."""
    volume = atomic_volume * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR

    assert diamond_eos.pressure(volume, temperature) == pytest.approx(
        expected_pressure, rel=5.0e-11
    )


def test_pressure_volume_and_temperature_round_trips(diamond_eos):
    expected_volume = 4.654270411587497 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    expected_temperature = 3000.0
    pressure = diamond_eos.pressure(expected_volume, expected_temperature)

    assert diamond_eos.volume(pressure, expected_temperature) == pytest.approx(
        expected_volume, rel=1.0e-11
    )
    assert diamond_eos.temperature(pressure, expected_volume) == pytest.approx(
        expected_temperature, rel=1.0e-11
    )


@pytest.mark.parametrize("f_dac", [0.0, 0.25, 0.7])
def test_dac_volume_pair_inversion_uses_reference_relative_pressure(diamond_eos, f_dac):
    heated_volume = 5.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    expected_temperature = 2000.0
    reference_pressure = diamond_eos.pressure(heated_volume, diamond_eos.Tr)
    thermal_increment = (
        diamond_eos.pressure(heated_volume, expected_temperature) - reference_pressure
    )
    ambient_pressure = (
        diamond_eos.pressure(heated_volume, expected_temperature)
        - f_dac * thermal_increment
    )
    ambient_volume = diamond_eos.volume(ambient_pressure, diamond_eos.Tr)

    recovered = diamond_eos.temperature_from_volumes(
        ambient_volume,
        heated_volume,
        f_dac=f_dac,
    )

    assert recovered == pytest.approx(expected_temperature, rel=1.0e-11)
    assert diamond_eos.dac_thermal_pressure(
        heated_volume, recovered, f_dac
    ) == pytest.approx(f_dac * thermal_increment, rel=1.0e-12)
    assert diamond_eos.pressure(heated_volume, recovered) == pytest.approx(
        diamond_eos.pressure(ambient_volume, diamond_eos.Tr)
        + diamond_eos.dac_thermal_pressure(heated_volume, recovered, f_dac),
        rel=1.0e-12,
    )


def test_dac_volume_pair_inversion_broadcasts(diamond_eos):
    expected_temperatures = np.array([1000.0, 2000.0, 3000.0])
    heated_volumes = np.array([5.2, 5.0, 4.8]) * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    f_dac = 0.25
    increments = diamond_eos.pressure(
        heated_volumes, expected_temperatures
    ) - diamond_eos.pressure(heated_volumes, diamond_eos.Tr)
    ambient_pressures = (
        diamond_eos.pressure(heated_volumes, expected_temperatures) - f_dac * increments
    )
    ambient_volumes = diamond_eos.volume(ambient_pressures, diamond_eos.Tr)

    recovered = diamond_eos.temperature_from_volumes(
        ambient_volumes, heated_volumes, f_dac=f_dac
    )

    assert np.allclose(recovered, expected_temperatures, rtol=1.0e-11)


def test_dac_forward_volume_uses_reference_relative_pressure(diamond_eos):
    cold_pressure = 60.0
    temperature = 2000.0
    f_dac = 0.25

    confined_volume = diamond_eos.volume_with_dac_confinement(
        cold_pressure,
        temperature,
        f_dac=f_dac,
    )
    thermal_increment = diamond_eos.pressure(
        confined_volume, temperature
    ) - diamond_eos.pressure(confined_volume, diamond_eos.Tr)

    assert diamond_eos.thermal_pressure_increment(
        confined_volume, temperature
    ) == pytest.approx(thermal_increment, rel=1.0e-12)
    assert diamond_eos.pressure(confined_volume, temperature) == pytest.approx(
        cold_pressure + f_dac * thermal_increment,
        rel=1.0e-11,
    )


def test_dac_volume_pair_rejects_nonheated_state(diamond_eos):
    with pytest.raises(ValueError, match="below the reference temperature"):
        diamond_eos.temperature_from_volumes(
            diamond_eos.rt_eos.V0,
            0.99 * diamond_eos.rt_eos.V0,
            f_dac=0.25,
        )


def test_parameter_reconstruction_and_fixed_cold_curve_fitting(diamond_eos):
    parameters = diamond_eos.parameter_values()
    rebuilt = diamond_eos.with_parameters(theta_a0=1900.0)

    assert parameters["rt_eos.K0"] == 432.4
    assert rebuilt.theta_a0 == 1900.0
    assert rebuilt.rt_eos.K0 == diamond_eos.rt_eos.K0

    volumes = np.array([5.4, 5.2, 5.0, 4.8]) * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    temperatures = np.array([500.0, 1000.0, 1500.0, 2500.0])
    pressures = diamond_eos.pressure(volumes, temperatures)
    own_parameters = diamond_eos.parameter_values(include_reference=False)
    result = fit_thermal_eos(
        DoubleDebyeHelmholtz,
        diamond_eos.rt_eos,
        volumes,
        temperatures,
        pressures,
        initial={"theta_b0": 1750.0},
        fixed={
            name: value for name, value in own_parameters.items() if name != "theta_b0"
        },
        absolute_sigma=True,
    )

    assert result.success
    assert result.model.theta_b0 == pytest.approx(1887.8, rel=2.0e-9)


@pytest.mark.parametrize("temperature", [-1.0, np.nan, np.inf])
def test_invalid_temperature_is_rejected(diamond_eos, temperature):
    with pytest.raises(ValueError, match="Temperature"):
        diamond_eos.pressure(diamond_eos.rt_eos.V0, temperature)


def test_non_vinet_cold_curve_is_rejected():
    with pytest.raises(TypeError, match="Vinet cold curve"):
        DoubleDebyeHelmholtz(
            BM3(1.0, 100.0, 4.0),
            1.0,
            100.0,
            0.0,
            1.0,
            200.0,
            0.0,
            1.0,
            150.0,
            0.0,
            1.0,
        )


def test_correa_logarithmic_moment_weights_follow_equation_13(
    correa_diamond_eos,
):
    volume = 4.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    theta_a, theta_b, theta_0 = correa_diamond_eos.debye_temperatures(volume)
    weight_a, weight_b = correa_diamond_eos.double_debye_weights(volume)

    assert weight_a == pytest.approx(
        np.log(theta_b / theta_0) / np.log(theta_b / theta_a)
    )
    assert weight_a + weight_b == pytest.approx(1.0)
    assert (weight_a * np.log(theta_a) + weight_b * np.log(theta_b)) == pytest.approx(
        np.log(theta_0)
    )


def test_correa_degenerate_reference_weight_has_finite_limit(correa_diamond_eos):
    weight_a, weight_b = correa_diamond_eos.double_debye_weights(correa_diamond_eos.Vp)

    gamma_a = -0.316 * 5.571 + 0.913
    gamma_b = 0.168 * 5.571 + 0.429
    gamma_0 = 0.131 * 5.571 + 0.202
    assert weight_a == pytest.approx((gamma_b - gamma_0) / (gamma_b - gamma_a))
    assert weight_a + weight_b == pytest.approx(1.0)
    assert np.isfinite(correa_diamond_eos.pressure(correa_diamond_eos.Vp, 300.0))


@pytest.mark.parametrize(
    "atomic_volume,temperature",
    [(5.785, 0.0), (5.571, 300.0), (4.43, 5000.0), (3.21, 9000.0)],
)
def test_correa_pressure_is_negative_free_energy_volume_derivative(
    correa_diamond_eos, atomic_volume, temperature
):
    volume = atomic_volume * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    step = 1.0e-6 * volume
    numerical_pressure = -(
        (
            correa_diamond_eos.helmholtz_free_energy(volume + step, temperature)
            - correa_diamond_eos.phi0
        )
        - (
            correa_diamond_eos.helmholtz_free_energy(volume - step, temperature)
            - correa_diamond_eos.phi0
        )
    ) / (2.0 * step * 1.0e4)

    assert correa_diamond_eos.pressure(volume, temperature) == pytest.approx(
        numerical_pressure, rel=1.0e-8, abs=5.0e-8
    )


def test_correa_anharmonic_term_changes_caloric_but_not_pressure(
    correa_diamond_eos,
):
    volume = 4.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    temperature = 5000.0
    harmonic = correa_diamond_eos.with_parameters(anharmonic_a=0.0)

    assert correa_diamond_eos.pressure(volume, temperature) == pytest.approx(
        harmonic.pressure(volume, temperature)
    )
    assert correa_diamond_eos.anharmonic_helmholtz_free_energy(
        volume, temperature
    ) == pytest.approx(-R * 3.8e-5 * temperature**2)
    assert (
        correa_diamond_eos.molar_heat_capacity_v(volume, temperature)
        - harmonic.molar_heat_capacity_v(volume, temperature)
    ) == pytest.approx(2.0 * R * 3.8e-5 * temperature)


def test_published_correa_benedict_anharmonic_factor_of_two_is_explicit(
    correa_diamond_eos, diamond_eos
):
    """The records reproduce the two papers' inconsistent normalizations."""
    volume = 4.0 * ATOMIC_ANGSTROM3_TO_MOLAR_J_PER_BAR
    temperature = 5000.0

    correa_coefficient = -correa_diamond_eos.anharmonic_helmholtz_free_energy(
        volume, temperature
    ) / (R * temperature**2)
    benedict_coefficient = -diamond_eos.anharmonic_helmholtz_free_energy(
        volume, temperature
    ) / (R * temperature**2)

    assert correa_coefficient == pytest.approx(3.8e-5)
    assert benedict_coefficient == pytest.approx(3.79e-5 / 2.0)


def test_correa_parameter_reconstruction(correa_diamond_eos):
    rebuilt = correa_diamond_eos.with_parameters(theta_0_0=1900.0)

    assert rebuilt.theta_0_0 == 1900.0
    assert rebuilt.theta_a0 == correa_diamond_eos.theta_a0
    assert rebuilt.rt_eos.K0 == correa_diamond_eos.rt_eos.K0
