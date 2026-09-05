"""Regression checks against numerical results in the source papers.

These tests intentionally distinguish the equations printed in the papers from
the calculation path distributed in the Sokolova et al. (2016) workbooks.
"""

import csv
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import R
from scipy.optimize import brentq

from peritheos.eos.rt import BM3, Holzapfel, Vinet
from peritheos.eos.thermal import (
    DorogokupetsOganov2007,
    MieGruneisenEinstein,
    MultiOscillatorGruneisenThermalEOS,
)
from peritheos.eos.thermal.multi_oscillator import I_gamV
from peritheos.materials import (
    AG_SOKOLOVA_2013,
    AL_SOKOLOVA_2013,
    AU_SOKOLOVA_2013,
    CU_SOKOLOVA_2013,
    DIAMOND_SOKOLOVA_2013,
    MGO_SOKOLOVA_2013,
    MO_SOKOLOVA_2013,
    NB_SOKOLOVA_2013,
    PT_SOKOLOVA_2013,
    TA_SOKOLOVA_2013,
    W_SOKOLOVA_2013,
)


def _dorogokupets_2007_eos(parameters):
    return DorogokupetsOganov2007(
        Vinet(parameters["V0"], parameters["K0"], parameters["K0_prime"]),
        Tr=298.15,
        theta_B1=parameters["theta_B1"],
        d_B1=parameters["d_B1"],
        m_B1=parameters["m_B1"],
        theta_B2=parameters.get("theta_B2", 1.0),
        d_B2=parameters.get("d_B2", 1.0),
        m_B2=parameters.get("m_B2", 1.0e-12),
        theta_E1=parameters["theta_E1"],
        m_E1=parameters["m_E1"],
        theta_E2=parameters["theta_E2"],
        m_E2=parameters["m_E2"],
        gamma0=parameters["gamma0"],
        gamma_inf=parameters["gamma_inf"],
        beta=parameters["beta"],
        anharmonic_a=parameters["anharmonic_a"],
        anharmonic_m=parameters["anharmonic_m"],
        electronic_e=parameters.get("electronic_e", 0.0),
        electronic_g=parameters.get("electronic_g", 0.0),
        defect_H=parameters.get("defect_H", 1.0e9),
        defect_S=parameters.get("defect_S", 0.0),
        n=parameters["n"],
    )


@pytest.mark.parametrize(
    ("parameters", "ratio", "temperature", "expected", "tolerance"),
    [
        (
            dict(
                V0=1.0272,
                K0=99.65,
                K0_prime=6.11,
                theta_B1=130.6,
                d_B1=8.572,
                m_B1=0.121,
                theta_B2=103.6,
                d_B2=5.326,
                m_B2=0.449,
                theta_E1=111.9,
                m_E1=0.766,
                theta_E2=189.12,
                m_E2=1.664,
                gamma0=2.376,
                gamma_inf=1.481,
                beta=2.507,
                anharmonic_a=6.70,
                anharmonic_m=3.44,
                electronic_e=25.9,
                electronic_g=0.666,
                defect_H=15239,
                defect_S=0.732,
                n=1,
            ),
            0.6,
            2500.0,
            233.648,
            0.003,
        ),
        (
            dict(
                V0=0.9999,
                K0=72.67,
                K0_prime=4.62,
                theta_B1=245.8,
                d_B1=5.575,
                m_B1=0.987,
                theta_E1=240.2,
                m_E1=1.000,
                theta_E2=356.2,
                m_E2=1.013,
                gamma0=2.144,
                gamma_inf=1.017,
                beta=3.942,
                anharmonic_a=5.14,
                anharmonic_m=3.44,
                electronic_e=54.1,
                electronic_g=1.8,
                defect_H=8679,
                defect_S=0.998,
                n=1,
            ),
            0.5,
            2500.0,
            230.556,
            0.17,
        ),
        (
            dict(
                V0=1.0215,
                K0=166.70,
                K0_prime=6.00,
                theta_B1=95.7,
                d_B1=8.290,
                m_B1=0.681,
                theta_B2=106.4,
                d_B2=3.239,
                m_B2=0.417,
                theta_E1=170.6,
                m_E1=1.063,
                theta_E2=105.2,
                m_E2=0.839,
                gamma0=2.965,
                gamma_inf=1.142,
                beta=3.030,
                anharmonic_a=25.33,
                anharmonic_m=3.79,
                electronic_e=18.92,
                electronic_g=0.66,
                defect_H=11690,
                defect_S=1.067,
                n=1,
            ),
            0.65,
            2500.0,
            256.35,
            0.01,
        ),
        (
            dict(
                V0=0.7113,
                K0=133.41,
                K0_prime=5.37,
                theta_B1=123.7,
                d_B1=3.776,
                m_B1=0.115,
                theta_B2=175.4,
                d_B2=10.372,
                m_B2=0.711,
                theta_E1=187.4,
                m_E1=0.756,
                theta_E2=286.9,
                m_E2=1.418,
                gamma0=1.974,
                gamma_inf=1.554,
                beta=4.647,
                anharmonic_a=3.50,
                anharmonic_m=3.46,
                electronic_e=27.698,
                electronic_g=0.666,
                defect_H=11687,
                defect_S=1.407,
                n=1,
            ),
            0.6,
            2500.0,
            265.581,
            0.23,
        ),
        (
            dict(
                V0=0.9091,
                K0=276.07,
                K0_prime=5.30,
                theta_B1=95.2,
                d_B1=8.199,
                m_B1=0.329,
                theta_B2=148.4,
                d_B2=4.005,
                m_B2=0.383,
                theta_E1=214.6,
                m_E1=1.211,
                theta_E2=140.8,
                m_E2=1.077,
                gamma0=2.802,
                gamma_inf=1.538,
                beta=5.550,
                anharmonic_a=160.9,
                anharmonic_m=4.06,
                electronic_e=260.0,
                electronic_g=2.4,
                defect_H=32572,
                defect_S=0.631,
                n=1,
            ),
            0.7,
            3000.0,
            262.523,
            0.02,
        ),
        (
            dict(
                V0=1.0851,
                K0=191.39,
                K0_prime=3.81,
                theta_B1=72.6,
                d_B1=5.536,
                m_B1=0.117,
                theta_B2=101.8,
                d_B2=24.513,
                m_B2=0.396,
                theta_E1=144.0,
                m_E1=1.118,
                theta_E2=214.9,
                m_E2=1.369,
                gamma0=1.714,
                gamma_inf=1.241,
                beta=6.825,
                anharmonic_a=61.9,
                anharmonic_m=4.00,
                electronic_e=167.0,
                electronic_g=1.3,
                defect_H=36278,
                defect_S=4.910,
                n=1,
            ),
            0.6,
            3000.0,
            258.537,
            0.16,
        ),
        (
            dict(
                V0=0.9545,
                K0=306.00,
                K0_prime=4.17,
                theta_B1=182.8,
                d_B1=13.270,
                m_B1=0.513,
                theta_B2=172.5,
                d_B2=3.305,
                m_B2=0.174,
                theta_E1=287.6,
                m_E1=1.166,
                theta_E2=213.8,
                m_E2=1.145,
                gamma0=1.553,
                gamma_inf=0.694,
                beta=3.698,
                anharmonic_a=-39.3,
                anharmonic_m=2.67,
                electronic_e=40.4,
                electronic_g=0.2,
                defect_H=14714,
                defect_S=0.672,
                n=1,
            ),
            0.7,
            3000.0,
            232.448,
            0.061,
        ),
        (
            dict(
                V0=1.1248,
                K0=160.31,
                K0_prime=4.18,
                theta_B1=447.3,
                d_B1=11.248,
                m_B1=1.429,
                theta_B2=384.0,
                d_B2=3.593,
                m_B2=0.276,
                theta_E1=703.8,
                m_E1=2.570,
                theta_E2=466.0,
                m_E2=1.725,
                gamma0=1.522,
                gamma_inf=1.111,
                beta=4.509,
                anharmonic_a=13.56,
                anharmonic_m=5.23,
                n=2,
            ),
            0.6,
            3000.0,
            243.935,
            0.145,
        ),
        (
            dict(
                V0=0.3417,
                K0=443.16,
                K0_prime=3.777,
                theta_B1=1202.1,
                d_B1=9.604,
                m_B1=1.163,
                theta_B2=1135.1,
                d_B2=3.380,
                m_B2=0.218,
                theta_E1=1687.2,
                m_E1=1.396,
                theta_E2=1033.7,
                m_E2=0.223,
                gamma0=0.820,
                gamma_inf=0.615,
                beta=10.121,
                anharmonic_a=-23.85,
                anharmonic_m=1.22,
                n=1,
            ),
            0.7,
            3000.0,
            316.383,
            0.015,
        ),
    ],
)
def test_dorogokupets_oganov_2007_tables_cover_every_material(
    parameters, ratio, temperature, expected, tolerance
):
    eos = _dorogokupets_2007_eos(parameters)

    assert eos.pressure(ratio * eos.rt_eos.V0, temperature) == pytest.approx(
        expected, abs=tolerance
    )


def test_dorogokupets_2010_fit_1_reproduces_published_mgo_pressures():
    eos = MieGruneisenEinstein(
        BM3(1.1248, 160.2, 3.99),
        Tr=298.15,
        theta0=599.0,
        gamma0=1.524,
        q=1.65,
        n=2.0,
    )

    for ratio, temperature, expected in (
        (1.0, 3000.0, 17.69),
        (0.85, 3000.0, 51.69),
        (0.64, 3663.0, 190.72),
    ):
        assert eos.pressure(ratio * eos.rt_eos.V0, temperature) == pytest.approx(
            expected, abs=0.004
        )


@pytest.mark.parametrize(
    ("record", "ratio", "temperature", "expected", "tolerance"),
    [
        (DIAMOND_SOKOLOVA_2013, 0.70, 3500.0, 326.865, 0.004),
        (AL_SOKOLOVA_2013, 0.50, 2000.0, 233.635, 0.003),
        (CU_SOKOLOVA_2013, 0.50, 3000.0, 556.037, 0.003),
        (NB_SOKOLOVA_2013, 0.50, 3500.0, 431.788, 0.006),
        (MO_SOKOLOVA_2013, 0.60, 3500.0, 396.307, 0.16),
        (AG_SOKOLOVA_2013, 0.50, 3000.0, 513.628, 0.003),
        (TA_SOKOLOVA_2013, 0.50, 3500.0, 496.098, 0.008),
        (W_SOKOLOVA_2013, 0.60, 3500.0, 456.855, 0.003),
        (PT_SOKOLOVA_2013, 0.60, 3000.0, 550.358, 0.004),
        (AU_SOKOLOVA_2013, 0.60, 3000.0, 380.613, 0.004),
    ],
)
def test_dorogokupets_2012_appendix_tables_cover_every_material(
    record, ratio, temperature, expected, tolerance
):
    eos = record.eos

    assert eos.pressure(ratio * eos.rt_eos.V0, temperature) == pytest.approx(
        expected, abs=tolerance
    )


def _paper_pressure(model, ratio, temperature):
    """Literal Sokolova et al. (2013) equations (6), (9), and (14)."""
    rt_eos = model.rt_eos
    volume = ratio * rt_eos.V0
    reference_pressure = rt_eos.pressure(volume)
    bulk_modulus = rt_eos.bulk_modulus(volume)
    bulk_modulus_derivative = rt_eos.bulk_modulus_derivative(volume)
    generalized_t = model.t - model.beta * np.cbrt(ratio)
    gamma = (
        -3.0 * bulk_modulus
        + 2.0 * reference_pressure * generalized_t
        + 9.0 * bulk_modulus * bulk_modulus_derivative
        - 6.0 * generalized_t * bulk_modulus
    ) / (6.0 * (3.0 * bulk_modulus - 2.0 * reference_pressure * generalized_t))
    gamma += model.delta
    volume_theta_factor = math.exp(
        I_gamV(ratio, model.delta, model.t, rt_eos, model.beta)
    )
    anharmonicity = model.a_0 * 1.0e-6 * ratio**model.m

    def absolute_thermal_pressure(current_temperature):
        pressure_factor = gamma - 0.5 * model.m * anharmonicity * current_temperature
        pressure_bar = 0.0
        for multiplicity, reference_theta in (
            (model.mE1, model.QE1o),
            (model.mE2, model.QE2o),
        ):
            theta = (
                reference_theta
                * volume_theta_factor
                * math.exp(0.5 * anharmonicity * current_temperature)
            )
            pressure_bar += (
                multiplicity
                * R
                * theta
                / math.expm1(theta / current_temperature)
                * pressure_factor
                / volume
            )
        pressure_bar += (
            1.5
            * model.n
            * R
            * model.e_0
            * 1.0e-6
            * ratio**model.g
            * model.g
            * current_temperature**2
            / volume
        )
        return pressure_bar / 10000.0

    return (
        reference_pressure
        + absolute_thermal_pressure(temperature)
        - absolute_thermal_pressure(model.Tr)
    )


def _corrected_2013_eos(material):
    if material == "Mo":
        return MultiOscillatorGruneisenThermalEOS(
            Holzapfel(0.9369, 260.5, 4.05, 1.0, 42.0),
            Tr=298.15,
            QE1o=356.0,
            mE1=1.5,
            QE2o=218.0,
            mE2=1.5,
            delta=-0.755,
            t=-0.735,
            a_0=0.0,
            m=0.0,
            e_0=123.9,
            g=3.5,
            n=1.0,
        )
    return MultiOscillatorGruneisenThermalEOS(
        Holzapfel(1.0215, 167.0, 5.75, 1.0, 79.0),
        Tr=298.15,
        QE1o=176.0,
        mE1=1.5,
        QE2o=84.0,
        mE2=1.5,
        delta=0.045,
        t=-0.463,
        a_0=0.0,
        m=0.0,
        e_0=0.0,
        g=0.0,
        n=1.0,
    )


@pytest.mark.parametrize(
    ("material", "ratio", "temperature", "expected"),
    [
        ("Mo", 1.0, 3500.0, 18.436),
        ("Mo", 0.6, 3500.0, 383.282),
        ("Au", 1.0, 3000.0, 19.151),
        ("Au", 0.5, 3000.0, 788.541),
    ],
)
def test_sokolova_2013_corrected_mo_and_au_tables(
    material, ratio, temperature, expected
):
    eos = _corrected_2013_eos(material)

    assert eos.pressure(ratio * eos.rt_eos.V0, temperature) == pytest.approx(
        expected, abs=0.004
    )


def test_sokolova_2013_corrected_mgo_parameters_do_not_reproduce_table_6():
    eos = MultiOscillatorGruneisenThermalEOS(
        Holzapfel(1.1248, 160.3, 4.25, 2.0, 10.34),
        Tr=298.15,
        QE1o=748.0,
        mE1=3.0,
        QE2o=401.0,
        mE2=3.0,
        delta=-0.25,
        t=0.583,
        a_0=17.4,
        m=5.5,
        e_0=0.0,
        g=0.0,
        n=2.0,
    )

    calculated = _paper_pressure(eos, 1.0, 3500.0)
    assert calculated == pytest.approx(18.4998333829, abs=1.0e-9)
    assert 19.100 - calculated == pytest.approx(0.6001666171, abs=1.0e-9)


def test_2016_workbook_and_printed_anharmonic_equation_diverge_for_mgo():
    eos = MGO_SOKOLOVA_2013.eos
    volume = eos.rt_eos.V0

    workbook_pressure = eos.pressure(volume, 3000.0)
    printed_equation_pressure = _paper_pressure(eos, 1.0, 3000.0)

    assert workbook_pressure == pytest.approx(16.2788956245, abs=1.0e-9)
    assert printed_equation_pressure == pytest.approx(19.5843501770, abs=1.0e-9)
    assert printed_equation_pressure - workbook_pressure == pytest.approx(
        3.3054545525, abs=1.0e-9
    )


def _digitized_rows(filename):
    path = Path(__file__).parents[1] / "docs" / "data" / filename
    with path.open(newline="") as handle:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def _dorogokupets_2010_fit_2_pressure(ratio, temperature):
    rt_eos = BM3(1.1248, 160.2, 3.99)
    volume = ratio * rt_eos.V0
    gamma_inf = 1.325
    gamma0 = 1.524
    exponent = 11.8
    gamma = gamma_inf + (gamma0 - gamma_inf) * ratio**exponent
    theta = (
        599.0
        * ratio ** (-gamma_inf)
        * math.exp((gamma0 - gamma_inf) / exponent * (1.0 - ratio**exponent))
    )

    def thermal_energy(current_temperature):
        return 6.0 * R * theta / math.expm1(theta / current_temperature)

    thermal_pressure = (
        gamma
        / volume
        * (thermal_energy(temperature) - thermal_energy(298.15))
        / 10000.0
    )
    return rt_eos.pressure(volume) + thermal_pressure


def test_dorogokupets_2010_figure_1_quantifies_both_model_mismatches():
    rows = _digitized_rows("dorogokupets-2010-figure1-digitized.csv")
    fit_1 = MieGruneisenEinstein(
        BM3(1.1248, 160.2, 3.99),
        Tr=298.15,
        theta0=599.0,
        gamma0=1.524,
        q=1.65,
        n=2.0,
    )
    fit_1_residuals = []
    fit_2_residuals = []

    for row in rows:
        pressure = row["pressure_gpa"]
        observed = row["speziale_fig6_curve_volume_ratio"]
        fit_1_ratio = fit_1.calculate_volume(pressure, 3663.0) / fit_1.rt_eos.V0
        fit_2_ratio = brentq(
            lambda ratio: _dorogokupets_2010_fit_2_pressure(ratio, 3663.0) - pressure,
            0.5,
            0.9,
        )
        assert fit_1_ratio == pytest.approx(row["fit1_volume_ratio"], abs=1.0e-9)
        assert fit_2_ratio == pytest.approx(row["fit2_volume_ratio"], abs=1.0e-9)
        fit_1_residuals.append(fit_1_ratio - observed)
        fit_2_residuals.append(fit_2_ratio - observed)

    fit_1_rmse = float(np.sqrt(np.mean(np.square(fit_1_residuals))))
    fit_2_rmse = float(np.sqrt(np.mean(np.square(fit_2_residuals))))
    assert fit_1_rmse == pytest.approx(0.0149405, abs=1.0e-7)
    assert fit_2_rmse == pytest.approx(0.00363754, abs=1.0e-8)
    assert fit_2_rmse < fit_1_rmse / 4.0

    # Figure 1 contains calculated curves, not experimental points. The paper
    # explicitly reports that neither fit reaches the graphically transcribed
    # Speziale endpoint near 209 GPa at x=0.64 and 3663 K.
    ratio = 0.64
    fit_1_pressure = fit_1.pressure(ratio * fit_1.rt_eos.V0, 3663.0)
    fit_2_pressure = _dorogokupets_2010_fit_2_pressure(ratio, 3663.0)
    assert fit_1_pressure == pytest.approx(190.72, abs=0.004)
    assert fit_2_pressure == pytest.approx(203.34, abs=0.004)
    assert 209.0 - fit_2_pressure == pytest.approx(5.6633, abs=0.01)


def _literal_zero_pressure_bulk_modulus(model, temperature):
    ratio = brentq(
        lambda current_ratio: _paper_pressure(model, current_ratio, temperature),
        1.0,
        1.2,
    )
    step = ratio * 1.0e-5
    pressure_derivative = (
        _paper_pressure(model, ratio + step, temperature)
        - _paper_pressure(model, ratio - step, temperature)
    ) / (2.0 * step)
    return -ratio * pressure_derivative


def test_sokolova_2016_figure_3_digitization_prefers_workbook_form():
    rows = _digitized_rows("sokolova-2016-figure3-kt-digitized.csv")
    eos = MGO_SOKOLOVA_2013.eos
    workbook_residuals = []
    printed_residuals = []

    for row in rows:
        temperature = row["temperature_k"]
        observed = row["isothermal_bulk_modulus_gpa"]
        volume = eos.calculate_volume(0.0, temperature)
        workbook_value = eos.bulk_modulus(volume, temperature)
        printed_value = _literal_zero_pressure_bulk_modulus(eos, temperature)
        assert workbook_value == pytest.approx(row["workbook_model_gpa"], abs=0.005)
        assert printed_value == pytest.approx(
            row["printed_equation_model_gpa"], abs=0.005
        )
        workbook_residuals.append(workbook_value - observed)
        printed_residuals.append(printed_value - observed)

    workbook_rmse = float(np.sqrt(np.mean(np.square(workbook_residuals))))
    printed_rmse = float(np.sqrt(np.mean(np.square(printed_residuals))))
    assert workbook_rmse == pytest.approx(0.4665953, abs=5.0e-6)
    assert printed_rmse == pytest.approx(9.7740234, abs=1.0e-4)
    assert workbook_rmse < 1.0
    assert printed_rmse > 20.0 * workbook_rmse


def test_2016_workbook_and_printed_equation_agree_when_anharmonicity_is_inactive():
    records = (
        DIAMOND_SOKOLOVA_2013,
        AL_SOKOLOVA_2013,
        CU_SOKOLOVA_2013,
        AG_SOKOLOVA_2013,
        AU_SOKOLOVA_2013,
        PT_SOKOLOVA_2013,
        NB_SOKOLOVA_2013,
        TA_SOKOLOVA_2013,
        MO_SOKOLOVA_2013,
        W_SOKOLOVA_2013,
    )
    for record in records:
        eos = record.eos
        for ratio in (1.0, 0.8, 0.6):
            volume = ratio * eos.rt_eos.V0
            assert eos.pressure(volume, 3000.0) == pytest.approx(
                _paper_pressure(eos, ratio, 3000.0), abs=5.0e-9
            )
