import csv
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import Avogadro, gas_constant
from scipy.optimize import least_squares

from peritheos import Material, load_eosmat, validate_eosmat_document

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos" / "data" / "materials"
DATASETS = ROOT / "peritheos" / "data" / "datasets"
DOI = "10.1002/2015JB012108"
PARAMETER_ORDER = ("V0", "K0", "K0_prime", "gamma0", "q")
_DEBYE_NODES, _DEBYE_WEIGHTS = np.polynomial.legendre.leggauss(64)

CASES = {
    "mg087fe013sio3_bridgmanite": {
        "record": "mg087fe013sio3_bridgmanite_wolf_2015_vinet_mgd_1",
        "formula": "Mg0.87Fe0.13SiO3",
        "parameters": (163.16, 243.8, 4.160, 1.400, 0.56),
        "errors": (0.19, 4.3, 0.110, 0.110, 0.37),
        "correlation": (
            (1.00, -0.85, 0.62, 0.03, 0.03),
            (-0.85, 1.00, -0.93, -0.32, -0.32),
            (0.62, -0.93, 1.00, 0.48, 0.50),
            (0.03, -0.32, 0.48, 1.00, 0.97),
            (0.03, -0.32, 0.50, 0.97, 1.00),
        ),
        "dataset": "mg087fe013sio3_bridgmanite_wolf_2015_table1_pvt",
        "resource": "mg087fe013sio3-bridgmanite-wolf-2015-table1-pvt.csv",
        "sha256": "dee6828c62a0682550d1a7a0b68c4d3bce6bc5c8b3f09b9c7417c36709eda070",
        "source_sha256": "6fa55c4e36e108301695b7a599f3ed6e104534363effe5b8cff3d78f5bad9958",
        "rows": 153,
        "groups": {1: 49, 2: 104},
        "range": (33.15, 132.70, 300.0, 2455.0, 122.42, 146.59),
        "first": (1, 33.15, 300.0, 146.59, 34.707),
        "last": (153, 132.70, 2062.0, 124.59, 25.010),
        "marker": "neon_unit_cell_volume_a3",
        "calibration_status": "partially_resolved",
        "recalculation_status": "reference_eos_not_bundled",
        "refit": (163.16078947, 243.83083910, 4.15944883, 1.40147771, 0.55884896),
        "published_rmse": 1.9496990046,
        "table5": (2387.0, 104.4872, 5.40, 604.5, 1.23, 0.9838),
    },
    "bridgmanite": {
        "record": "bridgmanite_wolf_2015_vinet_mgd_4",
        "formula": "MgSiO3",
        "parameters": (162.12, 262.3, 4.044, 1.675, 1.39),
        "errors": (0.13, 3.2, 0.075, 0.045, 0.16),
        "correlation": (
            (1.00, -0.96, 0.85, -0.12, -0.16),
            (-0.96, 1.00, -0.96, -0.07, -0.02),
            (0.85, -0.96, 1.00, 0.27, 0.24),
            (-0.12, -0.07, 0.27, 1.00, 0.95),
            (-0.16, -0.02, 0.24, 0.95, 1.00),
        ),
        "dataset": "bridgmanite_wolf_2015_table2_pvt",
        "resource": "bridgmanite-wolf-2015-table2-pvt.csv",
        "sha256": "9bf59dee8dcd228402e9bd640fec0c5bf926b3b567139ad11065708a50e3e015",
        "source_sha256": "30e76256dbf22d014425109cf6c4611e8d6012598ffd38be2babec682bd8e8c9",
        "rows": 42,
        "groups": {1: 27, 2: 6, 3: 9},
        "range": (27.96, 108.51, 300.0, 2430.0, 128.10, 150.26),
        "first": (1, 30.43, 1300.0, 150.260, 66.4150),
        "last": (42, 88.84, 300.0, 130.140, 55.5370),
        "marker": "mgo_unit_cell_volume_a3",
        "calibration_status": "resolved",
        "recalculation_status": "ready",
        "refit": (162.12798543, 262.21990748, 4.04745932, 1.67933146, 1.41547507),
        "published_rmse": 0.3027410708,
        "table5": (2429.0, 100.3870, 5.17, 603.7, 1.22, 0.9839),
    },
}


def _document(material_identifier):
    return load_eosmat(MATERIALS / f"{material_identifier}.eosmat")


def _record_document(document, identifier):
    return next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )


def _rows(case):
    with (DATASETS / case["resource"]).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _debye_function_3(argument):
    argument = np.asarray(argument, dtype=float)
    integration_points = 0.5 * argument[..., None] * (_DEBYE_NODES + 1.0)
    integral = (
        0.5
        * argument
        * np.sum(
            _DEBYE_WEIGHTS * integration_points**3 / np.expm1(integration_points),
            axis=-1,
        )
    )
    return 3.0 * integral / argument**3


def _source_pressure(parameters, volume, temperature):
    """Wolf equations 1-6, independently evaluated in published cell units."""
    V0, K0, K0_prime, gamma0, q = parameters
    volume = np.asarray(volume, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    x = (volume / V0) ** (1.0 / 3.0)
    cold_pressure = (
        3.0 * K0 * (1.0 - x) / x**2 * np.exp(1.5 * (K0_prime - 1.0) * (1.0 - x))
    )
    gamma = gamma0 * (volume / V0) ** q
    theta = 1000.0 * np.exp((gamma0 - gamma) / q)

    def energy(at_temperature):
        return (
            3.0
            * 5.0
            * gas_constant
            * at_temperature
            * _debye_function_3(theta / at_temperature)
        )

    molar_volume_j_per_bar = volume * Avogadro * 1.0e-25 / 4.0
    thermal_pressure = (
        gamma
        * (energy(temperature) - energy(np.full_like(temperature, 300.0)))
        / molar_volume_j_per_bar
        / 1.0e4
    )
    return cold_pressure + thermal_pressure


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_wolf_material_identity_model_and_provenance(material_identifier, case):
    document = _document(material_identifier)
    validate_eosmat_document(document)
    record = _record_document(document, case["record"])

    assert document["formula"] == case["formula"]
    assert document["space_group"] == "Pbnm"
    assert document["space_group_number"] == 62
    assert document["formula_units_per_cell"] == 4
    assert document["units"]["volume"] == "angstrom^3/conventional_unit_cell"
    assert record["reference"]["doi"] == DOI
    assert record["equation_kind"] == "thermal"
    assert record["eos"]["model"] == "vinet"
    assert tuple(record["eos"]["parameters"].values()) == case["parameters"][:3]
    assert record["thermal"]["model"] == "mie_gruneisen_debye"
    assert record["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert record["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 1000.0,
        "gamma0": case["parameters"][3],
        "q": case["parameters"][4],
        "n": 5.0,
    }
    assert record["fixed_parameters"] == []
    assert record["thermal"]["fixed_parameters"] == ["Tr", "theta0", "n"]
    assert record["parameter_error_confidence"] == 0.68
    assert record["scientific_validation"]["status"] == "primary_source_validated"
    assert record["fit_datasets"] == [case["dataset"]]
    assert record["pressure_calibration"]["status"] == case["calibration_status"]
    assert (
        record["pressure_calibration"]["recalculation"]["status"]
        == case["recalculation_status"]
    )

    if material_identifier.startswith("mg087"):
        a_site = [site for site in document["atom_sites"] if site["wyckoff"] == "4c"][
            :2
        ]
        assert [(site["element"], site["occupancy"]) for site in a_site] == [
            ("Mg", 0.87),
            ("Fe", 0.13),
        ]
        assert "proxies are not used by the EOS" in document["notes"]


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_wolf_complete_pvt_table_transcription(material_identifier, case):
    document = _document(material_identifier)
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == case["dataset"]
    )
    path = DATASETS / case["resource"]
    rows = _rows(case)

    assert hashlib.sha256(path.read_bytes()).hexdigest() == case["sha256"]
    assert dataset["resource"]["sha256"] == case["sha256"]
    assert dataset["source_file_sha256"] == case["source_sha256"]
    assert dataset["used_by_eos_records"] == [case["record"]]
    assert len(rows) == case["rows"]
    assert [int(row["source_order"]) for row in rows] == list(
        range(1, case["rows"] + 1)
    )
    assert Counter(int(row["measurement_group"]) for row in rows) == case["groups"]

    def numeric(name):
        return np.array([float(row[name]) for row in rows])

    expected_range = case["range"]
    assert (
        numeric("pressure_gpa").min(),
        numeric("pressure_gpa").max(),
    ) == pytest.approx(expected_range[:2])
    assert (
        numeric("temperature_k").min(),
        numeric("temperature_k").max(),
    ) == pytest.approx(expected_range[2:4])
    assert (
        numeric("bridgmanite_unit_cell_volume_a3").min(),
        numeric("bridgmanite_unit_cell_volume_a3").max(),
    ) == pytest.approx(expected_range[4:])

    for row, expected in ((rows[0], case["first"]), (rows[-1], case["last"])):
        actual = tuple(
            float(row[name])
            for name in (
                "source_order",
                "pressure_gpa",
                "temperature_k",
                "bridgmanite_unit_cell_volume_a3",
                case["marker"],
            )
        )
        assert actual == pytest.approx(expected)


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_wolf_table4_covariance_preserves_all_correlations(material_identifier, case):
    document = _document(material_identifier)
    source = _record_document(document, case["record"])
    executable = Material.from_eosmat(
        document, record_identifiers=[case["record"]]
    ).eos_records[0]

    assert source["parameter_covariance"]["parameter_order"] == [
        "rt_eos.V0",
        "rt_eos.K0",
        "rt_eos.K0_prime",
        "gamma0",
        "q",
    ]
    covariance = np.asarray(executable.parameter_covariance)
    standard_deviations = np.asarray(case["errors"], dtype=float)
    standard_deviations[0] *= executable.volume_scale
    recovered_correlation = covariance / np.outer(
        standard_deviations, standard_deviations
    )
    assert recovered_correlation == pytest.approx(
        np.asarray(case["correlation"]), abs=5.0e-12
    )


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_wolf_thermal_records_execute_invert_and_reproduce_table5(
    material_identifier, case
):
    record = Material.from_eosmat(
        _document(material_identifier), record_identifiers=[case["record"]]
    ).eos_records[0]

    assert record.pressure(record.reference_volume, 300.0) == pytest.approx(
        0.0, abs=1.0e-12
    )
    assert record.thermal_pressure_increment(140.0, 300.0) == pytest.approx(
        0.0, abs=1.0e-12
    )
    volume = record.volume(100.0, 2000.0, check_validity=True)
    assert record.pressure(volume, 2000.0, check_validity=True) == pytest.approx(
        100.0, rel=1.0e-11
    )

    temperature, molar_mass, density, bulk_modulus_s, gamma, cv_fraction = case[
        "table5"
    ]
    volume = record.volume(108.7, temperature)
    internal_volume = volume * record.volume_scale
    calculated_density = 4.0 * molar_mass / (Avogadro * volume * 1.0e-24)
    calculated_cv_fraction = record.eos.molar_heat_capacity_v(
        internal_volume, temperature
    ) / (3.0 * 5.0 * gas_constant)
    assert calculated_density == pytest.approx(density, abs=0.006)
    assert record.eos.adiabatic_bulk_modulus(
        internal_volume, temperature
    ) == pytest.approx(bulk_modulus_s, abs=0.12)
    assert record.eos.gruneisen_parameter(internal_volume) == pytest.approx(
        gamma, abs=0.006
    )
    assert calculated_cv_fraction == pytest.approx(cv_fraction, abs=5.0e-5)


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_wolf_rounded_tables_independently_refit_to_published_coefficients(
    material_identifier, case
):
    rows = _rows(case)
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    pressure_sigma = np.array(
        [float(row["pressure_total_uncertainty_gpa"]) for row in rows]
    )
    temperature = np.array([float(row["temperature_k"]) for row in rows])
    volume = np.array([float(row["bridgmanite_unit_cell_volume_a3"]) for row in rows])
    published = np.asarray(case["parameters"])

    def residual(parameters):
        data_residual = (
            _source_pressure(parameters, volume, temperature) - pressure
        ) / pressure_sigma
        prior_residual = (
            parameters[[0, 3, 4]]
            - np.array(
                [163.2 if material_identifier.startswith("mg087") else 162.5, 1.0, 1.0]
            )
        ) / np.array([0.2, 1.0, 1.0])
        return np.concatenate((data_residual, prior_residual))

    fit = least_squares(
        residual,
        published,
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
        max_nfev=2000,
    )
    published_rmse = np.sqrt(
        np.mean((_source_pressure(published, volume, temperature) - pressure) ** 2)
    )

    assert fit.success
    assert published_rmse == pytest.approx(case["published_rmse"], abs=2.0e-10)
    assert fit.x == pytest.approx(case["refit"], abs=6.0e-6)
    assert np.all(np.abs(fit.x - published) < np.asarray(case["errors"]))
    stored_refit = _record_document(_document(material_identifier), case["record"])[
        "scientific_validation"
    ]["independent_refit"]["refit_parameters"]
    assert tuple(stored_refit[name] for name in PARAMETER_ORDER) == pytest.approx(
        fit.x, abs=6.0e-6
    )
