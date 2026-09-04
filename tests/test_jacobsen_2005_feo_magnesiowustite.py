"""Primary-source checks for Jacobsen et al. (2005) FeO and (Mg,Fe)O."""

import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document
from peritheos.eosmat import validate_pressure_calibration_references

ROOT = Path(__file__).resolve().parents[1]
DOI = "10.1107/S0909049505022326"
DATASET_ID = "jacobsen_2005_table1_compression"
DATASET_SHA256 = "afe24343cd64bd3b17e1a0d708dbd79a5511d6d43755af0cd082c0cde533820f"
CASES = {
    "mg073fe027o": {
        "record": "mg073fe027o_jacobsen_2005_bm3_1",
        "parameters": (77.30, 154.0, 4.0),
        "errors": (0.09, 3.0, 0.1),
        "fixed": [],
        "flag": "mg073fe027o_fit_included",
        "volume": "mg073fe027o_derived_volume_a3",
        "volume_sigma": "mg073fe027o_derived_volume_sd_a3",
        "lattice": "mg073fe027o_a_angstrom",
        "lattice_sigma": "mg073fe027o_a_sd_angstrom",
        "observations": 16,
        "range": (0.0, 51.1),
        "refit": (77.30368820, 154.85488377, 3.97082563),
        "published_rmse": 0.1939114407,
        "published_max": 0.4161484547,
    },
    "feo": {
        "record": "fe093o_b1_jacobsen_2005_bm3_1",
        "parameters": (79.41, 146.0, 4.0),
        "errors": (0.04, 2.0, None),
        "fixed": ["K0_prime"],
        "flag": "fe093o_b1_fit_included",
        "volume": "fe093o_b1_derived_volume_a3",
        "volume_sigma": "fe093o_b1_derived_volume_sd_a3",
        "lattice": "fe093o_b1_a_angstrom",
        "lattice_sigma": "fe093o_b1_a_sd_angstrom",
        "observations": 7,
        "range": (0.0, 22.8),
        "refit": (79.40783326, 145.80527709),
        "published_rmse": 0.1199987974,
        "published_max": 0.1874971095,
    },
    "fe093o_rhombohedral": {
        "record": "fe093o_rhombohedral_jacobsen_2005_bm3_1",
        "parameters": (59.72, 134.0, 4.0),
        "errors": (0.22, 4.0, None),
        "fixed": ["K0_prime"],
        "flag": "fe093o_rhombohedral_fit_included",
        "volume": "fe093o_rhombohedral_volume_a3",
        "volume_sigma": "fe093o_rhombohedral_volume_sd_a3",
        "observations": 9,
        "range": (27.7, 51.1),
        "refit": (59.75515644, 133.40916064),
        "published_rmse": 0.4415688144,
        "published_max": 0.7741052573,
    },
}


def _document_record_dataset(material_identifier):
    document = get_material_document(material_identifier)
    record = next(
        item
        for item in document["eos_records"]
        if item["identifier"] == CASES[material_identifier]["record"]
    )
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    return document, record, dataset


def _rows():
    path = (
        ROOT
        / "peritheos"
        / "data"
        / "datasets"
        / "jacobsen-2005-table1-compression.csv"
    )
    payload = path.read_bytes()
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    return path, payload, rows


def _bm3(parameters, volume, fixed_k0_prime=None):
    v0, k0 = parameters[:2]
    k0_prime = fixed_k0_prime if fixed_k0_prime is not None else parameters[2]
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5
        * k0
        * (eta**7 - eta**5)
        * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def _selected_observations(case, rows):
    selected = [row for row in rows if row[case["flag"]] == "1"]
    pressure = np.asarray([float(row["pressure_gpa"]) for row in selected])
    if "lattice" in case:
        lattice = np.asarray([float(row[case["lattice"]]) for row in selected])
        lattice_sigma = np.asarray(
            [float(row[case["lattice_sigma"]]) for row in selected]
        )
        volume = lattice**3
        volume_sigma = 3.0 * lattice**2 * lattice_sigma
    else:
        volume = np.asarray([float(row[case["volume"]]) for row in selected])
        volume_sigma = np.asarray(
            [float(row[case["volume_sigma"]]) for row in selected]
        )
    return pressure, volume, volume_sigma


def _effective_variance_refit(case, rows):
    pressure, volume, volume_sigma = _selected_observations(case, rows)
    fixed_k0_prime = 4.0 if case["fixed"] else None

    def residuals(parameters):
        step = 1.0e-5
        derivative = (
            _bm3(parameters, volume + step, fixed_k0_prime)
            - _bm3(parameters, volume - step, fixed_k0_prime)
        ) / (2.0 * step)
        # Table 1 supplies a range, not row-wise pressure errors. Its midpoint is
        # a deterministic audit choice and is explicitly recorded as such.
        effective_sigma = np.sqrt(0.075**2 + (derivative * volume_sigma) ** 2)
        return (
            _bm3(parameters, volume, fixed_k0_prime) - pressure
        ) / effective_sigma

    initial = case["parameters"][:2] if case["fixed"] else case["parameters"]
    result = least_squares(residuals, initial, x_scale="jac")
    assert result.success
    return result.x


def test_jacobsen_2005_preserves_composition_phase_and_volume_bases():
    mg, mg_record, _ = _document_record_dataset("mg073fe027o")
    rhombo, rhombo_record, _ = _document_record_dataset("fe093o_rhombohedral")

    assert mg["formula"] == "Mg0.73Fe0.27O"
    assert mg["space_group"] == "Fm-3m"
    assert mg["formula_units_per_cell"] == 4
    assert mg["lattice"]["a"] ** 3 == pytest.approx(77.303332, abs=5.0e-7)
    assert [site["occupancy"] for site in mg["atom_sites"][:2]] == [0.73, 0.27]

    assert rhombo["formula"] == "Fe0.93O"
    assert rhombo["formula_units_per_cell"] == 3
    assert "equivalent hexagonal cell" in rhombo["phase"]
    assert "space_group" not in rhombo
    molar_volume = (
        rhombo_record["eos"]["parameters"]["V0"] * 0.602214076 / 3.0
    )
    assert molar_volume == pytest.approx(11.99, abs=0.005)
    assert mg_record["temperature_ref"] == rhombo_record["temperature_ref"] == 300.0


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_jacobsen_2005_registers_only_source_supported_bm3_variants(
    material_identifier, case
):
    document, record, dataset = _document_record_dataset(material_identifier)
    parameters = record["eos"]["parameters"]

    assert record["reference"]["doi"] == DOI
    assert record["eos"]["type"] == "BM3"
    assert tuple(parameters.values()) == case["parameters"]
    assert tuple(record["parameter_errors"].values()) == case["errors"]
    assert record["fixed_parameters"] == case["fixed"]
    assert record["parameter_covariance"] is None
    assert tuple(record["experimental_pressure_range_gpa"]) == case["range"]
    assert record["fit_datasets"] == [DATASET_ID]
    assert dataset["used_by_eos_records"] == [case["record"]]
    assert record["scientific_validation"]["primary_data_check"]["status"] == (
        "bundled"
    )

    executable = Material.from_eosmat(
        document, record_identifiers=[case["record"]]
    ).eos_records[0]
    assert executable.identifier == case["record"]


def test_jacobsen_2005_table1_is_complete_checksummed_and_unambiguous():
    path, payload, rows = _rows()
    assert hashlib.sha256(payload).hexdigest() == DATASET_SHA256
    assert len(rows) == 16
    assert rows[0]["pressure_gpa"] == "0.00"
    assert rows[0]["pressure_precision_min_gpa"] == ""
    assert rows[-1]["pressure_gpa"] == "51.1"
    assert rows[-1]["fe093o_rhombohedral_volume_a3"] == "47.00"
    assert sum(row["mg073fe027o_fit_included"] == "1" for row in rows) == 16
    assert sum(row["fe093o_b1_fit_included"] == "1" for row in rows) == 7
    assert sum(row["fe093o_rhombohedral_fit_included"] == "1" for row in rows) == 9

    for material_identifier in CASES:
        _, _, dataset = _document_record_dataset(material_identifier)
        assert dataset["resource"]["sha256"] == DATASET_SHA256
        assert [column["name"] for column in dataset["columns"]] == list(rows[0])

    for row in rows:
        lattice = float(row["mg073fe027o_a_angstrom"])
        lattice_sigma = float(row["mg073fe027o_a_sd_angstrom"])
        assert float(row["mg073fe027o_derived_volume_a3"]) == pytest.approx(
            lattice**3, abs=5.0e-7
        )
        assert float(row["mg073fe027o_derived_volume_sd_a3"]) == pytest.approx(
            3.0 * lattice**2 * lattice_sigma, abs=5.0e-7
        )
    assert path.name == "jacobsen-2005-table1-compression.csv"


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_jacobsen_2005_published_curves_and_independent_refits_have_parity(
    material_identifier, case
):
    _, record, _ = _document_record_dataset(material_identifier)
    _, _, rows = _rows()
    pressure, volume, _ = _selected_observations(case, rows)
    residual = _bm3(case["parameters"], volume) - pressure

    assert math.sqrt(float(np.mean(residual**2))) == pytest.approx(
        case["published_rmse"], abs=5.0e-10
    )
    assert float(np.max(np.abs(residual))) == pytest.approx(
        case["published_max"], abs=5.0e-10
    )

    refit = _effective_variance_refit(case, rows)
    assert refit == pytest.approx(case["refit"], abs=2.0e-6)
    published = np.asarray(case["parameters"][: len(refit)])
    errors = np.asarray(case["errors"][: len(refit)], dtype=float)
    assert np.all(np.abs(refit - published) < errors)

    stored = record["scientific_validation"]["numerical_reproduction"][
        "independent_refit"
    ]
    assert tuple(stored["parameters"].values())[: len(refit)] == pytest.approx(
        refit, abs=2.0e-6
    )
    assert "midpoint" in stored["objective"]


def test_jacobsen_2005_rhombohedral_alternatives_and_calibration_are_explicit():
    _, record, _ = _document_record_dataset("fe093o_rhombohedral")
    parameterizations = record["scientific_validation"][
        "reported_parameterizations"
    ]
    assert [item["role"] for item in parameterizations] == [
        "finite_pressure_fit",
        "finite_pressure_fit_back_calculated_to_zero_pressure",
        "preferred_ambient_reference_cross_check",
    ]
    assert parameterizations[0]["reference_pressure_gpa"] == 28.0
    assert parameterizations[-1]["fixed_parameters"] == ["K0_prime"]

    for material_identifier in CASES:
        _, source_record, _ = _document_record_dataset(material_identifier)
        calibration = source_record["pressure_calibration"]
        assert calibration["status"] == "resolved"
        assert calibration["recalculation"]["status"] == "ready"
        assert calibration["methods"][0]["reference_calibration_record"] == (
            "ruby_mao_1986"
        )
    validate_pressure_calibration_references()
