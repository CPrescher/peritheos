import csv
import hashlib
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document, validate_eosmat_document
from peritheos.eosmat import validate_pressure_calibration_references

ROOT = Path(__file__).resolve().parents[1]
MATERIAL_ID = "mg095al010si095o3_bridgmanite"
BM2_ID = "mg095al010si095o3_bridgmanite_daniel_2004_bm2_1"
BM3_ID = "mg095al010si095o3_bridgmanite_daniel_2004_bm3_2"
DATASET_ID = "mg095al010si095o3_bridgmanite_daniel_2004_table1_compression"
DOI = "10.1029/2004GL020213"
DATASET_SHA256 = "5c259a2b56bb8d40dd4e05cf8fa98c5101a009fb830771fd5b6f0f1897f7cfcb"


def _document_records_dataset_rows():
    document = get_material_document(MATERIAL_ID)
    validate_eosmat_document(document)
    source_records = {
        record["identifier"]: record for record in document["eos_records"]
    }
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    payload = path.read_bytes()
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    material = Material.from_eosmat(document, record_identifiers=[BM2_ID, BM3_ID])
    records = {record.identifier: record for record in material.eos_records}
    return document, source_records, dataset, payload, rows, records


def _bm3_pressure(volumes, parameters):
    v0, k0, k0_prime = parameters
    eta = (v0 / volumes) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def _weighted_refit(rows, *, bm2):
    pressures = np.asarray([float(row["pressure_gpa"]) for row in rows])
    pressure_errors = np.asarray(
        [float(row["pressure_uncertainty_gpa"]) for row in rows]
    )
    volumes = np.asarray([float(row["volume_a3_conventional_cell"]) for row in rows])

    def residuals(parameters):
        if bm2:
            parameters = np.asarray([parameters[0], parameters[1], 4.0])
        return (_bm3_pressure(volumes, parameters) - pressures) / pressure_errors

    initial = [163.234, 251.5] if bm2 else [163.23, 251.5, 4.0]
    fit = least_squares(
        residuals,
        x0=initial,
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    )
    dof = len(rows) - len(fit.x)
    reduced_chi_square = float(np.sum(fit.fun**2) / dof)
    covariance = np.linalg.inv(fit.jac.T @ fit.jac) * reduced_chi_square
    standard_errors = np.sqrt(np.diag(covariance))
    correlation = covariance / np.sqrt(
        np.outer(np.diag(covariance), np.diag(covariance))
    )
    return fit.x, standard_errors, correlation, reduced_chi_square


def test_daniel_material_identity_structure_and_cell_basis():
    document, source_records, _, _, rows, _ = _document_records_dataset_rows()

    assert document["formula"] == "Mg0.95Al0.10Si0.95O3"
    assert document["phase"] == (
        "orthorhombic aluminous silicate perovskite, Al.05-Pv specimen"
    )
    assert document["space_group"] == "Pbnm"
    assert document["space_group_number"] == 62
    assert document["formula_units_per_cell"] == 4
    assert set(source_records) == {BM2_ID, BM3_ID}

    lattice = document["lattice"]
    assert lattice["a"] * lattice["b"] * lattice["c"] == pytest.approx(
        float(rows[0]["volume_a3_conventional_cell"]), abs=0.002
    )
    assert rows[0]["pressure_gpa"] == "0.37"
    assert "topology proxy" in document["cell_contents"]
    assert "do not encode" in document["cell_contents"]
    assert document["source"]["same_batch_composition_reference"]["doi"] == (
        "10.1029/2004GL019918"
    )
    assert document["source"]["topology_proxy_reference"]["doi"] == (
        "10.1007/BF00308114"
    )

    proxy_contents = Counter()
    for site in document["atom_sites"]:
        multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
        proxy_contents[site["element"]] += multiplicity * site["occupancy"]
    assert proxy_contents == {"Mg": 4.0, "Si": 4.0, "O": 12.0}


def test_daniel_published_bm2_and_bm3_parameters_are_distinct_and_executable():
    _, sources, _, _, _, records = _document_records_dataset_rows()
    bm2_source = sources[BM2_ID]
    bm3_source = sources[BM3_ID]

    assert bm2_source["reference"]["doi"] == DOI
    assert bm2_source["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 163.234, "K0": 251.5},
    }
    assert bm2_source["parameter_errors"] == {"V0": 0.008, "K0": 1.3}
    assert bm2_source["fixed_parameters"] == []
    assert bm2_source["parameter_covariance"] == {
        "parameter_order": ["V0", "K0"],
        "matrix": [[0.000064, -0.00468], [-0.00468, 1.69]],
    }

    assert bm3_source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 163.23, "K0": 251.5, "K0_prime": 4.0},
    }
    assert bm3_source["parameter_errors"] == {
        "V0": 0.01,
        "K0": 2.7,
        "K0_prime": 0.3,
    }
    assert bm3_source["fixed_parameters"] == []
    assert bm3_source["parameter_covariance"] is None

    for record_id, v0 in ((BM2_ID, 163.234), (BM3_ID, 163.23)):
        record = records[record_id]
        assert record.pressure(v0, 300.0) == pytest.approx(0.0, abs=1.0e-12)
        assert record.eos.bulk_modulus(v0) == pytest.approx(251.5)
        for pressure in (0.37, 20.0, 39.97):
            volume = record.volume(pressure, 300.0, check_validity=True)
            assert record.pressure(volume, 300.0, check_validity=True) == pytest.approx(
                pressure, rel=1.0e-11
            )
        with pytest.raises(ValueError, match="outside the published calibration"):
            record.volume(40.0, 300.0, check_validity=True)


def test_daniel_official_table_is_complete_checksummed_and_lossless():
    _, sources, dataset, payload, rows, _ = _document_records_dataset_rows()

    assert hashlib.sha256(payload).hexdigest() == DATASET_SHA256
    assert dataset["resource"]["sha256"] == DATASET_SHA256
    assert dataset["used_by_eos_records"] == [BM2_ID, BM3_ID]
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == 42
    assert [int(row["source_order"]) for row in rows] == list(range(1, 43))
    assert all(row["used_in_published_fit"] == "1" for row in rows)
    assert Counter(row["pressure_calibration"] for row in rows) == {
        "mao_1986_ruby": 33,
        "hemley_1989_neon": 9,
    }
    assert sum(int(row["laser_annealed"]) for row in rows) == 9
    assert all(
        (float(row["pressure_gpa"]) > 26.0) == bool(int(row["laser_annealed"]))
        for row in rows
    )

    assert rows[0] == {
        "source_order": "1",
        "pressure_gpa": "0.37",
        "pressure_uncertainty_gpa": "0.01",
        "lattice_a_angstrom": "4.7757",
        "lattice_a_uncertainty_angstrom": "0.0001",
        "lattice_b_angstrom": "4.9349",
        "lattice_b_uncertainty_angstrom": "0.0001",
        "lattice_c_angstrom": "6.9166",
        "lattice_c_uncertainty_angstrom": "0.0001",
        "volume_a3_conventional_cell": "163.008",
        "volume_uncertainty_a3": "0.008",
        "pressure_calibration": "mao_1986_ruby",
        "laser_annealed": "0",
        "used_in_published_fit": "1",
    }
    assert rows[-1]["pressure_gpa"] == "38.86"
    assert rows[-1]["volume_a3_conventional_cell"] == "144.517"
    assert max(float(row["pressure_gpa"]) for row in rows) == pytest.approx(39.97)

    for row in rows:
        lattice_volume = math.prod(
            float(row[name])
            for name in (
                "lattice_a_angstrom",
                "lattice_b_angstrom",
                "lattice_c_angstrom",
            )
        )
        assert lattice_volume == pytest.approx(
            float(row["volume_a3_conventional_cell"]), abs=0.0035
        )

    for source in sources.values():
        check = source["scientific_validation"]["primary_source_check"]
        assert check["official_auxiliary_readme_sha256"] == (
            "816470aa2c394a62ace9406bf2677abeac0ffe0e42ead7009c3421dde8d0c10d"
        )
        assert check["official_auxiliary_table_sha256"] == (
            "492cfcfc883014fba90295a6cfd8c3fad5e3a0f240abee46fa93a86728e739d4"
        )


@pytest.mark.parametrize(
    ("record_id", "expected_rmse", "expected_max"),
    (
        (BM2_ID, 0.5986176369, 1.8190284258),
        (BM3_ID, 0.5990070581, 1.8092065727),
    ),
)
def test_daniel_published_curves_reproduce_complete_table(
    record_id, expected_rmse, expected_max
):
    _, sources, _, _, rows, records = _document_records_dataset_rows()
    pressures = np.asarray([float(row["pressure_gpa"]) for row in rows])
    volumes = np.asarray([float(row["volume_a3_conventional_cell"]) for row in rows])
    calculated = np.asarray(records[record_id].pressure(volumes), dtype=float)
    residuals = calculated - pressures

    assert math.sqrt(float(np.mean(residuals**2))) == pytest.approx(
        expected_rmse, abs=5.0e-10
    )
    assert float(np.max(np.abs(residuals))) == pytest.approx(expected_max, abs=5.0e-10)
    stored = sources[record_id]["scientific_validation"]["numerical_reproduction"]
    assert stored["published_curve_pressure_rmse_gpa"] == pytest.approx(expected_rmse)
    assert stored["published_curve_max_abs_pressure_residual_gpa"] == pytest.approx(
        expected_max
    )


@pytest.mark.parametrize(
    ("record_id", "bm2", "parameters", "errors", "reduced_chi_square"),
    (
        (
            BM2_ID,
            True,
            [163.23191791, 251.98496449],
            [0.00847082, 1.30420339],
            2.1776535235,
        ),
        (
            BM3_ID,
            False,
            [163.23220727, 251.83536437, 4.01758767],
            [0.00967836, 2.66099490, 0.27172321],
            2.2332506076,
        ),
    ),
)
def test_daniel_pressure_weighted_refits_have_coefficient_and_error_parity(
    record_id, bm2, parameters, errors, reduced_chi_square
):
    _, sources, _, _, rows, _ = _document_records_dataset_rows()
    fitted, fitted_errors, correlation, fitted_reduced_chi_square = _weighted_refit(
        rows, bm2=bm2
    )

    assert fitted == pytest.approx(parameters, abs=5.0e-7)
    assert fitted_errors == pytest.approx(errors, abs=5.0e-7)
    assert fitted_reduced_chi_square == pytest.approx(reduced_chi_square, abs=5.0e-9)

    reproduction = sources[record_id]["scientific_validation"][
        "numerical_reproduction"
    ]["independent_refit"]
    assert reproduction["conclusion"] == "parity"
    assert list(reproduction["parameters"].values()) == pytest.approx(fitted)
    assert list(reproduction["scaled_standard_errors"].values()) == pytest.approx(
        fitted_errors
    )

    published = sources[record_id]["eos"]["parameters"]
    published_errors = sources[record_id]["parameter_errors"]
    names = ["V0", "K0"] if bm2 else ["V0", "K0", "K0_prime"]
    assert all(
        abs(fitted[index] - published[name]) <= published_errors[name]
        for index, name in enumerate(names)
    )
    if bm2:
        assert correlation[0, 1] == pytest.approx(-0.46622237, abs=5.0e-7)
        assert correlation[0, 1] == pytest.approx(-0.45, abs=0.02)


def test_daniel_calibration_scope_and_recalculation_limits_are_explicit():
    _, sources, _, _, _, _ = _document_records_dataset_rows()

    for source in sources.values():
        calibration = source["pressure_calibration"]
        assert calibration["status"] == "partially_resolved"
        assert calibration["methods"][0]["reference_calibration_record"] == (
            "ruby_mao_1986"
        )
        assert calibration["methods"][0]["reference"]["doi"] == (
            "10.1029/JB091iB05p04673"
        )
        assert calibration["methods"][1]["reference"]["doi"] == (
            "10.1103/PhysRevB.39.11820"
        )
        assert calibration["recalculation"]["status"] == (
            "missing_calibrant_observations"
        )
        assert "neon lattice parameters" in calibration["recalculation"]["notes"]
        assert source["scientific_validation"]["status"] == ("primary_source_validated")
        assert source["fit_datasets"] == [DATASET_ID]

    validate_pressure_calibration_references()
