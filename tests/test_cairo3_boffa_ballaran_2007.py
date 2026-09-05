import csv
import hashlib
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document
from peritheos.eosmat import validate_pressure_calibration_references

ROOT = Path(__file__).parents[1]
DOI = "10.2138/am.2007.2715"
CASES = {
    "cairo3_perovskite": {
        "record": "cairo3_perovskite_boffa_ballaran_2007_bm3_1",
        "dataset": "cairo3_perovskite_boffa_ballaran_2007_table1_compression",
        "resource": "cairo3-perovskite-boffa-ballaran-2007-table1-compression.csv",
        "sha256": "fc1b5f7c0458a08e64c689c6c631b29f6711938d63d12cb5d19aca600903d41a",
        "space_group": "Pbnm",
        "space_group_number": 62,
        "lattice_volume": 229.526,
        "parameters": {"V0": 229.463, "K0": 198.0, "K0_prime": 1.2},
        "errors": {"V0": 0.008, "K0": 3.0, "K0_prime": 0.8},
        "high_pressure_volume": 220.801,
        "high_pressure_calculated": 7.779735,
        "pressure_rmse": 0.0277249632,
        "maximum_residual": 0.0524508168,
        "figure_2_slope": (-825.0, 244.0),
        "refit_parameters": (229.46343467, 198.37055335, 1.23310114),
    },
    "cairo3_post_perovskite": {
        "record": "cairo3_post_perovskite_boffa_ballaran_2007_bm3_1",
        "dataset": "cairo3_post_perovskite_boffa_ballaran_2007_table1_compression",
        "resource": "cairo3-post-perovskite-boffa-ballaran-2007-table1-compression.csv",
        "sha256": "c6c182484b36c270c2782a124562a43f1c3fbefd2635c1706f9f1ffae554b7d9",
        "space_group": "Cmcm",
        "space_group_number": 63,
        "lattice_volume": 226.75,
        "parameters": {"V0": 226.38, "K0": 181.0, "K0_prime": 2.3},
        "errors": {"V0": 0.01, "K0": 3.0, "K0_prime": 0.8},
        "high_pressure_volume": 217.28,
        "high_pressure_calculated": 7.776392,
        "pressure_rmse": 0.0239824255,
        "maximum_residual": 0.04679018,
        "figure_2_slope": (-433.0, 214.0),
        "refit_parameters": (226.37736314, 181.46697214, 2.27631675),
    },
}


def _rows(case):
    path = ROOT / "peritheos" / "data" / "datasets" / case["resource"]
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _source_bm3(beta, volume):
    """Standard BM3 in the Eulerian convention defined by source Figure 2."""
    v0, k0, k0_prime = beta
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_cairo3_polymorph_identity_and_published_record(material_identifier, case):
    document = get_material_document(material_identifier)
    assert document["formula"] == "CaIrO3"
    assert document["space_group"] == case["space_group"]
    assert document["space_group_number"] == case["space_group_number"]
    assert document["formula_units_per_cell"] == 4

    lattice = document["lattice"]
    calculated_lattice_volume = lattice["a"] * lattice["b"] * lattice["c"]
    assert calculated_lattice_volume == pytest.approx(case["lattice_volume"], abs=0.002)

    cell_contents = Counter()
    for site in document["atom_sites"]:
        cell_contents[site["element"]] += site["site_multiplicity"] * site["occupancy"]
    assert cell_contents == {"Ca": 4.0, "Ir": 4.0, "O": 12.0}

    assert len(document["eos_records"]) == 1
    record = document["eos_records"][0]
    assert record["identifier"] == case["record"]
    assert record["reference"]["doi"].lower() == DOI
    assert record["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": case["parameters"],
    }
    assert record["parameter_errors"] == case["errors"]
    assert record["parameter_error_confidence"] is None
    assert record["parameter_covariance"] is None
    assert record["fixed_parameters"] == []
    assert record["experimental_pressure_range_gpa"] == [0.0001, 7.79]
    assert (
        record["pressure_calibration"]["methods"][0]["reference_calibration_record"]
        == "ruby_mao_1986"
    )
    assert record["pressure_calibration"]["recalculation"]["status"] == "ready"
    assert record["scientific_validation"]["status"] == "primary_source_validated"


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_cairo3_table1_transcription_and_checksum(material_identifier, case):
    document = get_material_document(material_identifier)
    dataset = document["datasets"][0]
    assert dataset["identifier"] == case["dataset"]
    assert dataset["used_by_eos_records"] == [case["record"]]

    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == case["sha256"]
    assert dataset["resource"]["sha256"] == case["sha256"]

    rows = _rows(case)
    assert len(rows) == 11
    assert len(rows[0]) == 12
    assert rows[0]["pressure_gpa"] == "0.00010"
    assert float(rows[0]["pressure_standard_deviation_gpa"]) == 0.00001
    assert float(rows[9]["pressure_gpa"]) == 7.79
    assert float(rows[9]["volume_a3"]) == case["high_pressure_volume"]
    assert rows[10]["loading_path"] == "2"
    assert rows[10]["fit_included"] == "0"
    assert sum(int(row["fit_included"]) for row in rows) == 10

    for row in rows:
        lattice_volume = (
            float(row["a_angstrom"])
            * float(row["b_angstrom"])
            * float(row["c_angstrom"])
        )
        volume_sigma = float(row["volume_standard_deviation_a3"])
        assert lattice_volume == pytest.approx(
            float(row["volume_a3"]), abs=max(0.006, 2.0 * volume_sigma)
        )


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_cairo3_bm3_reproduces_table1_and_figure2(material_identifier, case):
    document = get_material_document(material_identifier)
    material = Material.from_eosmat(document, record_identifiers=[case["record"]])
    record = material.get_eos_record(case["record"])
    rows = [row for row in _rows(case) if row["fit_included"] == "1"]
    volumes = np.array([float(row["volume_a3"]) for row in rows])
    source_pressures = np.array([float(row["pressure_gpa"]) for row in rows])

    calculated_pressures = np.asarray(record.pressure(volumes), dtype=float)
    residuals = calculated_pressures - source_pressures
    assert math.sqrt(float(np.mean(residuals**2))) == pytest.approx(
        case["pressure_rmse"], abs=5.0e-10
    )
    assert float(np.max(np.abs(residuals))) == pytest.approx(
        case["maximum_residual"], abs=5.0e-9
    )
    assert record.pressure(case["high_pressure_volume"]) == pytest.approx(
        case["high_pressure_calculated"], abs=5.0e-7
    )
    assert record.pressure(record.reference_volume) == pytest.approx(0.0, abs=1.0e-12)
    assert record.eos.bulk_modulus(record.reference_volume) == pytest.approx(
        case["parameters"]["K0"], rel=1.0e-12
    )

    source_slope, source_sigma = case["figure_2_slope"]
    implied_slope = (
        1.5 * case["parameters"]["K0"] * (case["parameters"]["K0_prime"] - 4.0)
    )
    assert implied_slope == pytest.approx(source_slope, abs=source_sigma)

    pressure_grid = np.linspace(0.0001, 7.79, 9)
    volume_grid = record.volume(pressure_grid)
    assert record.pressure(volume_grid) == pytest.approx(pressure_grid, rel=1.0e-10)


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_cairo3_independent_refit_has_parameter_parity(material_identifier, case):
    rows = [row for row in _rows(case) if row["fit_included"] == "1"]
    pressures = np.array([float(row["pressure_gpa"]) for row in rows])
    pressure_sigmas = np.array(
        [float(row["pressure_standard_deviation_gpa"]) for row in rows]
    )
    volumes = np.array([float(row["volume_a3"]) for row in rows])
    volume_sigmas = np.array(
        [float(row["volume_standard_deviation_a3"]) for row in rows]
    )

    def weighted_residuals(beta):
        step = 1.0e-5
        pressure_derivative = (
            _source_bm3(beta, volumes + step) - _source_bm3(beta, volumes - step)
        ) / (2.0 * step)
        effective_sigmas = np.sqrt(
            pressure_sigmas**2 + (pressure_derivative * volume_sigmas) ** 2
        )
        return (_source_bm3(beta, volumes) - pressures) / effective_sigmas

    result = least_squares(
        weighted_residuals,
        x0=tuple(case["parameters"].values()),
    )

    assert result.x == pytest.approx(case["refit_parameters"], abs=2.0e-4)
    published = np.array(tuple(case["parameters"].values()))
    published_errors = np.array(tuple(case["errors"].values()))
    assert np.all(np.abs(result.x - published) < published_errors)

    validation = get_material_document(material_identifier)["eos_records"][0][
        "scientific_validation"
    ]
    stored_refit = tuple(validation["independent_refit"]["parameters"].values())
    assert stored_refit == pytest.approx(result.x, abs=2.0e-4)


def test_cairo3_pressure_calibration_links_resolve_globally():
    validate_pressure_calibration_references()
