import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, load_eosmat

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos" / "data" / "materials"
DATASETS = ROOT / "peritheos" / "data" / "datasets"

CASES = {
    "calcium_titanate_perovskite": {
        "record": "calcium_titanate_perovskite_ross_1999_bm3_1",
        "dataset": "calcium-titanate-perovskite-ross-1999-table1-compression.csv",
        "sha256": "a1f9acdcd0bbf98839dc9d0b85ae398fd6869a7b2bfa06df513edb52448e1680",
        "rows": 11,
        "parameters": (223.764, 170.9, 6.6),
        "refit": (223.76464009, 170.85581942, 6.57739326),
        "reduced_chi_square": 1.3808862174,
        "published_rmse": 0.0189459986,
        "published_max_residual": 0.0296793188,
        "last_row": (9.700, 213.207, 0.020),
        "cell_atoms": {"Ca": 4, "Ti": 4, "O": 12},
    },
    "calcium_germanate_perovskite": {
        "record": "calcium_germanate_perovskite_ross_1999_bm3_1",
        "dataset": "calcium-germanate-perovskite-ross-1999-table1-compression.csv",
        "sha256": "fc407eca7e0a841c1a2ceb5b995fb8df22596eba81238221924cf30992520ace",
        "rows": 9,
        "parameters": (206.490, 194.0, 6.1),
        "refit": (206.48942690, 194.03279818, 6.09858913),
        "reduced_chi_square": 1.1028596230,
        "published_rmse": 0.0191108763,
        "published_max_residual": 0.0304503891,
        "last_row": (8.553, 198.577, 0.015),
        "cell_atoms": {"Ca": 4, "Ge": 4, "O": 12},
    },
}


def _document(material_id):
    return load_eosmat(MATERIALS / f"{material_id}.eosmat")


def _rows(filename):
    with (DATASETS / filename).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _bm3_pressure(volume, V0, K0, K0_prime):
    compression = (V0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5
        * K0
        * (compression**7 - compression**5)
        * (1.0 + 0.75 * (K0_prime - 4.0) * (compression**2 - 1.0))
    )


@pytest.mark.parametrize("material_id", CASES)
def test_ross_1999_records_load_execute_and_invert(material_id):
    case = CASES[material_id]
    document = _document(material_id)
    record_document = document["eos_records"][0]
    record = Material.from_eosmat(
        document, record_identifiers=[case["record"]]
    ).eos_records[0]

    assert record_document["reference"]["doi"] == "10.2138/am-1999-0309"
    assert record_document["eos"]["parameters"] == dict(
        zip(("V0", "K0", "K0_prime"), case["parameters"], strict=True)
    )
    assert document["space_group"] == "Pbnm"
    assert document["space_group_number"] == 62
    assert document["formula_units_per_cell"] == 4
    assert [site["wyckoff"] for site in document["atom_sites"]] == [
        "4c",
        "4b",
        "4c",
        "8d",
    ]
    cell_atoms = {}
    for site, multiplicity in zip(document["atom_sites"], (4, 4, 4, 8), strict=True):
        cell_atoms[site["element"]] = cell_atoms.get(site["element"], 0) + int(
            multiplicity * site["occupancy"]
        )
    assert cell_atoms == case["cell_atoms"]
    assert record_document["volume_basis"]["formula_units"] == 4.0
    assert record_document["fixed_parameters"] == []
    assert record_document["parameter_covariance"] is None
    calibration = record_document["pressure_calibration"]
    assert calibration["status"] == "resolved"
    assert calibration["methods"][0]["reference_eos_record"] == (
        "alpha_quartz_angel_1997_bm3_1"
    )
    assert calibration["recalculation"]["status"] == ("missing_calibrant_observations")

    pressure, volume, _ = case["last_row"]
    calculated = record.pressure(volume, check_validity=True)
    assert abs(calculated - pressure) <= 0.031
    assert record.volume(calculated, check_validity=True) == pytest.approx(
        volume, abs=1.0e-8
    )


@pytest.mark.parametrize("material_id", CASES)
def test_ross_1999_table1_transcription_and_checksum(material_id):
    case = CASES[material_id]
    path = DATASETS / case["dataset"]
    rows = _rows(case["dataset"])

    assert hashlib.sha256(path.read_bytes()).hexdigest() == case["sha256"]
    assert len(rows) == case["rows"]
    assert list(rows[0]) == [
        "source_row",
        "pressure_gpa",
        "pressure_sd_gpa",
        "lattice_a_angstrom",
        "lattice_a_sd_angstrom",
        "lattice_b_angstrom",
        "lattice_b_sd_angstrom",
        "lattice_c_angstrom",
        "lattice_c_sd_angstrom",
        "unit_cell_volume_a3",
        "unit_cell_volume_sd_a3",
    ]
    assert rows[0]["pressure_sd_gpa"] == ""
    document = _document(material_id)
    lattice_volume = np.prod([document["lattice"][axis] for axis in ("a", "b", "c")])
    assert lattice_volume == pytest.approx(
        float(rows[0]["unit_cell_volume_a3"]), abs=0.002
    )
    assert tuple(
        float(rows[-1][name])
        for name in (
            "pressure_gpa",
            "unit_cell_volume_a3",
            "unit_cell_volume_sd_a3",
        )
    ) == pytest.approx(case["last_row"])


@pytest.mark.parametrize("material_id", CASES)
def test_ross_1999_published_curve_and_independent_effective_variance_refit(
    material_id,
):
    case = CASES[material_id]
    rows = _rows(case["dataset"])
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["unit_cell_volume_a3"]) for row in rows])
    volume_sigma = np.array([float(row["unit_cell_volume_sd_a3"]) for row in rows])
    pressure_sigma = np.array([float(row["pressure_sd_gpa"] or 0.0) for row in rows])

    published_pressure = _bm3_pressure(volume, *case["parameters"])
    published_residual = published_pressure - pressure
    assert np.sqrt(np.mean(published_residual**2)) == pytest.approx(
        case["published_rmse"], abs=1.0e-10
    )
    assert np.max(np.abs(published_residual)) == pytest.approx(
        case["published_max_residual"], abs=1.0e-10
    )

    # Ross and Angel minimized pressure residuals weighted by the sample-volume
    # e.s.d. transformed to pressure and the quartz-derived pressure e.s.d.
    # The near-ambient rows have no printed pressure e.s.d.; zero here means
    # that only their reported volume e.s.d. contributes to the weight.
    def effective_variance_residual(parameters):
        V0, K0, K0_prime = parameters
        compression = (V0 / volume) ** (1.0 / 3.0)
        correction = 1.0 + 0.75 * (K0_prime - 4.0) * (compression**2 - 1.0)
        model_pressure = _bm3_pressure(volume, *parameters)
        derivative_compression = (
            1.5
            * K0
            * (
                (7.0 * compression**6 - 5.0 * compression**4) * correction
                + (compression**7 - compression**5)
                * 1.5
                * (K0_prime - 4.0)
                * compression
            )
        )
        derivative = derivative_compression * -compression / (3.0 * volume)
        sigma = np.sqrt(pressure_sigma**2 + (derivative * volume_sigma) ** 2)
        return (model_pressure - pressure) / sigma

    fit = least_squares(
        effective_variance_residual,
        np.asarray(case["parameters"]),
        bounds=(0.0, np.inf),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
        max_nfev=2000,
    )
    reduced_chi_square = np.sum(fit.fun**2) / (len(rows) - 3)

    assert fit.success
    assert fit.x == pytest.approx(case["refit"], abs=2.0e-6)
    assert reduced_chi_square == pytest.approx(case["reduced_chi_square"], abs=2.0e-8)
    assert round(reduced_chi_square, 1) == (
        1.4 if material_id == "calcium_titanate_perovskite" else 1.1
    )
