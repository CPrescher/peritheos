"""Primary-source reproduction for Jacobsen et al. (2002) (Mg,Fe)O."""

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
DOI = "10.1029/2001JB000490"

CASES = {
    "ferropericlase_fe27": {
        "record": "ferropericlase_fe27_jacobsen_2002_bm3_1",
        "dataset": "ferropericlase_fe27_jacobsen_2002_table5_compression",
        "resource": "ferropericlase-fe27-jacobsen-2002-table5-compression.csv",
        "sha256": "b14a733b90ad310ceebe566653cc9c80b4f586edcda3b106321c1ad6a14bc16e",
        "rows": 10,
        "composition": (0.725, 0.268, 0.007),
        "parameters": (76.336, 158.4, 5.49),
        "errors": (0.002, 0.4, 0.11),
        "refit": (76.33633739, 158.41974587, 5.48425764),
        "published_rmse": 0.0057383539,
        "published_max_residual": 0.0090052559,
        "quartz_rmse": 0.0056029636,
        "quartz_max_residual": 0.0081263206,
        "last_row": (9.292, 4.17061, 0.95033),
    },
    "magnesiowustite_fe56": {
        "record": "magnesiowustite_fe56_jacobsen_2002_bm3_1",
        "dataset": "magnesiowustite_fe56_jacobsen_2002_table4_compression",
        "resource": "magnesiowustite-fe56-jacobsen-2002-table4-compression.csv",
        "sha256": "0597a3737f66275eb1e86361510360dc99641c64bd1559a60594bf3a2322cd61",
        "rows": 11,
        "composition": (0.423, 0.541, 0.036),
        "parameters": (77.453, 155.8, 5.53),
        "errors": (0.004, 0.9, 0.23),
        "refit": (77.45691695, 155.19063404, 5.67628583),
        "published_rmse": 0.0128727844,
        "published_max_residual": 0.0193735635,
        "quartz_rmse": 0.0014032461,
        "quartz_max_residual": 0.0043130544,
        "last_row": (8.923, 4.19247, 0.95142),
    },
    "magnesiowustite_fe75": {
        "record": "magnesiowustite_fe75_jacobsen_2002_bm3_1",
        "dataset": "magnesiowustite_fe75_jacobsen_2002_table5_compression",
        "resource": "magnesiowustite-fe75-jacobsen-2002-table5-compression.csv",
        "sha256": "61eef1bb96b389d97c67867a314c01ed8c5b453da80a4d40f86856241fbae42a",
        "rows": 7,
        "composition": (0.240, 0.717, 0.043),
        "parameters": (78.082, 151.3, 5.55),
        "errors": (0.003, 0.6, 0.19),
        "refit": (78.08178497, 151.33110627, 5.54968606),
        "published_rmse": 0.0050573088,
        "published_max_residual": 0.0083438849,
        "quartz_rmse": 0.0049481673,
        "quartz_max_residual": 0.0066362842,
        "last_row": (7.179, 4.21444, 0.95864),
    },
}


def _document(material_identifier):
    return load_eosmat(MATERIALS / f"{material_identifier}.eosmat")


def _rows(case):
    with (DATASETS / case["resource"]).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _bm3(parameters, volume):
    v0, k0, k0_prime = parameters
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_jacobsen_2002_material_identity_and_published_record(
    material_identifier, case
):
    document = _document(material_identifier)
    record = document["eos_records"][0]
    magnesium, iron, vacancy = case["composition"]

    assert document["space_group"] == "Fm-3m"
    assert document["space_group_number"] == 225
    assert document["formula_units_per_cell"] == 4
    assert magnesium + iron + vacancy == pytest.approx(1.0)
    assert document["atom_sites"][0]["occupancy"] == magnesium
    assert document["atom_sites"][1]["occupancy"] == iron
    assert document["lattice"]["a"] ** 3 == pytest.approx(
        float(_rows(case)[0]["lattice_a_angstrom"]) ** 3
    )

    assert record["identifier"] == case["record"]
    assert record["reference"]["doi"] == DOI
    assert tuple(record["eos"]["parameters"].values()) == case["parameters"]
    assert tuple(record["parameter_errors"].values()) == case["errors"]
    assert record["fixed_parameters"] == []
    assert record["parameter_covariance"] is None
    assert record["volume_basis"]["formula_units"] == 4.0
    assert record["pressure_calibration"]["status"] == "resolved"
    assert (
        record["pressure_calibration"]["methods"][0]["reference_eos_record"]
        == "alpha_quartz_angel_1997_bm3_1"
    )
    assert record["scientific_validation"]["status"] == "primary_source_validated"


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_jacobsen_2002_tables_are_complete_lossless_and_checksummed(
    material_identifier, case
):
    document = _document(material_identifier)
    dataset = document["datasets"][0]
    path = DATASETS / case["resource"]
    rows = _rows(case)

    assert dataset["identifier"] == case["dataset"]
    assert dataset["used_by_eos_records"] == [case["record"]]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == case["sha256"]
    assert dataset["resource"]["sha256"] == case["sha256"]
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == case["rows"]
    assert rows[0]["pressure_gpa"] == "0.0001"
    assert rows[0]["pressure_sd_gpa"] == ""
    assert rows[0]["volume_ratio_sd"] == ""

    pressure, lattice_a, volume_ratio = case["last_row"]
    assert (
        float(rows[-1]["pressure_gpa"]),
        float(rows[-1]["lattice_a_angstrom"]),
        float(rows[-1]["volume_ratio"]),
    ) == pytest.approx((pressure, lattice_a, volume_ratio))

    ambient_volume = float(rows[0]["lattice_a_angstrom"]) ** 3
    for row in rows:
        lattice_volume = float(row["lattice_a_angstrom"]) ** 3
        ratio_sigma = float(row["volume_ratio_sd"] or 6.0e-6)
        assert lattice_volume / ambient_volume == pytest.approx(
            float(row["volume_ratio"]), abs=ratio_sigma
        )


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_jacobsen_2002_published_curves_execute_and_reproduce_tables(
    material_identifier, case
):
    document = _document(material_identifier)
    rows = _rows(case)
    volume = np.array([float(row["lattice_a_angstrom"]) ** 3 for row in rows])
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    calculated = _bm3(case["parameters"], volume)
    residual = calculated - pressure

    assert np.sqrt(np.mean(residual**2)) == pytest.approx(
        case["published_rmse"], abs=1.0e-10
    )
    assert np.max(np.abs(residual)) == pytest.approx(
        case["published_max_residual"], abs=1.0e-10
    )

    loaded = Material.from_eosmat(document, record_identifiers=[case["record"]])
    record = loaded.get_eos_record(case["record"])
    assert np.asarray(record.pressure(volume), dtype=float) == pytest.approx(calculated)
    assert record.volume(record.pressure(volume)) == pytest.approx(volume, rel=1.0e-10)


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_jacobsen_2002_independent_error_weighted_refits_have_parity(
    material_identifier, case
):
    rows = _rows(case)
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    pressure_sigma = np.array([float(row["pressure_sd_gpa"] or 0.0) for row in rows])
    lattice_a = np.array([float(row["lattice_a_angstrom"]) for row in rows])
    lattice_a_sigma = np.array([float(row["lattice_a_sd_angstrom"]) for row in rows])
    volume = lattice_a**3
    volume_sigma = 3.0 * lattice_a**2 * lattice_a_sigma

    def effective_variance_residual(parameters):
        step = 1.0e-5
        pressure_derivative = (
            _bm3(parameters, volume + step) - _bm3(parameters, volume - step)
        ) / (2.0 * step)
        sigma = np.sqrt(pressure_sigma**2 + (pressure_derivative * volume_sigma) ** 2)
        return (_bm3(parameters, volume) - pressure) / sigma

    result = least_squares(
        effective_variance_residual,
        x0=case["parameters"],
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
        max_nfev=10000,
    )

    assert result.success
    # Tight optimizer tolerances expose small SciPy/platform differences because
    # the effective variance itself depends on the fitted parameters.
    assert result.x == pytest.approx(case["refit"], abs=5.0e-4)
    assert np.all(
        np.abs(result.x - np.asarray(case["parameters"])) < np.asarray(case["errors"])
    )
    stored = _document(material_identifier)["eos_records"][0]["scientific_validation"][
        "independent_refit"
    ]["parameters"]
    assert tuple(stored.values()) == pytest.approx(result.x, abs=5.0e-4)


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_jacobsen_2002_quartz_calibration_is_numerically_auditable(
    material_identifier, case
):
    rows = _rows(case)
    quartz_volume = np.array([float(row["quartz_volume_a3"]) for row in rows])
    source_pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    # The run's near-ambient quartz observation removes the cell-specific
    # reference-volume offset; K0 and K0' remain those of Angel et al. (1997).
    recalculated = _bm3((quartz_volume[0], 37.12, 5.99), quartz_volume)
    residual = recalculated - source_pressure

    assert np.sqrt(np.mean(residual**2)) == pytest.approx(
        case["quartz_rmse"], abs=1.0e-10
    )
    assert np.max(np.abs(residual)) == pytest.approx(
        case["quartz_max_residual"], abs=1.0e-10
    )
    validation = _document(material_identifier)["eos_records"][0][
        "scientific_validation"
    ]["numerical_reproduction"]
    assert validation["quartz_recalculation_max_abs_residual_gpa"] == pytest.approx(
        np.max(np.abs(residual)), abs=1.0e-10
    )
