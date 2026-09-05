import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_joint_eos

ROOT = Path(__file__).parents[1]
DOI = "10.1029/2023JB028026"
RECORD_IDENTIFIER = "magnesite_yu_2024_bm3_mgd_3"
DATASET_IDENTIFIER = "magnesite_yu_2024_table_s1_pvt"
DATASET_SHA256 = "a3c5ef81403988ed368d62e03712323032276e078a70826fcdd161d26f1b5d72"


def _source_dataset_rows():
    document = get_material_document("magnesite")
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == RECORD_IDENTIFIER
    )
    dataset = next(
        item
        for item in document["datasets"]
        if item["identifier"] == DATASET_IDENTIFIER
    )
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, source, dataset, path, rows


def test_yu_2024_record_preserves_model_identity_and_provenance():
    document, source, dataset, _, _ = _source_dataset_rows()

    assert document["formula"] == "MgCO3"
    assert document["phase"] == "magnesite, rhombohedral R-3c"
    assert document["space_group"] == "R-3c"
    assert document["formula_units_per_cell"] == 6
    assert source["reference"]["doi"] == DOI
    assert source["record_kind"] == "published"
    assert source["default"] is True
    assert source["phase_scope"] == (
        "rhombohedral magnesite, space group R-3c, conventional hexagonal cell with Z=6"
    )
    assert source["sample_composition"]["nominal_formula"] == "MgCO3"
    assert (
        "less than 0.5 mol.% Mn and Fe"
        in source["sample_composition"]["reported_impurities"]
    )

    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 280.71, "K0": 105.0, "K0_prime": 4.49},
    }
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": 2.0,
        "K0_prime": 0.05,
    }
    assert source["fixed_parameters"] == ["V0"]
    assert source["parameter_error_confidence"] is None
    assert source["parameter_covariance"] is None
    assert source["temperature_ref"] == 300.0
    assert source["thermal"] == {
        "type": "MieGruneisenDebye",
        "model": "mie_gruneisen_debye",
        "debye_temperature_law": "integrated_gruneisen",
        "parameters": {
            "Tr": 300.0,
            "theta0": 850.0,
            "gamma0": 1.7,
            "q": 0.9,
            "n": 5.0,
        },
        "parameter_errors": {
            "Tr": None,
            "theta0": 200.0,
            "gamma0": 0.2,
            "q": 0.3,
            "n": None,
        },
        "fixed_parameters": ["Tr", "n"],
    }
    assert source["fit_datasets"] == [DATASET_IDENTIFIER]
    assert source["experimental_pressure_range_gpa"] == [0.0, 119.5]
    assert source["experimental_temperature_range_k"] == [300.0, 2200.0]

    pressure_method = source["pressure_calibration"]["methods"][0]
    assert pressure_method["material"] == "Pt"
    assert pressure_method["reference"]["doi"] == "10.1073/pnas.0609013104"
    assert "variable-exponent" in pressure_method["notes"]
    assert (
        source["pressure_calibration"]["recalculation"]["status"]
        == "missing_calibrant_observations"
    )
    assert dataset["used_by_eos_records"] == [RECORD_IDENTIFIER]


def test_yu_2024_table_s1_transcription_checksum_and_source_anomaly():
    _, source, dataset, path, rows = _source_dataset_rows()

    assert hashlib.sha256(path.read_bytes()).hexdigest() == DATASET_SHA256
    assert dataset["resource"]["sha256"] == DATASET_SHA256
    assert dataset["license"] == "CC BY 4.0"
    assert dataset["source_url"] == "https://doi.org/10.5281/zenodo.8399354"
    assert dataset["source_artifact"] == {
        "filename": "data.docx",
        "download_url": (
            "https://zenodo.org/api/records/8399354/files/data.docx/content"
        ),
        "sha256": ("9eed9b3a3c0573e1f9c7fcd29d4689bbd4c479259e2ea92e918ca688aceae1f0"),
        "repository_checksum": "md5:b41936980d505f5d6e337806f476fa57",
        "retrieved": "2026-09-05",
    }
    assert len(rows) == 74
    assert len(rows[0]) == 11
    assert rows[0] == {
        "source_order": "1",
        "pressure_gpa": "0",
        "pressure_uncertainty_gpa": "",
        "temperature_k": "300",
        "temperature_uncertainty_k": "",
        "a_angstrom": "4.6431",
        "a_uncertainty_angstrom": "",
        "c_angstrom": "15.0351",
        "c_uncertainty_angstrom": "",
        "volume_a3": "280.71",
        "volume_uncertainty_a3": "",
    }
    assert rows[-1]["source_order"] == "74"
    assert rows[-1]["pressure_gpa"] == "95"
    assert rows[-1]["temperature_k"] == "2200"
    assert rows[-1]["volume_a3"] == "198.53"

    lattice_volume_residuals = []
    for row in rows:
        lattice_volume = (
            math.sqrt(3.0)
            / 2.0
            * float(row["a_angstrom"]) ** 2
            * float(row["c_angstrom"])
        )
        lattice_volume_residuals.append(lattice_volume - float(row["volume_a3"]))
    anomalous = np.flatnonzero(np.abs(lattice_volume_residuals) > 0.3)
    assert anomalous.tolist() == [4]
    assert lattice_volume_residuals[4] == pytest.approx(4.8877528799)
    assert rows[4]["source_order"] == "5"
    assert rows[4]["c_angstrom"] == "14.310"
    assert rows[4]["volume_a3"] == "246.96"
    assert "preserved verbatim" in dataset["notes"]
    assert (
        "source-order row 5"
        in source["scientific_validation"]["primary_data_check"]["finding"]
    )


def test_yu_2024_published_mgd_curve_and_volume_conversion():
    document, _, _, _, rows = _source_dataset_rows()
    material = Material.from_eosmat(document, record_identifiers=[RECORD_IDENTIFIER])
    record = material.eos_records[0]
    volumes = np.array([float(row["volume_a3"]) for row in rows])
    temperatures = np.array([float(row["temperature_k"]) for row in rows])
    pressures = np.array([float(row["pressure_gpa"]) for row in rows])

    assert isinstance(record.eos, MieGruneisenDebye)
    assert isinstance(record.eos.rt_eos, BM3)
    assert record.eos.debye_temperature_law == "integrated_gruneisen"
    assert record.reference_volume == pytest.approx(280.71)
    assert record.volume_scale == pytest.approx(0.010036901266666667)
    assert record.eos.rt_eos.V0 == pytest.approx(280.71 * record.volume_scale)
    assert record.pressure(280.71, 300.0) == pytest.approx(0.0, abs=1.0e-12)

    calculated = np.asarray(
        record.pressure(volumes, temperatures, check_validity=False), dtype=float
    )
    residuals = calculated - pressures
    assert np.sqrt(np.mean(residuals**2)) == pytest.approx(
        1.0304556630476074, abs=2.0e-10
    )
    assert np.max(np.abs(residuals)) == pytest.approx(2.8019396360871553, abs=2.0e-9)

    volume = record.volume(80.0, 1800.0, check_validity=True)
    assert record.pressure(volume, 1800.0, check_validity=True) == pytest.approx(
        80.0, abs=1.0e-9
    )


def test_yu_2024_independent_peritheos_refit_has_parameter_parity():
    document, source, _, _, rows = _source_dataset_rows()
    record = Material.from_eosmat(
        document, record_identifiers=[RECORD_IDENTIFIER]
    ).eos_records[0]
    volumes = np.array([float(row["volume_a3"]) for row in rows])
    temperatures = np.array([float(row["temperature_k"]) for row in rows])
    pressures = np.array([float(row["pressure_gpa"]) for row in rows])

    result = fit_joint_eos(
        MieGruneisenDebye,
        BM3,
        volumes * record.volume_scale,
        temperatures,
        pressures,
        initial={
            "rt_eos.K0": 105.0,
            "rt_eos.K0_prime": 4.49,
            "theta0": 850.0,
            "gamma0": 1.7,
            "q": 0.9,
        },
        fixed={
            "rt_eos.V0": 280.71 * record.volume_scale,
            "Tr": 300.0,
            "n": 5.0,
        },
        bounds={
            "rt_eos.K0": (50.0, 200.0),
            "rt_eos.K0_prime": (0.0, 10.0),
            "theta0": (100.0, 3000.0),
            "gamma0": (0.1, 4.0),
            "q": (-2.0, 5.0),
        },
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        max_nfev=5000,
    )

    expected = {
        "rt_eos.K0": 105.63974787587165,
        "rt_eos.K0_prime": 4.477101605441763,
        "theta0": 657.5642321101462,
        "gamma0": 1.9735760309784676,
        "q": 1.3229215647013208,
    }
    assert result.success
    assert result.parameters == pytest.approx(
        {
            "rt_eos.V0": 280.71 * record.volume_scale,
            "Tr": 300.0,
            "n": 5.0,
            **expected,
        },
        abs=5.0e-5,
    )
    stored = source["scientific_validation"]["independent_refit"]
    assert stored["implementation"] == "peritheos.fitting.fit_joint_eos"
    assert stored["selection"] == (
        "all 74 official Table S1 rows, using the printed volume column"
    )
    assert stored["objective"] == "unweighted pressure residuals"
    assert stored["parameters"] == pytest.approx(expected, abs=2.0e-6)
    assert stored["observations"] == 74
    assert stored["degrees_of_freedom"] == 69

    published = {
        "rt_eos.K0": 105.0,
        "rt_eos.K0_prime": 4.49,
        "theta0": 850.0,
        "gamma0": 1.7,
        "q": 0.9,
    }
    published_errors = {
        "rt_eos.K0": 2.0,
        "rt_eos.K0_prime": 0.05,
        "theta0": 200.0,
        "gamma0": 0.2,
        "q": 0.3,
    }
    normalized_differences = {
        name: abs(result.parameters[name] - published[name]) / published_errors[name]
        for name in published
    }
    assert max(normalized_differences.values()) == pytest.approx(
        1.4097385490, abs=2.0e-8
    )
    assert all(value < 2.0 for value in normalized_differences.values())

    refit_pressures = np.asarray(
        result.model.pressure(volumes * record.volume_scale, temperatures),
        dtype=float,
    )
    refit_residuals = refit_pressures - pressures
    assert np.sqrt(np.mean(refit_residuals**2)) == pytest.approx(
        stored["refit_pressure_rmse_gpa"], abs=2.0e-9
    )
    assert np.max(np.abs(refit_residuals)) == pytest.approx(
        stored["refit_maximum_absolute_residual_gpa"], abs=2.0e-8
    )
    assert result.correlation[0, 1] == pytest.approx(-0.9769620368, abs=2.0e-8)
