import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.eos.rt import BM2
from peritheos.eosmat import validate_pressure_calibration_references
from peritheos.fitting import fit_rt_eos

ROOT = Path(__file__).parents[1]
MATERIAL_ID = "mg090fe009al0005ca0005sio3_post_perovskite"
RECORD_ID = f"{MATERIAL_ID}_shieh_2006_bm2_1"
DATASET_ID = f"{MATERIAL_ID}_shieh_2006_figure2_digitized"
DOI = "10.1073/pnas.0506811103"
RESOURCE = "mg090fe009al0005ca0005sio3-post-perovskite-shieh-2006-figure2-digitized.csv"
SHA256 = "9f1e1ce35ee0f56b950dd1f8b4875681db947ba61f69e91a2f38a005a4792382"


def _source_and_rows():
    document = get_material_document(MATERIAL_ID)
    source = document["eos_records"][0]
    dataset = document["datasets"][0]
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, source, dataset, path, rows


def _record(document):
    return Material.from_eosmat(
        document, record_identifiers=[RECORD_ID]
    ).get_eos_record(RECORD_ID)


def test_shieh_natural_post_perovskite_is_composition_specific_and_distinct():
    document, source, _, _, _ = _source_and_rows()
    pure = get_material_document("mgsio3_post_perovskite")

    assert document["identifier"] == MATERIAL_ID
    assert document["formula"] == "Mg0.90Fe0.09Al0.005Ca0.005SiO3"
    assert "not measured separately" in document["phase"]
    assert document["space_group"] == "Cmcm"
    assert document["space_group_number"] == 63
    assert document["formula_units_per_cell"] == 4
    assert "atom_sites" not in document
    assert source["sample_composition"].startswith(
        "Natural orthopyroxene starting material"
    )
    assert source["reference"]["doi"].lower() == DOI
    assert not any(
        record["reference"]["doi"].lower() == DOI for record in pure["eos_records"]
    )

    lattice = document["lattice"]
    assert lattice["a"] * lattice["b"] * lattice["c"] == pytest.approx(124.569397008)
    assert "bulk-composition label" in document["cell_contents"]


def test_shieh_figure2_complete_digitization_checksum_and_pixel_calibration():
    _, source, dataset, path, rows = _source_and_rows()

    assert dataset["identifier"] == DATASET_ID
    assert path.name == RESOURCE
    assert hashlib.sha256(path.read_bytes()).hexdigest() == SHA256
    assert dataset["resource"]["sha256"] == SHA256
    assert len(rows) == 14
    assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
    assert all(row["used_in_published_fit"] == "1" for row in rows)
    assert [float(rows[index]["pressure_gpa"]) for index in (0, 2, 3, -1)] == [
        12.036,
        32.046,
        83.028,
        105.911,
    ]

    axes = dataset["provenance"]["axis_calibration_pixels"]
    pressure_mapping = np.polyfit(*np.asarray(axes["pressure"], dtype=float).T, 1)
    volume_mapping = np.polyfit(*np.asarray(axes["volume"], dtype=float).T, 1)
    for row in rows:
        pressure = np.polyval(pressure_mapping, float(row["source_marker_x_pixel"]))
        volume = np.polyval(volume_mapping, float(row["source_marker_y_pixel"]))
        assert float(row["pressure_gpa"]) == pytest.approx(pressure, abs=5.1e-4)
        assert float(row["unit_cell_volume_a3"]) == pytest.approx(volume, abs=5.1e-4)

    highest = source["scientific_validation"]["numerical_reproduction"][
        "highest_pressure_cross_check"
    ]
    assert highest["figure_2_digitized_volume_a3"] == float(
        rows[-1]["unit_cell_volume_a3"]
    )
    assert highest["absolute_difference_a3"] < 0.002


def test_shieh_published_bm2_parameters_and_derived_bulk_moduli():
    document, source, _, _, rows = _source_and_rows()
    record = _record(document)

    assert source["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 164.9, "K0": 219.0},
    }
    assert source["parameter_errors"] == {"V0": 0.6, "K0": 5.0}
    assert source["parameter_error_confidence"] == pytest.approx(0.682689492137)
    assert source["parameter_covariance"] is None
    assert source["fixed_parameters"] == []
    assert source["experimental_pressure_range_gpa"] == [12.0, 106.0]
    assert source["fit_datasets"] == [DATASET_ID]

    volume = np.array([float(row["unit_cell_volume_a3"]) for row in rows])
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    residual = np.asarray(record.pressure(volume), dtype=float) - pressure
    assert math.sqrt(float(np.mean(residual**2))) == pytest.approx(
        1.2881240230, abs=5.0e-10
    )
    assert float(np.max(np.abs(residual))) == pytest.approx(2.7496079338, abs=5.0e-10)

    checks = {
        100.0: (575.0, 15.0, 575.2079434686),
        120.0: (641.0, 17.0, 641.1518061148),
    }
    stored = source["scientific_validation"]["numerical_reproduction"]
    for pressure_gpa, (published, uncertainty, calculated) in checks.items():
        volume_a3 = record.volume(pressure_gpa)
        peritheos_bulk_modulus = record.eos.bulk_modulus(volume_a3)
        assert peritheos_bulk_modulus == pytest.approx(calculated, abs=5.0e-10)
        assert abs(peritheos_bulk_modulus - published) < uncertainty
        key = f"bulk_modulus_at_{int(pressure_gpa)}_gpa"
        assert stored[key]["peritheos_curve_gpa"] == pytest.approx(
            peritheos_bulk_modulus
        )


def test_shieh_digitized_bm2_refit_has_published_parameter_parity():
    _, source, _, _, rows = _source_and_rows()
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["unit_cell_volume_a3"]) for row in rows])

    fit = fit_rt_eos(
        BM2,
        volume,
        pressure,
        {"V0": 164.9, "K0": 219.0},
        bounds={"V0": (150.0, 180.0), "K0": (100.0, 400.0)},
    )
    assert fit.success
    assert fit.parameters["V0"] == pytest.approx(164.60842068, abs=5.0e-8)
    assert fit.parameters["K0"] == pytest.approx(221.62625098, abs=5.0e-8)
    assert fit.standard_errors["V0"] == pytest.approx(0.73432093, abs=5.0e-8)
    assert fit.standard_errors["K0"] == pytest.approx(6.01516185, abs=5.0e-8)
    assert fit.correlation[0, 1] == pytest.approx(-0.9872470560, abs=5.0e-10)
    assert abs(fit.parameters["V0"] - 164.9) < 0.6
    assert abs(fit.parameters["K0"] - 219.0) < 5.0

    stored = source["scientific_validation"]["independent_refit"]
    assert stored["classification"] == "parity"
    assert stored["parameters"]["V0"] == pytest.approx(fit.parameters["V0"])
    assert stored["parameters"]["K0"] == pytest.approx(fit.parameters["K0"])
    assert stored["pressure_rmse_gpa"] == pytest.approx(
        math.sqrt(float(np.mean(fit.residuals**2)))
    )


def test_shieh_pressure_scale_and_unresolved_source_limits_are_explicit():
    _, source, dataset, _, _ = _source_and_rows()
    calibration = source["pressure_calibration"]
    method = calibration["methods"][0]

    assert calibration["status"] == "resolved"
    assert method["reference_eos_record"] == "platinum_holmes_1989_vinet_1"
    assert calibration["recalculation"]["status"] == ("missing_calibrant_observations")
    assert source["scientific_validation"]["primary_data_check"]["status"] == (
        "plot_only"
    )
    assert dataset["uncertainty"]["digitization_only"] is True
    assert len(source["scientific_validation"]["unresolved_issues"]) == 5
    validate_pressure_calibration_references()
