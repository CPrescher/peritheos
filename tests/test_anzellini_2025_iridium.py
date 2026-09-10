import csv
import hashlib
import io
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import HollandPowellThermalPressure
from peritheos.fitting import fit_joint_eos
from peritheos.materials import Material

RECORD_ID = "iridium_anzellini_2025_bm3_1"
DATASET_ID = "iridium_anzellini_2025_tables_s1_s3_pvt"
DATASET_SHA256 = "2133efc7b402be8d804d098ef8badfe0e46d8a7ffc91351ca71bffaf6fc85250"
PLOT_DATASET_ID = "iridium_anzellini_2025_figure4_300k_830k_vector_digitized"
PLOT_DATASET_SHA256 = "097f35f293d7f6d93d286e86ed72166e7660817132e481b43e3301cac2d7bd8d"


def _source_record_and_rows():
    document = get_material_document("iridium")
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == RECORD_ID
    )
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_bytes()
    )
    assert hashlib.sha256(payload).hexdigest() == DATASET_SHA256
    assert dataset["resource"]["sha256"] == DATASET_SHA256
    rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
    return document, source, dataset, rows


def _source_pressure(volume, temperature):
    reference = BM3(V0=56.62, K0=327.0, K0_prime=5.46)
    temperature = np.asarray(temperature, dtype=float)
    theta = 298.0
    ratio = theta / 300.0
    zeta0 = ratio**2 * np.exp(ratio) / np.expm1(ratio) ** 2
    thermal = (
        1.87e-5
        * 327.0
        * theta
        / zeta0
        * (1.0 / np.expm1(theta / temperature) - 1.0 / np.expm1(ratio))
    )
    return reference.pressure(volume) + thermal


def _plot_dataset_and_rows(document):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == PLOT_DATASET_ID
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_bytes()
    )
    assert hashlib.sha256(payload).hexdigest() == PLOT_DATASET_SHA256
    assert dataset["resource"]["sha256"] == PLOT_DATASET_SHA256
    return dataset, list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))


def test_iridium_record_is_the_published_bm3_holland_powell_surface():
    document, source, _, _ = _source_record_and_rows()
    assert source["equation_kind"] == "thermal"
    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 56.62, "K0": 327.0, "K0_prime": 5.46},
    }
    assert source["thermal"] == {
        "type": "HollandPowellThermalPressure",
        "model": "holland_powell_thermal_pressure",
        "parameters": {"Tr": 300.0, "theta": 298.0, "alpha0": 1.87e-5, "n": 1.0},
        "parameter_errors": {"alpha0": 8e-7},
        "fixed_parameters": ["Tr", "theta", "n"],
    }

    executable = Material.from_eosmat(
        document, record_identifiers=[RECORD_ID]
    ).eos_records[0]
    assert isinstance(executable.eos, HollandPowellThermalPressure)
    volume = np.array([56.62, 52.0, 48.0])
    temperature = np.array([300.0, 2500.0, 5000.0])
    assert executable.pressure(volume, temperature) == pytest.approx(
        _source_pressure(volume, temperature), abs=2.0e-12
    )


def test_iridium_supplement_preserves_both_temperatures_and_source_mean():
    _, _, dataset, rows = _source_record_and_rows()
    assert len(rows) == 122
    assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
    upstream = np.array([float(row["upstream_temperature_k"]) for row in rows])
    downstream = np.array([float(row["downstream_temperature_k"]) for row in rows])
    mean = np.array([float(row["mean_temperature_k"]) for row in rows])
    assert mean == pytest.approx((upstream + downstream) / 2.0, abs=0.0)
    assert np.count_nonzero(mean == 300.0) == 2
    assert (float(np.min(mean)), float(np.max(mean))) == (300.0, 5560.5)
    assert "temperature_uncertainty" not in dataset
    assert "no uncertainty is invented" in dataset["notes"]


def test_iridium_conditional_hot_data_fit_recovers_alpha0():
    _, _, _, rows = _source_record_and_rows()
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    temperature = np.array([float(row["mean_temperature_k"]) for row in rows])
    lattice = np.array([float(row["iridium_lattice_a_angstrom"]) for row in rows])
    volume = lattice**3
    result = fit_joint_eos(
        HollandPowellThermalPressure,
        BM3,
        volume,
        temperature,
        pressure,
        initial={"alpha0": 1.87e-5},
        fixed={
            "rt_eos.V0": 56.62,
            "rt_eos.K0": 327.0,
            "rt_eos.K0_prime": 5.46,
            "Tr": 300.0,
            "theta": 298.0,
            "n": 1.0,
        },
        bounds={"alpha0": (-3.74e-4, 3.74e-4)},
        absolute_sigma=False,
        max_nfev=5000,
    )
    residuals = np.asarray(result.residuals, dtype=float)
    combined_two_sigma = 2.0 * np.hypot(8.0e-7, result.standard_errors["alpha0"])

    assert result.parameters["alpha0"] == pytest.approx(
        1.7740295689829463e-5, abs=2.0e-12
    )
    assert np.sqrt(np.mean(residuals**2)) == pytest.approx(
        3.0325674077570715, abs=2.0e-9
    )
    assert abs(result.parameters["alpha0"] - 1.87e-5) <= combined_two_sigma


def test_iridium_figure4_vector_digitization_is_calibrated_and_separate():
    document, _, _, _ = _source_record_and_rows()
    dataset, rows = _plot_dataset_and_rows(document)
    temperatures = np.array([float(row["temperature_k"]) for row in rows])

    assert len(rows) == 107
    assert np.count_nonzero(temperatures == 300.0) == 91
    assert np.count_nonzero(temperatures == 830.0) == 16
    assert dataset["provenance"]["type"] == "digitized_from_figure"
    assert dataset["provenance"]["source_pdf_page"] == 3
    assert "not author-reported" in dataset["provenance"]["digitization_uncertainty"]

    calibration = dataset["provenance"]["axis_calibration_pdf_points"]
    first = rows[0]
    x0, p0 = calibration["pressure"][0]
    x1, p1 = calibration["pressure"][1]
    pressure = p0 + (float(first["plot_x_pt"]) - x0) * (p1 - p0) / (x1 - x0)
    y0, v0 = calibration["volume"][0]
    y1, v1 = calibration["volume"][1]
    volume = v0 + (float(first["plot_y_pt"]) - y0) * (v1 - v0) / (y1 - y0)
    assert float(first["pressure_gpa"]) == pytest.approx(pressure, abs=1.0e-6)
    assert float(first["volume_a3_conventional_cell"]) == pytest.approx(
        volume, abs=1.0e-6
    )


def test_iridium_plot_derived_combined_fit_approaches_published_coefficients():
    document, _, _, hot_rows = _source_record_and_rows()
    _, plot_rows = _plot_dataset_and_rows(document)

    plot_pressure = np.array([float(row["pressure_gpa"]) for row in plot_rows])
    plot_volume = np.array(
        [float(row["volume_a3_conventional_cell"]) for row in plot_rows]
    )
    plot_temperature = np.array([float(row["temperature_k"]) for row in plot_rows])
    hot_pressure = np.array([float(row["pressure_gpa"]) for row in hot_rows])
    hot_volume = np.array(
        [float(row["iridium_lattice_a_angstrom"]) ** 3 for row in hot_rows]
    )
    hot_temperature = np.array([float(row["mean_temperature_k"]) for row in hot_rows])
    centers = np.array([1700.0, 2000.0, 2200.0, 2600.0, 3300.0, 3800.0, 4200.0])
    selected = (
        np.min(np.abs(hot_temperature[:, None] - centers[None, :]), axis=1) <= 100.0
    )
    assert np.count_nonzero(selected) == 65

    result = fit_joint_eos(
        HollandPowellThermalPressure,
        BM3,
        np.concatenate([plot_volume, hot_volume[selected]]),
        np.concatenate([plot_temperature, hot_temperature[selected]]),
        np.concatenate([plot_pressure, hot_pressure[selected]]),
        initial={
            "rt_eos.V0": 56.62,
            "rt_eos.K0": 327.0,
            "rt_eos.K0_prime": 5.46,
            "alpha0": 1.87e-5,
        },
        fixed={"Tr": 300.0, "theta": 298.0, "n": 1.0},
        bounds={
            "rt_eos.V0": (14.155, 113.24),
            "rt_eos.K0": (16.35, 1635.0),
            "rt_eos.K0_prime": (0.0, 20.0),
            "alpha0": (-3.74e-4, 3.74e-4),
        },
        absolute_sigma=False,
        max_nfev=5000,
    )
    residuals = np.asarray(result.residuals, dtype=float)
    assert result.parameters["rt_eos.V0"] == pytest.approx(
        56.62669063796952, abs=2.0e-8
    )
    assert result.parameters["rt_eos.K0"] == pytest.approx(
        326.15545648530644, abs=2.0e-6
    )
    assert result.parameters["rt_eos.K0_prime"] == pytest.approx(
        5.581046756619618, abs=5.0e-8
    )
    assert result.parameters["alpha0"] == pytest.approx(
        1.829165666059566e-5, abs=2.0e-12
    )
    assert np.sqrt(np.mean(residuals**2)) == pytest.approx(
        1.9858345098349754, abs=2.0e-9
    )


def test_iridium_provenance_does_not_claim_a_complete_combined_refit():
    _, source, _, _ = _source_record_and_rows()
    provenance = source["fit_provenance"]
    assert provenance["software"] == {"name": "EosFit7", "version": "not reported"}
    assert provenance["refined_parameters"] == ["V0", "K0", "K0_prime", "alpha0"]
    assert provenance["reproduction"]["refined_parameters"] == ["alpha0"]
    assert (
        "not an exact reproduction of the complete combined fit"
        in provenance["reproduction"]["status"]
    )
    assert provenance["reproduction"]["plot_derived_combined_fit"][
        "refined_parameters"
    ] == [
        "V0",
        "K0",
        "K0_prime",
        "alpha0",
    ]
    assert any(
        "original numerical rows" in item for item in provenance["unavailable_details"]
    )
