import csv
import hashlib
import io
import re
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS
from peritheos.fitting import fit_joint_eos
from peritheos.materials import Material

MATERIAL_ID = "e_feooh"
RECORD_ID = "e_feooh_suzuki_2016_bm3_thermal_2"
DATASET_ID = "epsilon_feooh_suzuki_2016_table1_pvt"


def _document_source_and_executable():
    document = get_material_document(MATERIAL_ID)
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == RECORD_ID
    )
    executable = Material.from_eosmat(
        document, record_identifiers=[RECORD_ID]
    ).eos_records[0]
    return document, source, executable


def _dataset_rows(document):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_bytes()
    )
    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    return dataset, list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))


def _source_pressure(volume, temperature):
    volume = np.asarray(volume, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    delta = temperature - 300.0
    v0_t = 66.278 * np.exp(2.6e-5 * delta + 0.5e-7 * delta**2)
    k0_t = 135.0 - 0.05 * delta
    eta = (v0_t / volume) ** (1.0 / 3.0)
    return 1.5 * k0_t * (eta**7 - eta**5) * (1.0 + 0.75 * (6.1 - 4.0) * (eta**2 - 1.0))


def test_suzuki_material_identity_is_distinct_and_diffraction_ready():
    document, source, _ = _document_source_and_executable()
    high_pressure = get_material_document("e_feooh_hc_low_spin")

    assert document["phase"] == "hydrogen-off-center P21nm epsilon-FeOOH"
    assert document["space_group"] == "P21nm"
    assert document["space_group_number"] == 31
    assert document["formula_units_per_cell"] == 2
    assert source["reference"]["doi"] == "10.2465/jmps.160719c"
    assert high_pressure["phase"] == "hydrogen-centered Pnnm low-spin epsilon-FeOOH"

    counts = {"Fe": 0.0, "O": 0.0, "H": 0.0}
    for site in document["atom_sites"]:
        multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
        counts[site["element"]] += multiplicity * site["occupancy"]
    assert counts == pytest.approx({"Fe": 2.0, "O": 4.0, "H": 2.0})
    assert (
        "e_feooh_hc_low_spin"
        in source["scientific_validation"]["phase_boundary"]["finding"]
    )


def test_suzuki_published_parameters_and_reference_temperature_convention():
    _, source, record = _document_source_and_executable()

    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 66.278, "K0": 135.0, "K0_prime": 6.1},
    }
    assert source["parameter_errors"] == {
        "V0": 0.006,
        "K0": 3.0,
        "K0_prime": 0.9,
    }
    assert source["fixed_parameters"] == ["V0"]
    assert source["parameter_error_confidence"] is None
    assert source["thermal"] == {
        "type": "AlphaKT",
        "model": "thermal_reference_state",
        "thermal_expansion_law": "linear_reference_temperature",
        "reference_volume_law": "integrated_expansivity",
        "parameters": {
            "Tr": 300.0,
            "alpha0": 2.6e-5,
            "alpha1": 1.0e-7,
            "dK_dT": -0.05,
        },
        "parameter_errors": {
            "Tr": None,
            "alpha0": 0.7e-5,
            "alpha1": 0.3e-7,
            "dK_dT": 0.02,
        },
        "fixed_parameters": ["Tr"],
    }

    model = record.eos
    assert model.pressure(66.278, 300.0) == pytest.approx(0.0, abs=1.0e-12)
    assert model.bulk_modulus(66.278, 300.0) == pytest.approx(135.0)
    assert model.pressure(62.63, 700.0) == pytest.approx(10.833562803310029)
    for pressure, temperature in ((0.68, 300.0), (9.07, 600.0), (11.07, 700.0)):
        volume = record.volume(pressure, temperature, check_validity=True)
        assert record.pressure(
            volume, temperature, check_validity=True
        ) == pytest.approx(pressure, rel=1.0e-11)
    for temperature in (300.0, 500.0, 700.0):
        delta = temperature - 300.0
        zero_pressure_volume = 66.278 * np.exp(2.6e-5 * delta + 0.5e-7 * delta**2)
        assert model.pressure(zero_pressure_volume, temperature) == pytest.approx(
            0.0, abs=2.0e-12
        )


def test_suzuki_table1_transcription_and_published_curve_reproduction():
    document, source, record = _document_source_and_executable()
    dataset, rows = _dataset_rows(document)

    assert len(rows) == 33
    assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
    assert rows[0] == {
        "source_row": "1",
        "pressure_gpa": "0.0001",
        "pressure_sigma_gpa": "",
        "temperature_k": "300",
        "lattice_a_angstrom": "4.9544",
        "lattice_a_sigma_angstrom": "0.0002",
        "lattice_b_angstrom": "4.4594",
        "lattice_b_sigma_angstrom": "0.0003",
        "lattice_c_angstrom": "2.9999",
        "lattice_c_sigma_angstrom": "0.0001",
        "cell_volume_a3": "66.278",
        "cell_volume_sigma_a3": "0.006",
    }
    assert rows[21]["pressure_gpa"] == "9.07"
    assert rows[21]["temperature_k"] == "600"
    assert rows[21]["cell_volume_a3"] == "63.08"
    assert rows[-1]["pressure_gpa"] == "11.07"
    assert rows[-1]["temperature_k"] == "700"
    assert rows[-1]["cell_volume_a3"] == "62.63"
    assert dataset["uncertainty"]["confidence"] == "one_standard_deviation"
    assert "No open data reuse license" in dataset["license"]

    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    temperature = np.array([float(row["temperature_k"]) for row in rows])
    volume = np.array([float(row["cell_volume_a3"]) for row in rows])
    independent = _source_pressure(volume, temperature)
    predicted = np.asarray(record.pressure(volume, temperature), dtype=float)
    residual = predicted - pressure

    assert predicted == pytest.approx(independent, abs=1.0e-12)
    assert np.sqrt(np.mean(residual**2)) == pytest.approx(
        0.10201978513163779, abs=1.0e-12
    )
    assert np.max(np.abs(residual)) == pytest.approx(0.24320535477927407, abs=1.0e-12)
    reproduction = source["scientific_validation"]["numerical_reproduction"]
    assert reproduction["published_curve_on_table1_rows"][
        "pressure_rmse_gpa"
    ] == pytest.approx(0.10201978513163779)


def test_suzuki_peritheos_errors_in_variables_refit_has_parameter_parity():
    document, source, _ = _document_source_and_executable()
    _, rows = _dataset_rows(document)
    selected = [row for row in rows if row["pressure_sigma_gpa"]]

    def values(name):
        return np.array([float(row[name]) for row in selected])

    result = fit_joint_eos(
        ThermalReferenceStateEOS,
        BM3,
        values("cell_volume_a3"),
        values("temperature_k"),
        values("pressure_gpa"),
        initial={
            "rt_eos.K0": 135.0,
            "rt_eos.K0_prime": 6.1,
            "alpha0": 2.6e-5,
            "alpha1": 1.0e-7,
            "dK_dT": -0.05,
        },
        fixed={"rt_eos.V0": 66.278, "Tr": 300.0},
        configuration={
            "thermal_expansion_law": "linear_reference_temperature",
            "reference_volume_law": "integrated_expansivity",
        },
        bounds={
            "rt_eos.K0": (1.0, 500.0),
            "rt_eos.K0_prime": (-5.0, 30.0),
            "alpha0": (-1.0e-4, 2.0e-4),
            "alpha1": (-1.0e-6, 1.0e-6),
            "dK_dT": (-0.5, 0.5),
        },
        pressure_sigma=values("pressure_sigma_gpa"),
        volume_sigma=values("cell_volume_sigma_a3"),
        absolute_sigma=True,
    )
    expected = source["scientific_validation"]["numerical_reproduction"][
        "diagnostic_peritheos_errors_in_variables_refit"
    ]
    assert result.parameters["rt_eos.K0"] == pytest.approx(
        expected["K0_gpa"], abs=2.0e-6
    )
    assert result.parameters["rt_eos.K0_prime"] == pytest.approx(
        expected["K0_prime"], abs=2.0e-6
    )
    assert result.parameters["alpha0"] == pytest.approx(
        expected["alpha0_per_k"], abs=2.0e-11
    )
    assert result.parameters["alpha1"] == pytest.approx(
        expected["alpha1_per_k2"], abs=2.0e-13
    )
    assert result.parameters["dK_dT"] == pytest.approx(
        expected["dK_dT_gpa_per_k"], abs=2.0e-7
    )
    assert result.reduced_chi_square == pytest.approx(
        expected["reduced_chi_square"], abs=2.0e-8
    )

    published = {
        "rt_eos.K0": (135.0, 3.0),
        "rt_eos.K0_prime": (6.1, 0.9),
        "alpha0": (2.6e-5, 0.7e-5),
        "alpha1": (1.0e-7, 0.3e-7),
        "dK_dT": (-0.05, 0.02),
    }
    for name, (value, sigma) in published.items():
        assert abs(result.parameters[name] - value) < sigma
