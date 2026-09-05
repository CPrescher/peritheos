import csv
import hashlib
import io
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import LinearThermalPressure, MieGruneisenDebye
from peritheos.eosmat import validate_pressure_calibration_references
from peritheos.fitting import fit_thermal_eos
from peritheos.materials import Material
from peritheos.units import cell_volume_to_molar_volume

WANG_RECORD_ID = "ca_perovskite_wang_1996_thermal_pressure_5"
COMBINED_RECORD_ID = "ca_perovskite_shim_2000_combined_bm3_mgd_6"
DAC_ONLY_RECORD_ID = "ca_perovskite_shim_2000_dac_only_bm3_mgd_4"
WANG_DATASET_ID = "ca_perovskite_wang_1996_table1_pvt"
SHIM_DATASET_ID = "ca_perovskite_shim_2000_jgr_table1_dac_pvt"


def _document_and_record(identifier):
    document = get_material_document("ca_perovskite")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == identifier
    )
    return document, record


def _dataset_rows(document, identifier):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == identifier
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_bytes()
    )
    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
    return dataset, rows


def test_wang_1996_is_a_distinct_published_thermal_pressure_record():
    document, record = _document_and_record(WANG_RECORD_ID)
    matching = [
        item
        for item in document["eos_records"]
        if item["reference"].get("doi", "").lower() == "10.1029/95jb03254"
    ]

    assert [item["identifier"] for item in matching] == [WANG_RECORD_ID]
    assert record["record_kind"] == "published"
    assert record["equation_kind"] == "thermal"
    assert record["default"] is False
    assert record["fit_datasets"] == [WANG_DATASET_ID]
    assert record["phase_selection"]["cell_basis"] == (
        "one-formula-unit cubic cell (Z=1)"
    )
    assert record["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 45.58, "K0": 232.0, "K0_prime": 4.8},
    }
    assert record["parameter_errors"] == {
        "V0": 0.04,
        "K0": 8.0,
        "K0_prime": 0.3,
    }
    assert record["fixed_parameters"] == ["K0_prime"]
    assert record["thermal"] == {
        "type": "LinearThermalPressure",
        "model": "linear_thermal_pressure",
        "parameters": {"Tr": 300.0, "alpha_KT": 0.0070992},
        "parameter_errors": {"Tr": None, "alpha_KT": 0.0001},
        "fixed_parameters": ["Tr"],
    }

    executable = Material.from_eosmat(
        document, record_identifiers=[WANG_RECORD_ID]
    ).eos_records[0]
    assert executable.pressure(45.58, 300.0) == pytest.approx(0.0, abs=1.0e-12)
    assert executable.thermal_pressure_increment(44.5, 1300.0) == pytest.approx(7.0992)


def test_wang_1996_table1_is_complete_checksummed_and_excludes_only_flagged_rows():
    document, _ = _document_and_record(WANG_RECORD_ID)
    dataset, rows = _dataset_rows(document, WANG_DATASET_ID)

    assert dataset["used_by_eos_records"] == [WANG_RECORD_ID, COMBINED_RECORD_ID]
    assert len(rows) == 66
    assert [row["source_order"] for row in rows] == [str(i) for i in range(1, 67)]
    assert {
        run: sum(row["run"] == run for row in rows) for run in ("13", "3", "4", "5")
    } == {
        "13": 34,
        "3": 15,
        "4": 12,
        "5": 5,
    }
    assert sum(row["source_lineage"] == "Wang_et_al_1996" for row in rows) == 34
    assert (
        sum(
            row["source_lineage"]
            == "Wang_and_Weidner_1994_reprinted_in_Wang_et_al_1996"
            for row in rows
        )
        == 32
    )
    included = [row for row in rows if row["include_in_eos_fit"] == "true"]
    excluded = [row for row in rows if row["include_in_eos_fit"] == "false"]
    assert len(included) == 64
    assert [row["source_order"] for row in excluded] == ["33", "34"]
    assert excluded[0]["exclusion_reason"] == "high_nonhydrostatic_stress"
    assert excluded[1]["exclusion_reason"] == (
        "amorphization_and_high_nonhydrostatic_stress"
    )
    assert rows[0] == {
        "source_order": "1",
        "run": "13",
        "pressure_gpa": "12.04",
        "temperature_k": "1368",
        "volume_a3": "44.76",
        "volume_standard_deviation_a3": "0.01",
        "differential_stress_gpa": "0.001",
        "include_in_eos_fit": "true",
        "exclusion_reason": "",
        "source_lineage": "Wang_et_al_1996",
    }
    assert rows[-1]["pressure_gpa"] == "11.69"
    assert rows[-1]["temperature_k"] == "1594"
    assert rows[-1]["volume_standard_deviation_a3"] == "0.05"


def test_wang_1996_unweighted_thermal_pressure_refit_has_parameter_parity():
    document, record = _document_and_record(WANG_RECORD_ID)
    _, rows = _dataset_rows(document, WANG_DATASET_ID)
    fit_rows = [
        row
        for row in rows
        if row["run"] == "13" and row["include_in_eos_fit"] == "true"
    ]
    volumes = np.array([float(row["volume_a3"]) for row in fit_rows])
    temperatures = np.array([float(row["temperature_k"]) for row in fit_rows])
    pressures = np.array([float(row["pressure_gpa"]) for row in fit_rows])

    def residuals(parameters):
        model = LinearThermalPressure(
            BM3(V0=parameters[0], K0=parameters[1], K0_prime=4.8),
            Tr=300.0,
            alpha_KT=parameters[2],
        )
        return np.asarray(model.pressure(volumes, temperatures)) - pressures

    result = least_squares(residuals, x0=(45.58, 232.0, 0.0070992))
    degrees_of_freedom = len(pressures) - len(result.x)
    covariance = np.linalg.inv(result.jac.T @ result.jac) * (
        np.sum(result.fun**2) / degrees_of_freedom
    )
    standard_errors = np.sqrt(np.diag(covariance))

    assert result.success
    assert result.x == pytest.approx(
        [45.5812586, 231.372767, 0.00707418205], rel=1.0e-8
    )
    assert standard_errors == pytest.approx(
        [0.025982499, 4.49884021, 0.000131288893], rel=1.0e-7
    )
    assert abs(result.x[0] - 45.58) < 0.04
    assert abs(result.x[1] - 232.0) < 8.0
    assert abs(result.x[2] - 0.0071) < 0.0001

    executable = Material.from_eosmat(
        document, record_identifiers=[WANG_RECORD_ID]
    ).eos_records[0]
    published_residuals = (
        np.asarray(executable.pressure(volumes, temperatures)) - pressures
    )
    assert np.sqrt(np.mean(published_residuals**2)) == pytest.approx(0.1979129256)
    assert np.sqrt(np.mean(result.fun**2)) == pytest.approx(0.1969957061)
    assert record["scientific_validation"]["independent_refit"]["result"] == ("parity")


def test_shim_2000_preferred_combined_fit_reproduces_from_both_primary_tables():
    document, record = _document_and_record(COMBINED_RECORD_ID)
    _, wang_rows = _dataset_rows(document, WANG_DATASET_ID)
    _, shim_rows = _dataset_rows(document, SHIM_DATASET_ID)
    wang_rows = [row for row in wang_rows if row["include_in_eos_fit"] == "true"]
    rows = wang_rows + shim_rows
    volumes = np.array([float(row["volume_a3"]) for row in rows])
    temperatures = np.array([float(row["temperature_k"]) for row in rows])
    pressures = np.array([float(row["pressure_gpa"]) for row in rows])
    volume_scale = cell_volume_to_molar_volume(1.0, 1.0)
    reference = BM3(V0=2.745, K0=236.0, K0_prime=3.9)

    result = fit_thermal_eos(
        MieGruneisenDebye,
        reference,
        volume=volumes * volume_scale,
        temperature=temperatures,
        pressure=pressures,
        initial={"gamma0": 1.92, "q": 0.6},
        fixed={"Tr": 300.0, "theta0": 1000.0, "n": 5.0},
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds={"gamma0": (1.0e-9, 10.0), "q": (-20.0, 20.0)},
        max_nfev=5000,
    )

    assert result.success
    assert len(rows) == 98
    assert [result.parameters[name] for name in ("gamma0", "q")] == pytest.approx(
        [1.9188965478, 0.5889618399], rel=1.0e-9
    )
    assert [result.standard_errors[name] for name in ("gamma0", "q")] == pytest.approx(
        [0.0327741346, 0.1963589263], rel=1.0e-8
    )
    assert abs(result.parameters["gamma0"] - 1.92) < 0.05
    assert abs(result.parameters["q"] - 0.6) < 0.3

    executable = Material.from_eosmat(
        document, record_identifiers=[COMBINED_RECORD_ID]
    ).eos_records[0]
    published_residuals = (
        np.asarray(executable.pressure(volumes, temperatures)) - pressures
    )
    refit_residuals = (
        np.asarray(result.model.pressure(volumes * volume_scale, temperatures))
        - pressures
    )
    assert np.sqrt(np.mean(published_residuals**2)) == pytest.approx(0.7108254828)
    assert np.sqrt(np.mean(refit_residuals**2)) == pytest.approx(0.7108117647)
    assert record["scientific_validation"]["independent_refit"]["result"] == ("parity")
    assert (
        record["scientific_validation"]["related_alternative_fit"]["identifier"]
        == DAC_ONLY_RECORD_ID
    )


def test_wang_and_combined_pressure_calibration_links_resolve_globally():
    validate_pressure_calibration_references()
