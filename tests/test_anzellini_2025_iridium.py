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


def test_iridium_provenance_does_not_claim_a_complete_combined_refit():
    _, source, _, _ = _source_record_and_rows()
    provenance = source["fit_provenance"]
    assert provenance["software"] == {"name": "EosFit7", "version": "not reported"}
    assert provenance["refined_parameters"] == ["V0", "K0", "K0_prime", "alpha0"]
    assert provenance["reproduction"]["refined_parameters"] == ["alpha0"]
    assert (
        "not a reproduction of the complete combined fit"
        in provenance["reproduction"]["status"]
    )
    assert any("300 K" in item for item in provenance["unavailable_details"])
