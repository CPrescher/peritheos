import csv
import hashlib
import io
import math
import re
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import get_material_document
from peritheos.materials import Material

MATERIAL_ID = "e_feooh_hc_low_spin"
RECORD_ID = "e_feooh_hc_low_spin_thompson_2017_bm3_1"
DATASET_ID = "epsilon_feooh_thompson_2017_supplement_table_s1"


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
    rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
    return dataset, rows


def test_thompson_hc_material_is_phase_specific_and_diffraction_ready():
    document, source, _ = _document_source_and_executable()

    assert document["phase"] == "hydrogen-centered Pnnm low-spin epsilon-FeOOH"
    assert document["space_group"] == "Pnnm"
    assert document["space_group_number"] == 58
    assert document["formula_units_per_cell"] == 2
    assert source["reference"]["doi"] == "10.1002/2017JB014168"

    counts = {"Fe": 0.0, "O": 0.0, "H": 0.0}
    for site in document["atom_sites"]:
        multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
        counts[site["element"]] += multiplicity * site["occupancy"]
    assert counts == pytest.approx({"Fe": 2.0, "O": 4.0, "H": 2.0})

    lattice = document["lattice"]
    assert lattice["a"] * lattice["b"] * lattice["c"] == pytest.approx(47.975187204)
    assert "Subtracting 0.25 from y" in document["source"]["coordinate_transform"]


def test_thompson_hc_bm3_parameters_reference_state_and_validity():
    _, source, record = _document_source_and_executable()

    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 58.62, "K0": 223.0, "K0_prime": 4.07},
    }
    assert source["parameter_errors"] == {
        "V0": 0.05,
        "K0": 2.0,
        "K0_prime": 0.03,
    }
    assert source["parameter_error_confidence"] is None
    assert source["fixed_parameters"] == []
    assert source["temperature_ref"] == 0.0
    assert source["pressure_range_status"] == "theoretical"
    assert source["experimental_pressure_range_gpa"] == [30.0, 140.0]
    assert source["fit_datasets"] == [DATASET_ID]

    assert record.pressure(58.62) == pytest.approx(0.0, abs=1.0e-12)
    assert record.eos.bulk_modulus(58.62) == pytest.approx(223.0)
    for pressure in (30.0, 60.0, 90.0, 139.99):
        volume = record.volume(pressure, check_validity=True)
        assert record.pressure(volume, check_validity=True) == pytest.approx(
            pressure, rel=1.0e-11
        )


def test_thompson_supplement_table_s1_transcription_and_hc_curve():
    document, source, record = _document_source_and_executable()
    dataset, rows = _dataset_rows(document)

    assert len(rows) == 17
    assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
    assert rows[0] == {
        "pressure_gpa": "0.01",
        "unit_cell_volume_a3": "59.88",
        "feo6_volume_a3": "10.18",
        "o_o_distance_a": "2.44",
        "o_h_distance_a": "1.35",
        "o_h_o_angle_deg": "179.37",
        "lattice_a_a": "4.83",
        "lattice_b_a": "4.24",
        "lattice_c_a": "2.91",
    }
    assert rows[-1]["pressure_gpa"] == "139.99"
    assert rows[-1]["unit_cell_volume_a3"] == "42.40"
    assert dataset["uncertainty"]["status"] == "not_reported"

    hc_rows = [row for row in rows if float(row["pressure_gpa"]) >= 30.0]
    assert len(hc_rows) == 12
    residuals = np.array(
        [
            record.pressure(float(row["unit_cell_volume_a3"]))
            - float(row["pressure_gpa"])
            for row in hc_rows
        ]
    )
    assert np.sqrt(np.mean(residuals**2)) == pytest.approx(
        0.10819370267208199, abs=1.0e-12
    )
    assert np.max(np.abs(residuals)) == pytest.approx(0.1696839117523723, abs=1.0e-12)
    reproduction = source["scientific_validation"]["numerical_reproduction"]
    assert reproduction["published_curve_on_rounded_hc_rows"][
        "pressure_rmse_gpa"
    ] == pytest.approx(0.10819370267208199)


def test_thompson_hc_diagnostic_refit_recovers_published_coefficients():
    document, source, _ = _document_source_and_executable()
    _, rows = _dataset_rows(document)
    hc_rows = [row for row in rows if float(row["pressure_gpa"]) >= 30.0]
    pressure = np.array([float(row["pressure_gpa"]) for row in hc_rows])
    volume = np.array([float(row["unit_cell_volume_a3"]) for row in hc_rows])

    def residuals(parameters):
        v0, k0, k0_prime = parameters
        compression = (v0 / volume) ** (1.0 / 3.0)
        predicted = (
            1.5
            * k0
            * (compression**7 - compression**5)
            * (1.0 + 0.75 * (k0_prime - 4.0) * (compression**2 - 1.0))
        )
        return predicted - pressure

    result = least_squares(
        residuals,
        x0=[58.62, 223.0, 4.07],
        bounds=([50.0, 100.0, 1.0], [70.0, 400.0, 8.0]),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    )
    expected = source["scientific_validation"]["numerical_reproduction"][
        "diagnostic_unweighted_pressure_refit"
    ]
    assert result.x == pytest.approx(
        [expected["V0"], expected["K0"], expected["K0_prime"]],
        abs=1.0e-6,
    )
    assert math.sqrt(np.mean(result.fun**2)) == pytest.approx(
        expected["pressure_rmse_gpa"], abs=1.0e-12
    )


def test_thompson_density_benchmarks_and_hoc_exclusion():
    document, source, record = _document_source_and_executable()
    molar_mass_g_mol = 55.845 + 2.0 * 15.999 + 1.008
    avogadro = 6.02214076e23

    for pressure, published in ((60.0, 6.055), (90.0, 6.428)):
        volume = record.volume(pressure)
        density = 2.0 * molar_mass_g_mol / (avogadro * volume * 1.0e-24)
        assert density == pytest.approx(published, abs=0.0021)

    validation = source["scientific_validation"]
    hoc, hc = validation["reported_parameterizations"]
    assert hoc["branch"] == "hydrogen off-center (HOC)"
    assert hoc["disposition"] == "excluded"
    assert "7.765 GPa" in hoc["reason"]
    assert hc["branch"] == "hydrogen centered (HC)"
    assert hc["disposition"] == "selected"
    assert (
        "not volumetric EOS coefficients"
        in validation["elastic_density_boundary"]["finding"]
    )
    assert all("_hoc_" not in item["identifier"] for item in document["eos_records"])
