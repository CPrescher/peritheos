import csv
import hashlib
import io
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_thermal_eos
from peritheos.materials import Material
from peritheos.units import cell_volume_to_molar_volume

RECORD_ID = "ca_perovskite_shim_2000_dac_only_bm3_mgd_4"
DATASET_ID = "ca_perovskite_shim_2000_jgr_table1_dac_pvt"
TARGET_DOI = "10.1029/2000jb900183"


def _source_record():
    document = get_material_document("ca_perovskite")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    return document, record


def _source_rows(document):
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


def test_shim_2000_jgr_record_is_distinct_nonpreferred_dac_only_fit():
    document, record = _source_record()
    matching = [
        item
        for item in document["eos_records"]
        if item["reference"].get("doi", "").lower() == TARGET_DOI
    ]

    assert [item["identifier"] for item in matching] == [RECORD_ID]
    assert record["record_kind"] == "published"
    assert record["equation_kind"] == "thermal"
    assert record["default"] is False
    assert record["fit_datasets"] == [DATASET_ID]
    assert document["space_group"] == "Pm-3m"
    assert document["formula_units_per_cell"] == 1
    assert record["phase_selection"]["cell_basis"] == (
        "one-formula-unit pseudo-cubic cell (Z=1)"
    )
    assert record["scientific_validation"]["excluded_preferred_combined_fit"] == {
        "parameters": {"gamma0": 1.92, "q": 0.6},
        "reason": (
            "The preferred result combines the 34 Table 1 DAC rows with Wang et "
            "al. (1996) LVP observations that are only plotted in this article and "
            "are not provided numerically in an official supplement. It therefore "
            "fails the complete-checksummed-dataset rule and is not represented by "
            "this record."
        ),
    }

    default_record = next(
        item for item in document["eos_records"] if item.get("default", False)
    )
    assert default_record["identifier"] == "ca_perovskite_shim_2000_bm3_1"
    assert default_record["reference"]["doi"].lower() == (
        "10.1016/s0031-9201(00)00154-0"
    )


def test_shim_2000_jgr_dac_only_parameters_and_fixed_free_protocol():
    document, record = _source_record()

    assert record["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {
            "V0": pytest.approx(45.5817973939),
            "K0": 236.0,
            "K0_prime": 3.9,
        },
    }
    assert record["parameter_errors"] == {
        "V0": pytest.approx(0.0332107813),
        "K0": 4.0,
        "K0_prime": 0.2,
    }
    assert record["fixed_parameters"] == ["V0", "K0", "K0_prime"]
    assert record["thermal"] == {
        "type": "MieGruneisenDebye",
        "model": "mie_gruneisen_debye",
        "debye_temperature_law": "integrated_gruneisen",
        "parameters": {
            "Tr": 300.0,
            "theta0": 1000.0,
            "gamma0": 2.0,
            "q": 0.9,
            "n": 5.0,
        },
        "parameter_errors": {
            "Tr": None,
            "theta0": None,
            "gamma0": 0.08,
            "q": 0.4,
            "n": None,
        },
        "fixed_parameters": ["Tr", "theta0", "n"],
    }
    assert record["experimental_pressure_range_gpa"] == [20.6, 68.5]
    assert record["experimental_temperature_range_k"] == [1272.0, 2380.0]

    executable = Material.from_eosmat(
        document, record_identifiers=[RECORD_ID]
    ).eos_records[0]
    assert executable.reference_volume == pytest.approx(45.5817973939)
    assert executable.pressure(45.5817973939, 300.0) == pytest.approx(0.0, abs=1.0e-10)


def test_shim_2000_jgr_table1_transcription_is_complete_and_checksummed():
    document, _ = _source_record()
    dataset, rows = _source_rows(document)

    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert len(rows) == 34
    assert [row["source_order"] for row in rows] == [str(i) for i in range(1, 35)]
    assert sum(row["pressure_medium"] == "nacl" for row in rows) == 17
    assert sum(row["pressure_medium"] == "argon" for row in rows) == 17
    assert rows[0] == {
        "source_order": "1",
        "experiment_group": "experiment_1",
        "pressure_medium": "nacl",
        "pressure_gpa": "35.5",
        "pressure_standard_deviation_gpa": "0.7",
        "volume_a3": "41.599",
        "volume_standard_deviation_a3": "0.062",
        "temperature_k": "1416",
        "temperature_standard_deviation_k": "54",
    }
    assert rows[16]["pressure_gpa"] == "47.2"
    assert rows[16]["temperature_k"] == "1585"
    assert rows[17]["pressure_gpa"] == "24.3"
    assert rows[17]["temperature_standard_deviation_k"] == "141"
    assert rows[-1]["pressure_gpa"] == "64.3"
    assert rows[-1]["volume_standard_deviation_a3"] == "0.083"
    assert rows[-1]["temperature_k"] == "1752"


def test_shim_2000_jgr_unweighted_dac_only_refit_has_parameter_parity():
    document, record = _source_record()
    _, rows = _source_rows(document)
    values = {
        name: np.array([float(row[name]) for row in rows])
        for name in ("pressure_gpa", "volume_a3", "temperature_k")
    }
    volume_scale = cell_volume_to_molar_volume(1.0, 1.0)
    reference = BM3(V0=2.745, K0=236.0, K0_prime=3.9)

    result = fit_thermal_eos(
        MieGruneisenDebye,
        reference,
        volume=values["volume_a3"] * volume_scale,
        temperature=values["temperature_k"],
        pressure=values["pressure_gpa"],
        initial={"gamma0": 2.0, "q": 0.9},
        fixed={"Tr": 300.0, "theta0": 1000.0, "n": 5.0},
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds={"gamma0": (1.0e-9, 10.0), "q": (1.0e-9, 10.0)},
        max_nfev=5000,
    )

    assert result.success
    assert result.free_parameters == ("gamma0", "q")
    assert [
        result.parameters[name] for name in result.free_parameters
    ] == pytest.approx([2.0006405181, 0.9404932591], rel=1.0e-6)
    assert [
        result.standard_errors[name] for name in result.free_parameters
    ] == pytest.approx([0.0875352616, 0.4223094684], rel=1.0e-6)
    assert abs(result.parameters["gamma0"] - 2.0) < 0.08
    assert abs(result.parameters["q"] - 0.9) < 0.4

    loaded = Material.from_eosmat(document, record_identifiers=[RECORD_ID])
    executable = loaded.eos_records[0]
    published_residuals = (
        np.asarray(executable.pressure(values["volume_a3"], values["temperature_k"]))
        - values["pressure_gpa"]
    )
    refit_residuals = (
        np.asarray(
            result.model.pressure(
                values["volume_a3"] * volume_scale, values["temperature_k"]
            )
        )
        - values["pressure_gpa"]
    )
    assert np.sqrt(np.mean(published_residuals**2)) == pytest.approx(1.05194403)
    assert np.sqrt(np.mean(refit_residuals**2)) == pytest.approx(1.0511631566)
    assert record["scientific_validation"]["independent_refit"]["result"] == ("parity")
