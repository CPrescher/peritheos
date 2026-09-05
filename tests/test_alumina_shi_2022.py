import csv
import hashlib
import io
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_joint_eos
from peritheos.materials import Material

CORUNDUM_RECORD = "alumina_shi_2022_bm3_mgd_2"
CORUNDUM_DATASET = "alumina_shi_2022_table_s1_pvt"
RH2O3_RECORD = "alumina_rh2o3_ii_shi_2022_bm3_mgd_1"


def _source_record():
    document = get_material_document("alumina")
    record = next(
        item
        for item in document["eos_records"]
        if item["identifier"] == CORUNDUM_RECORD
    )
    return document, record


def _source_rows(document):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == CORUNDUM_DATASET
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_bytes()
    )
    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    return dataset, list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))


def test_shi_2022_records_are_distinct_phase_specific_models():
    corundum_document, corundum = _source_record()
    rh2o3_document = get_material_document("alumina_rh2o3_ii")
    rh2o3 = next(
        item
        for item in rh2o3_document["eos_records"]
        if item["identifier"] == RH2O3_RECORD
    )

    assert corundum_document["phase"] == "corundum, trigonal R-3c (hexagonal setting)"
    assert corundum_document["space_group"] == "R-3c"
    assert corundum_document["formula_units_per_cell"] == 6
    assert corundum_document["cell_contents"].startswith("Al12O18")
    assert corundum["eos"]["parameters"]["V0"] == 255.1
    assert rh2o3_document["space_group"] == "Pbcn"
    assert rh2o3_document["formula_units_per_cell"] == 4
    assert rh2o3["eos"]["parameters"]["V0"] == 165.2

    doi_records = []
    for document in (corundum_document, rh2o3_document):
        doi_records.extend(
            item
            for item in document["eos_records"]
            if item["reference"].get("doi", "").lower() == "10.1029/2021jb023805"
        )
    assert [item["identifier"] for item in doi_records] == [
        CORUNDUM_RECORD,
        RH2O3_RECORD,
    ]


def test_shi_2022_selected_corundum_parameters_and_protocol():
    document, record = _source_record()

    assert record["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 255.1, "K0": 246.0, "K0_prime": 4.0},
    }
    assert record["parameter_errors"] == {
        "V0": None,
        "K0": 2.0,
        "K0_prime": None,
    }
    assert record["fixed_parameters"] == ["V0", "K0_prime"]
    assert record["thermal"] == {
        "type": "MieGruneisenDebye",
        "model": "mie_gruneisen_debye",
        "debye_temperature_law": "integrated_gruneisen",
        "parameters": {
            "Tr": 300.0,
            "theta0": 1100.0,
            "gamma0": 1.32,
            "q": 0.8,
            "n": 5.0,
        },
        "parameter_errors": {
            "Tr": None,
            "theta0": None,
            "gamma0": 0.07,
            "q": 0.4,
            "n": None,
        },
        "fixed_parameters": ["Tr", "theta0", "n"],
    }
    assert record["fit_datasets"] == [CORUNDUM_DATASET]
    assert record["experimental_pressure_range_gpa"] == [0.0, 97.7]
    assert record["experimental_temperature_range_k"] == [300.0, 2980.0]

    loaded = Material.from_eosmat(document, record_identifiers=[CORUNDUM_RECORD])
    executable = loaded.eos_records[0]
    assert executable.reference_volume == pytest.approx(255.1)
    assert executable.pressure(255.1, 300.0) == pytest.approx(0.0, abs=1.0e-12)
    assert executable.volume_scale == pytest.approx(0.010036901266666666)


def test_shi_2022_corundum_table_s1_transcription_is_complete():
    document, _ = _source_record()
    dataset, rows = _source_rows(document)

    assert dataset["used_by_eos_records"] == [CORUNDUM_RECORD]
    assert len(rows) == 75
    assert [row["source_order"] for row in rows] == [str(i) for i in range(1, 76)]
    assert sum(row["pressure_calibrant"] == "pt" for row in rows) == 59
    assert sum(row["pressure_calibrant"] == "nacl" for row in rows) == 16
    assert rows[0] == {
        "source_order": "1",
        "pressure_gpa": "0",
        "pressure_uncertainty_gpa": "",
        "temperature_k": "300",
        "temperature_uncertainty_k": "",
        "volume_a3": "255.1",
        "volume_uncertainty_a3": "",
        "pressure_calibrant": "pt",
    }
    assert rows[40] == {
        "source_order": "41",
        "pressure_gpa": "64.4",
        "pressure_uncertainty_gpa": "1.1",
        "temperature_k": "2650",
        "temperature_uncertainty_k": "270",
        "volume_a3": "218.8",
        "volume_uncertainty_a3": "0.2",
        "pressure_calibrant": "pt",
    }
    assert rows[59]["pressure_gpa"] == "69.6"
    assert rows[59]["pressure_calibrant"] == "nacl"
    assert rows[-1]["temperature_k"] == "2730"
    assert rows[-1]["volume_a3"] == "209.3"


def test_shi_2022_unweighted_corundum_refit_has_parameter_parity():
    document, record = _source_record()
    _, rows = _source_rows(document)
    values = {
        name: np.array([float(row[name]) for row in rows])
        for name in ("pressure_gpa", "temperature_k", "volume_a3")
    }
    loaded = Material.from_eosmat(document, record_identifiers=[CORUNDUM_RECORD])
    executable = loaded.eos_records[0]
    scale = executable.volume_scale

    result = fit_joint_eos(
        MieGruneisenDebye,
        BM3,
        volume=values["volume_a3"] * scale,
        temperature=values["temperature_k"],
        pressure=values["pressure_gpa"],
        initial={"rt_eos.K0": 246.0, "gamma0": 1.32, "q": 0.8},
        fixed={
            "rt_eos.V0": 255.1 * scale,
            "rt_eos.K0_prime": 4.0,
            "Tr": 300.0,
            "theta0": 1100.0,
            "n": 5.0,
        },
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds={
            "rt_eos.K0": (12.3, 1230.0),
            "gamma0": (1.0e-9, 13.2),
            "q": (1.0e-9, 10.0),
        },
        max_nfev=5000,
    )

    assert result.success
    assert result.free_parameters == ("rt_eos.K0", "gamma0", "q")
    assert [
        result.parameters[name] for name in result.free_parameters
    ] == pytest.approx([246.3082853542, 1.35931384348, 0.800867575499], rel=2.0e-6)
    assert [
        result.standard_errors[name] for name in result.free_parameters
    ] == pytest.approx([1.33306862429, 0.0672542431192, 0.352700002079], rel=2.0e-6)
    assert abs(result.parameters["rt_eos.K0"] - 246.0) < 2.0
    assert abs(result.parameters["gamma0"] - 1.32) < 0.07
    assert abs(result.parameters["q"] - 0.8) < 0.4

    published_residuals = (
        np.asarray(executable.pressure(values["volume_a3"], values["temperature_k"]))
        - values["pressure_gpa"]
    )
    refit_residuals = (
        np.asarray(
            result.model.pressure(values["volume_a3"] * scale, values["temperature_k"])
        )
        - values["pressure_gpa"]
    )
    assert np.sqrt(np.mean(published_residuals**2)) == pytest.approx(1.23304945511)
    assert np.sqrt(np.mean(refit_residuals**2)) == pytest.approx(1.19138719943)
