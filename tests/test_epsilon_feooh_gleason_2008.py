import csv
import hashlib
import io
import json
from importlib import resources
from pathlib import Path

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM2
from peritheos.materials import Material
from scripts.reproduce_gleason_2008_epsilon_feooh import (
    fit_fixed_v0,
    g_vs_g_volume,
    reduce_volume_to_reference,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIAL_ID = "e_feooh"
RECORD_ID = "e_feooh_gleason_2008_bm2_1"
DATASET_ID = "epsilon_feooh_gleason_2008_deposit_table2_pvt"


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


def _arrays(rows):
    return (
        np.array([float(row["pressure_gpa"]) for row in rows]),
        np.array([float(row["cell_volume_a3"]) for row in rows]),
        np.array([float(row["temperature_celsius"]) + 273.15 for row in rows]),
    )


def test_gleason_record_encodes_the_staged_source_protocol():
    document, source, executable = _document_source_and_executable()

    assert source["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 66.3, "K0": 158.0},
    }
    assert source["parameter_errors"] == {"V0": 0.5, "K0": 5.0}
    assert source["fixed_parameters"] == ["V0"]
    assert source["temperature_ref"] == 300.0
    assert source["experimental_temperature_range_k"] == [473.15, 673.15]
    assert source["fit_datasets"] == [DATASET_ID]
    assert source["fit_protocol"]["free_parameters"] == ["K0"]
    assert source["fit_protocol"]["fixed_parameters"] == {
        "V0": 66.3,
        "K0_prime": 4.0,
        "Tr": 300.0,
        "alpha0": 2.3e-5,
    }
    assert source["fit_protocol"]["exclusions"] == []
    assert source["fit_protocol"]["weights"] is None
    assert source["fit_protocol"]["temperature_handling"]["status"] == (
        "numerically_recovered_not_explicitly_printed"
    )
    assert (
        "not stored as an epsilon thermal EOS coefficient"
        in source["fit_protocol"]["temperature_handling"]["source_scope"]
    )

    assert document["phase"] == "hydrogen-off-center P21nm epsilon-FeOOH"
    assert source["scientific_validation"]["phase_boundary"]["finding"].startswith(
        "Only rows refined as the source's high-pressure P21mn epsilon phase"
    )
    assert isinstance(executable.eos, BM2)
    assert executable.reference_temperature == 300.0
    assert executable.eos.pressure(66.3) == pytest.approx(0.0, abs=1.0e-12)


def test_gleason_table2_is_complete_and_every_observation_is_hot():
    document, _, _ = _document_source_and_executable()
    dataset, rows = _dataset_rows(document)

    assert len(rows) == 49
    assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
    assert rows[0]["pressure_gpa"] == "19.55"
    assert rows[0]["temperature_celsius"] == "200"
    assert rows[0]["cell_volume_a3"] == "60.50"
    assert rows[-1]["pressure_gpa"] == "4.63"
    assert rows[-1]["temperature_celsius"] == "310"
    assert rows[-1]["cell_volume_a3"] == "64.6"
    assert min(float(row["temperature_celsius"]) for row in rows) == 200.0
    assert max(float(row["temperature_celsius"]) for row in rows) == 400.0
    assert all(
        float(row["temperature_celsius"]) != pytest.approx(26.85) for row in rows
    )
    assert dataset["uncertainty"]["confidence"] == "95_percent_confidence_limit"
    assert dataset["source_artifact"]["archive_sha256"] == (
        "11cb79a0961ae7664650265d72aadc6f1997caf8c6b2984950c6738aca3aef07"
    )


def test_gleason_temperature_reduction_recovers_published_bulk_modulus():
    document, source, _ = _document_source_and_executable()
    _, rows = _dataset_rows(document)
    pressure, volume, temperature = _arrays(rows)

    corrected_volume = reduce_volume_to_reference(volume, temperature)
    corrected_fit = fit_fixed_v0(corrected_volume, pressure)
    raw_fit = fit_fixed_v0(volume, pressure)
    published = BM2(V0=66.3, K0=158.0)
    corrected_published_residual = published.pressure(corrected_volume) - pressure
    raw_published_residual = published.pressure(volume) - pressure

    assert corrected_fit.parameters["K0"] == pytest.approx(158.0986983485346)
    assert corrected_fit.standard_errors["K0"] == pytest.approx(1.6723069894662843)
    assert np.sqrt(np.mean(corrected_published_residual**2)) == pytest.approx(
        0.9457320622338763
    )
    assert raw_fit.parameters["K0"] == pytest.approx(175.62426318239244)
    assert np.sqrt(np.mean(raw_published_residual**2)) == pytest.approx(
        1.4987955733085832
    )
    assert g_vs_g_volume(volume, pressure) == pytest.approx(66.16275734712724)

    reproduction = source["scientific_validation"]["numerical_reproduction"]
    assert reproduction["source_faithful_temperature_reduced_bm2"][
        "K0_refit_gpa"
    ] == pytest.approx(corrected_fit.parameters["K0"])
    assert reproduction["controls"][
        "raw_hot_volumes_fixed_V0_unweighted_K0_gpa"
    ] == pytest.approx(raw_fit.parameters["K0"])


def test_gleason_pressure_calibration_and_generated_refit_ledger():
    _, source, _ = _document_source_and_executable()
    calibration = source["pressure_calibration"]

    assert calibration["status"] == "partially_resolved"
    assert calibration["methods"][0]["reference"]["authors"] == [
        "Shim",
        "Duffy",
        "Takemura",
    ]
    assert calibration["recalculation"]["status"] == ("missing_calibrant_observations")

    ledger = json.loads(
        (ROOT / "docs" / "data" / "primary-eos-refits.json").read_text(encoding="utf-8")
    )
    result = next(
        row for row in ledger["records"] if row["record_identifier"] == RECORD_ID
    )
    assert result["status"] == "parity"
    assert result["observations"] == 49
    assert result["fit_kind"] == ("source_pooled_multitemperature_pv_reduced_to_300k")
    assert (
        result["temperature_reduction"][
            "hot_rows_misclassified_as_reference_temperature"
        ]
        == 0
    )
    assert result["parameters"][0]["refit"] == pytest.approx(158.0986983485346)
    assert result["raw_hot_volume_control"]["K0_refit_gpa"] == pytest.approx(
        175.62426318239244
    )
