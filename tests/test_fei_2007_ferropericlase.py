import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, search_eos_records
from scripts.reproduce_fei_2007_ferropericlase import (
    DOI,
    FP20_LS_REFIT,
    TABLES,
    bm3_pressure,
    normalized_csv,
    refit_fp20_ls,
    reproduce,
    selected_rows,
)
from scripts.validate_primary_eos_refits import _fit_record

ROOT = Path(__file__).resolve().parents[1]
# Figure captions 1-3, checked independently against paragraphs 12 and 14.
PUBLISHED = {
    "mg080fe020o_fei_2007_hs_b1_bm3": (76.16, 158, 4),
    "mg080fe020o_fei_2007_ls_b1_bm3": (74.2, 170, 4),
    "mg061fe039o_fei_2007_hs_b1_bm3": (77.48, 156, 4),
    "mg061fe039o_fei_2007_ls_b1_bm3": (73.6, 170, 4),
    "mg042fe058o_fei_2007_hs_b1_bm3": (79.36, 153, 4),
}


def documents():
    return [
        json.loads((ROOT / f"peritheos/data/materials/{name}.eosmat").read_text())
        for name in TABLES
    ]


def test_fei_source_owned_catalog_scope_and_published_parameters():
    found = search_eos_records(doi=DOI)
    assert {r.identifier for r in found} == set(PUBLISHED) | {FP20_LS_REFIT}
    for doc in documents():
        for record in doc["eos_records"]:
            if (
                record["reference"]["doi"] != DOI
                or record["record_kind"] != "published"
            ):
                continue
            assert (
                tuple(record["eos"]["parameters"].values())
                == PUBLISHED[record["identifier"]]
            )
            assert record["record_kind"] == "published"
            assert record["fixed_parameters"] == ["K0_prime"]
            assert record["volume_basis"]["formula_units"] == 4
            assert record["temperature_ref"] == 300
            assert record["spin_path_context"]["path"] == "compression"
            assert record["parameter_covariance"] is None
    fp39 = documents()[1]
    later = [
        r for r in fp39["eos_records"] if "solomatova_2016_fei_" in r["identifier"]
    ]
    assert len(later) == 2
    assert all(r["reference"]["doi"] == "10.2138/am-2016-5510" for r in later)
    assert [r["eos"]["parameters"]["V0"] for r in later] == [77.49, 74.83]


def test_official_tables_transcribe_every_row_and_preserve_checksums():
    for doc, count in zip(documents(), [46, 77, 27]):
        dataset = doc["datasets"][-1]
        for key in ("resource", "original_resource"):
            resource = dataset[key]
            assert (
                hashlib.sha256(
                    (ROOT / "peritheos/data" / resource["path"]).read_bytes()
                ).hexdigest()
                == resource["sha256"]
            )
        text = normalized_csv(doc["identifier"])
        assert (
            text == (ROOT / "peritheos/data" / dataset["resource"]["path"]).read_text()
        )
        rows = list(csv.DictReader(io.StringIO(text)))
        assert len(rows) == count
        assert [int(row["source_order"]) for row in rows] == list(range(1, count + 1))
        raw = (
            ROOT / "peritheos/data" / dataset["original_resource"]["path"]
        ).read_text()
        fields = [line.split() for line in raw.splitlines()[1:] if line.strip()]
        for row, original in zip(rows, fields):
            assert [
                row[k]
                for k in (
                    "sample_name",
                    "nacl_lattice_a_angstrom",
                    "nacl_lattice_error_angstrom",
                    "pressure_gpa",
                    "sample_lattice_a_angstrom",
                )
            ] == original[:5]
            assert row["sample_lattice_error_angstrom"] == (
                original[5] if len(original) == 6 else ""
            )


def test_paths_phase_boundaries_and_missing_errors_are_not_fabricated():
    fp39 = list(csv.DictReader(io.StringIO(normalized_csv("mg061fe039o"))))
    assert sum(r["compression_path"] == "compression" for r in fp39) == 38
    assert sum(r["compression_path"] == "decompression" for r in fp39) == 39
    assert fp39[0]["pressure_gpa"] == "-0.10"
    assert fp39[8]["sample_lattice_error_angstrom"] == "0.0000"
    fp20 = list(csv.DictReader(io.StringIO(normalized_csv("mg080fe020o"))))
    assert next(r for r in fp20 if r["pressure_gpa"] == "29.30")["nacl_phase"] == "B1"
    assert next(r for r in fp20 if r["pressure_gpa"] == "26.86")["nacl_phase"] == "B2"
    fp58 = list(csv.DictReader(io.StringIO(normalized_csv("mg042fe058o"))))
    assert all(r["sample_lattice_error_angstrom"] == "" for r in fp58)
    assert all(
        r["volume_a3_conventional_cell"] == ""
        for r in fp58
        if float(r["pressure_gpa"]) >= 44
    )
    assert float(fp58[0]["volume_a3_conventional_cell"]) == pytest.approx(4.2214**3)
    doc = documents()[2]
    selected, _ = selected_rows(fp58, doc["eos_records"][0])
    assert len(selected) == 7  # Explicitly assume the rounded endpoint includes 42.75.
    assert selected[-1]["sample_name"] == "mw58_033"
    assert max(float(row["pressure_gpa"]) for row in selected) == 42.75
    assert all(float(row["pressure_gpa"]) < 44 for row in selected)


def test_published_equations_execute_independently_on_correct_cell_basis():
    for doc in documents():
        material = Material.from_eosmat(doc)
        for record in doc["eos_records"]:
            if (
                record["reference"]["doi"] != DOI
                or record["record_kind"] != "published"
            ):
                continue
            v0, k0, kp = PUBLISHED[record["identifier"]]
            volumes = v0 * np.array([1, 0.95, 0.85, 0.75])
            assert material.get_eos_record(record["identifier"]).pressure(
                volumes
            ) == pytest.approx(bm3_pressure(volumes, v0, k0, kp))
    residuals = {r["record"]: r for r in reproduce()}
    assert residuals["mg061fe039o_fei_2007_ls_b1_bm3"][
        "pressure_rmse_gpa"
    ] == pytest.approx(4.0474956584)


def test_independent_refits_never_replace_published_parameters_or_include_decompression():
    for doc in documents():
        original = json.dumps(doc, sort_keys=True)
        for record in doc["eos_records"]:
            if (
                record["reference"]["doi"] != DOI
                or record["record_kind"] != "published"
            ):
                continue
            fit = _fit_record(doc, record, doc["datasets"][-1])
            assert fit["fit_kind"] == "independent_validation_refit"
            assert fit["absolute_sigma"] is False
            assert (
                "Compression and decompression remain separate" in fit["qualification"]
            )
            if record["identifier"] == "mg061fe039o_fei_2007_ls_b1_bm3":
                assert fit["observations"] == 24
                assert next(
                    p["refit"] for p in fit["parameters"] if p["parameter"] == "K0"
                ) == pytest.approx(176.65344, rel=1e-6)
            if record["identifier"] == "mg042fe058o_fei_2007_hs_b1_bm3":
                assert fit["observations"] == 7
                assert fit["status"] == "similar"
                assert fit["parameters"][0]["refit"] == pytest.approx(153.1773703)
                assert fit["parameters"][0]["refit_error"] == pytest.approx(1.4284155)
                assert fit["published_rmse_gpa"] == pytest.approx(0.4878806783)
                assert "rounded-endpoint interpretation" in fit["selection"]
            if record["identifier"] == "mg080fe020o_fei_2007_ls_b1_bm3":
                assert fit["status"] == "parity_not_achieved"
        assert json.dumps(doc, sort_keys=True) == original


def test_fp20_ls_stored_refit_is_distinct_reproducible_and_has_joint_uncertainty():
    doc = documents()[0]
    record = next(r for r in doc["eos_records"] if r["identifier"] == FP20_LS_REFIT)
    source = next(
        r
        for r in doc["eos_records"]
        if r["identifier"] == record["derived_from_record"]
    )
    result = refit_fp20_ls()
    assert record["record_kind"] == "refit"
    assert "Peritheos refit" in record["label"]
    assert "default_for" not in record
    assert source["record_kind"] == "published"
    assert source["eos"]["parameters"] == {"V0": 74.2, "K0": 170, "K0_prime": 4}
    assert source["parameter_covariance"] is None
    assert record["eos"]["parameters"] == pytest.approx(result["parameters"])
    assert record["parameter_errors"] == pytest.approx(result["parameter_errors"])
    assert result["observations"] == 19
    assert result["selected_source_rows"] == list(range(17, 27)) + list(range(38, 47))
    assert (
        record["fit_provenance"]["selection"]["included_sample_names"]
        == result["selected_sample_names"]
    )
    assert result["pressure_range_gpa"] == [40.95, 95.48]
    assert FP20_LS_REFIT in doc["datasets"][-1]["used_by_eos_records"]
    assert (
        record["fit_provenance"]["dataset_resource"] == doc["datasets"][-1]["resource"]
    )

    # Independent analytic BM2 Jacobian checks the covariance scaling and order.
    rows = list(csv.DictReader(io.StringIO(normalized_csv("mg080fe020o"))))
    rows = [r for r in rows if int(r["source_order"]) in result["selected_source_rows"]]
    volume = np.array([float(r["volume_a3_conventional_cell"]) for r in rows])
    pressure = np.array([float(r["pressure_gpa"]) for r in rows])
    v0, k0 = result["parameters"]["V0"], result["parameters"]["K0"]
    ratio = (v0 / volume) ** (1 / 3)
    jacobian = np.column_stack(
        (
            1.5 * k0 / v0 * (7 / 3 * ratio**7 - 5 / 3 * ratio**5),
            1.5 * (ratio**7 - ratio**5),
        )
    )
    residuals = bm3_pressure(volume, v0, k0, 4) - pressure
    covariance = np.linalg.inv(jacobian.T @ jacobian) * np.sum(residuals**2) / 17
    stored = record["parameter_covariance"]
    assert stored["parameter_order"] == ["V0", "K0"]
    assert np.asarray(stored["matrix"]) == pytest.approx(covariance, rel=1e-5)
    assert np.sqrt(np.diag(covariance)) == pytest.approx(
        [0.37404085, 4.1782002], rel=1e-5
    )
    executable = Material.from_eosmat(doc).get_eos_record(FP20_LS_REFIT)
    assert executable.pressure(volume) == pytest.approx(pressure + residuals)
    assert np.sqrt(np.mean(residuals**2)) < 0.87
    assert result["source_curve_pressure_rmse_gpa"] > 1.76
    check = _fit_record(doc, record, doc["datasets"][-1])
    assert check["status"] == "parity"
    assert check["fit_kind"] == "stored_refit_reproduction"
    assert "not Fei's published coefficients" in check["qualification"]
