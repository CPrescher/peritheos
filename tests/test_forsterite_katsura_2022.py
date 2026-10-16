"""Forsterite source transcription, final-table mapping and native execution."""

import csv
import hashlib
import json
import re

import jsonschema
import numpy as np
import pytest

from peritheos import Material, eosmat_schema, get_eos_record, get_material_document
from peritheos.errors import MaterialError
from scripts.forsterite_recover_katsura_2022 import (
    COLUMN_MAP,
    DATA,
    PREFIX,
    SOURCES,
    extract,
)
from scripts.forsterite_reproduce_katsura_2022 import (
    RECORD,
    SCRIPT_PARAMETERS,
    TABLE_RECORD,
    author_pressure,
    diagnostics,
    load_rows,
)


def test_forsterite_final_table_mapping_and_unchanged_static_default():
    document = get_material_document("forsterite")
    jsonschema.validate(document, eosmat_schema())
    material = Material.from_eosmat(document)
    assert material.default_record().identifier == "forsterite_finkelstein_2014_bm3_1"
    assert document["eos_records"][0]["eos"]["parameters"] == {
        "V0": 290.1,
        "K0": 130,
        "K0_prime": 4.12,
    }
    record = document["eos_records"][-1]
    assert record["identifier"] == TABLE_RECORD
    assert record["eos"]["parameters"] == dict(V0=289.58, K0=127.4, K0_prime=4.2)
    assert record["parameter_errors"] == dict(V0=None, K0=None, K0_prime=0.4)
    assert record["thermal"]["parameters"] == dict(
        Tr=300,
        theta0=768,
        gamma0=1,
        q=2.4,
        n=7,
    )
    assert record["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert record["thermal"]["thermal_pressure_reference"] == "reference_temperature"
    assert document["space_group"] == "Pbnm"
    assert document["formula_units_per_cell"] == 4
    contents = {}
    for site in document["atom_sites"]:
        contents[site["element"]] = contents.get(site["element"], 0) + (
            int(site["wyckoff"][:-1]) * site["occupancy"]
        )
    assert contents == dict(Mg=8, Si=4, O=16)
    np.testing.assert_allclose(
        np.prod([document["lattice"][axis] for axis in ["a", "b", "c"]]),
        290.1,
        atol=0.05,
    )
    script = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert script["eos"]["parameters"]["K0_prime"] == SCRIPT_PARAMETERS[0]
    assert script["thermal"]["parameters"]["theta0"] == 760
    assert script["thermal"]["parameters"]["gamma0"] == SCRIPT_PARAMETERS[1]
    assert script["thermal"]["parameters"]["q"] == SCRIPT_PARAMETERS[2]
    assert script["parameter_errors"]["K0_prime"] is None
    assert script["thermal"]["parameter_errors"]["gamma0"] is None
    assert document["eos_records"][1]["identifier"] == RECORD
    source = (SOURCES / f"{PREFIX}-author-geotherm.m").read_text(errors="replace")
    for name, expected in {
        "V00_lat_phase": script["eos"]["parameters"]["V0"],
        "KT0_phase": script["eos"]["parameters"]["K0"],
        "KT0P_phase": script["eos"]["parameters"]["K0_prime"],
        "Theta0_phase": script["thermal"]["parameters"]["theta0"],
        "gamma0_phase": script["thermal"]["parameters"]["gamma0"],
        "q_phase": script["thermal"]["parameters"]["q"],
        "n_phase": script["thermal"]["parameters"]["n"],
        "Z_phase": document["formula_units_per_cell"],
    }.items():
        match = re.search(rf"\b{name}\s*=\s*\[([^\]]+)\]", source)
        assert match is not None
        assert float(match.group(1).split()[0]) == expected


def test_forsterite_workbook_recovery_preserves_all_rows_columns_and_formulas():
    text, evidence = extract()
    assert text == (DATA / f"{PREFIX}-olivine-workbook.csv").read_text()
    assert evidence == json.loads((DATA / f"{PREFIX}-olivine-cells.json").read_text())
    rows = list(csv.DictReader(text.splitlines()))
    assert len(rows) == 51
    assert len(rows[0]) == 42
    assert [int(row["source_row"]) for row in rows] == list(range(4, 55))
    assert rows[0]["col_a"] == "M371008"
    assert rows[-1]["col_a"] == "M375135"
    assert float(rows[0]["col_f"]) == 1.03706  # MgO, not sample volume
    assert float(rows[0]["col_an"]) == pytest.approx(1.0300114906130446)
    assert evidence["cells"]["AD4"]["formula"] == "AVERAGE(V4,W4,Z4,AA4)"
    assert evidence["cells"]["AE4"]["formula"] == "STDEV(V4,W4,Z4,AA4)"
    doc = get_material_document("forsterite")
    dataset = next(
        d
        for d in doc["datasets"]
        if d["identifier"] == "forsterite_katsura_2022_olivine_pvt"
    )
    resource = DATA.parent / dataset["resource"]["path"]
    assert (
        hashlib.sha256(resource.read_bytes()).hexdigest()
        == dataset["resource"]["sha256"]
    )
    assert len(load_rows()) == 51
    for row in load_rows():
        index = row["source_row"]
        for name, column in COLUMN_MAP.items():
            assert row[name] == evidence["cells"][f"{column}{index}"]["value"]
    assert min(float(r["pressure_mean_gpa"]) for r in load_rows()) < 0
    assert max(float(r["temperature_corrected_k"]) for r in load_rows()) > 1958


def test_forsterite_source_code_equation_and_observation_checkpoints():
    record = get_eos_record(RECORD)
    rows = load_rows()
    ratio = np.array([float(r["forsterite_volume_ratio"]) for r in rows])
    temperature = np.array([float(r["temperature_corrected_k"]) for r in rows])
    native = record.pressure(289.58 * ratio, temperature)
    # Literal author polynomial/constants versus exact native Debye quadrature.
    np.testing.assert_allclose(
        native, author_pressure(ratio, temperature), atol=7e-5, rtol=0
    )
    result = diagnostics()
    assert result["pressure_rmse_gpa"] == pytest.approx(0.2130082, abs=1e-6)
    assert result["table_comparison"]["pressure_rmse_gpa"] == pytest.approx(
        0.2124994, abs=1e-6
    )
    for point in result["checkpoints"]:
        assert abs(
            point["observed_pressure_gpa"] - point["calculated_pressure_gpa"]
        ) < (3 * point["projected_spread_gpa"])
    # This diagnostic must not turn into a claim that every datum agrees at 3 sigma.
    assert result["max_residual_over_projected_spread"] > 6
    for fit in result["ordinary_refits"]:
        assert fit["success"]
        assert fit["parameters"]["K0_prime"] == pytest.approx(4.2, abs=0.4)
        assert fit["parameters"]["gamma0"] == pytest.approx(1.0, abs=0.02)
        assert fit["parameters"]["q"] == pytest.approx(2.4, abs=0.5)


def test_forsterite_mgo_ratio_identity_uses_retained_unadjusted_readings():
    # Workbook AG is the uncorrected-temperature pressure at the central MgO
    # ratio, independently resolving F as calibrant rather than sample V/V0.
    source = json.loads((DATA / f"{PREFIX}-olivine-cells.json").read_text())["cells"]
    marker = get_eos_record("mgo_b1_tange_2009_vinet")
    for index in [4, 54]:
        ratio = float(source[f"F{index}"]["value"])
        temperature = float(source[f"B{index}"]["value"])
        expected = float(source[f"AG{index}"]["value"])
        actual = marker.pressure(
            marker.reference_volume * ratio, temperature, check_validity=False
        )
        assert actual == pytest.approx(expected, abs=8e-5)


def test_forsterite_reference_derivative_inversion_and_metadata_roundtrip():
    record = get_eos_record(RECORD)
    assert record.pressure(289.58, 300, check_validity=False) == pytest.approx(
        0, abs=1e-12
    )
    assert record.eos.thermal_pressure_increment(
        record.eos.rt_eos.V0, 300
    ) == pytest.approx(0)
    pressure = np.array([[0.0, 5.0, 13.5], [0.0, 5.0, 13.5]])
    temperature = np.array([[310.0], [1900.0]])
    volume = record.volume(pressure, temperature)
    np.testing.assert_allclose(
        record.pressure(volume, temperature), pressure, atol=1e-8
    )
    step = volume * 1e-5
    numerical = (
        -volume
        * (
            record.pressure(volume + step, temperature)
            - record.pressure(volume - step, temperature)
        )
        / (2 * step)
    )
    np.testing.assert_allclose(
        numerical,
        record.eos.bulk_modulus(volume * record.volume_scale, temperature),
        rtol=1e-7,
    )
    document = get_material_document("forsterite")
    exported = Material.from_eosmat(document).to_eosmat()
    jsonschema.validate(exported, eosmat_schema())
    restored = Material.from_eosmat(exported).get_eos_record(RECORD)
    np.testing.assert_allclose(
        restored.pressure(volume, temperature), pressure, atol=1e-8
    )
    assert (
        exported["eos_records"][-1]["parameter_provenance"]
        == document["eos_records"][-1]["parameter_provenance"]
    )
    with pytest.raises(MaterialError):
        record.volume(20, 1900, check_validity=True)


def test_forsterite_candidate_fitting_input_is_source_mapped_and_keeps_missing_sigma_t():
    with (DATA / f"{PREFIX}-candidate-fitting-input.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    source = load_rows()
    assert len(rows) == len(source) == 51
    mapping = {
        "pressure_gpa": "pressure_mean_gpa",
        "pressure_sigma_gpa": "pressure_spread_gpa",
        "volume_ratio": "forsterite_volume_ratio",
        "volume_ratio_sigma": "forsterite_volume_ratio_sigma",
        "temperature_k": "temperature_corrected_k",
    }
    for row, original in zip(rows, source):
        for target, origin in mapping.items():
            assert row[target] == original[origin]
        assert row["temperature_sigma_k"] == ""
    metadata = json.loads((DATA / f"{PREFIX}-candidate-fitting-input.json").read_text())
    assert metadata["missing"] == ["temperature_sigma_k"]
    assert "not_original_author_export" in metadata["status"]


def test_forsterite_2009_complete_primary_table_and_workbook_join():
    with (DATA / "forsterite-katsura-2009-table1.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 52
    assert sum(r["run_label"] == "M371" for r in rows) == 11
    assert sum(r["run_label"] == "M375" for r in rows) == 41
    assert sum(bool(r["workbook_run_id"]) for r in rows) == 50
    assert "M375006" not in {r["workbook_run_id"] for r in rows}
    assert [r["table_row"] for r in rows if not r["workbook_run_id"]] == ["1", "12"]
    assert rows[0]["volume_a3_source_token"] == "290.09(10)"
    assert rows[-1]["volume_a3_source_token"] == "274.62(12)"
    assert rows[-1]["volume_a3_reported_uncertainty"] == "0.12"
    assert rows[0]["pressure_gpa_reported_uncertainty"] == "0.00"
    for row in rows:
        product = np.prod([float(row[f"{axis}_angstrom"]) for axis in "abc"])
        assert product == pytest.approx(float(row["volume_a3"]), abs=0.015)
    dataset = next(
        d
        for d in get_material_document("forsterite")["datasets"]
        if d["identifier"] == "forsterite_katsura_2009_table1_original"
    )
    assert dataset["used_by_eos_records"] == []
    assert (
        hashlib.sha256(
            (DATA.parent / dataset["resource"]["path"]).read_bytes()
        ).hexdigest()
        == dataset["resource"]["sha256"]
    )


def test_forsterite_actual_candidate_input_refit_and_author_optimizer_agree():
    result = diagnostics()["reconstructed_input_refits"]
    assert result["missing_temperature_sigma_count"] == 51
    assert result["bootstrap_performed"] is False
    reproduction = result["parameter_reproduction"]
    assert reproduction["status"] == "parity"
    assert reproduction["parity_basis"] == "within_reported_parameter_errors"
    assert reproduction["uncertainty_reference_record"] == TABLE_RECORD
    assert all(
        p["within_reported_error"] and p["similar"] for p in reproduction["parameters"]
    )
    script = next(
        r
        for r in get_material_document("forsterite")["eos_records"]
        if r["identifier"] == RECORD
    )
    assert script["scientific_validation"]["reproduction_status"] == "parity"
    assert (
        script["scientific_validation"]["numerical_reproduction"]["status"] == "parity"
    )
    all_rows, matched_rows = result["selections"]
    assert all_rows["row_count"] == 51 and all_rows["excluded_run_ids"] == []
    assert matched_rows["row_count"] == 50
    assert matched_rows["excluded_run_ids"] == ["M375006"]
    fitted = all_rows["fits"][0]
    assert fitted["parameters"] == dict(
        K0_prime=pytest.approx(4.249221, abs=1e-5),
        gamma0=pytest.approx(0.99927258, abs=1e-7),
        q=pytest.approx(2.4677172, abs=1e-5),
    )
    for selection in [all_rows, matched_rows]:
        independent, native, author = selection["fits"]
        for name in ["K0_prime", "gamma0", "q"]:
            assert independent["parameters"][name] == pytest.approx(
                author["parameters"][name], abs=1e-5
            )
            assert independent["parameters"][name] == pytest.approx(
                native["parameters"][name], abs=1e-4
            )
        assert (
            independent["pressure_rmse_gpa"]
            <= independent["author_script_pressure_rmse_gpa"]
        )
    assert matched_rows["fits"][0]["pressure_rmse_gpa"] == pytest.approx(
        0.162259167, abs=1e-7
    )


def test_forsterite_global_ledgers_preserve_parity_and_script_recommendation(
    monkeypatch,
):
    from scripts import validate_primary_eos_refits as validator
    from scripts.forsterite_reproduce_katsura_2022 import ROOT
    from scripts.generate_paper_investigation_ledger import build_papers

    document = get_material_document("forsterite")
    document = {
        **document,
        "eos_records": [
            r
            for r in document["eos_records"]
            if r["identifier"] in {RECORD, TABLE_RECORD}
        ],
    }
    monkeypatch.setattr(validator, "list_material_documents", lambda: ["forsterite"])
    monkeypatch.setattr(validator, "get_material_document", lambda _: document)
    current = validator.validate_all()
    saved = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    saved_by_id = {r["record_identifier"]: r for r in saved["records"]}
    for result in current["records"]:
        assert result == saved_by_id[result["record_identifier"]]
        assert result["status"] == "parity"
        assert result["recommended_thermal_record"] == RECORD
        assert result["original_author_workflow_parity"] is False
        assert result["recommended"] == (result["record_identifier"] == RECORD)
    script = saved_by_id[RECORD]
    assert all(p["published_error"] is None for p in script["parameters"])
    papers, _ = build_papers()
    paper = next(p for p in papers if p.get("doi") == "10.1029/2021jb023562")
    assert paper["outcome"] == "reproduced"
    assert paper["statuses"] == {"parity": 2}
    original = next(p for p in papers if p.get("doi") == "10.1016/j.pepi.2008.08.002")
    assert original["outcome"] == "direct_refit_unavailable"
