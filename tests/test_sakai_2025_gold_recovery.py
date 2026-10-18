"""Au observations, pressure ancestry, source precision and replay limits."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import (
    Material,
    get_material_document,
    get_pressure_calibration,
    validate_eosmat_document,
)
from scripts.recover_sakai_2025_gold import paired_rows, write_csv
from scripts.reproduce_sakai_2025_gold import (
    DATA,
    read_rows,
    reproduce,
    ruby2020_from_mao1986,
)

ROOT = Path(__file__).resolve().parents[1]
RECORD = "gold_sakai_2025_rydberg_stacey_1"


def test_sakai_gold_observation_datasets_load_without_changing_coefficients():
    document = get_material_document("gold")
    validate_eosmat_document(document)
    record = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert record["eos"]["parameters"] == {
        "V0": 67.716,
        "K0": 167,
        "K0_prime": 5.78,
        "K_infinity_prime": 5 / 3,
    }
    assert record["parameter_errors"]["K0_prime"] == 0.004
    check = record["scientific_validation"]["primary_data_check"]
    assert check["status"] == "bundled"
    assessment = check["reconstruction_assessment"]
    assert assessment["status"] == "conditional_numerical_parity_achieved"
    assert assessment["conditional_numerical_parity_achieved"] is True
    assert assessment["uncertainty_parity_achieved"] is False
    assert assessment["full_sequential_calibration_reproduced"] is False
    assert assessment["author_fit_parity"] is False
    material = Material.from_eosmat(document, record_identifiers=[RECORD])
    expected = {
        "gold_sakai_2025_calibration_paired_volumes": 326,
        "gold_sakai_2025_table_s4_au_w_volumes": 19,
        "gold_sakai_2025_dewaele_2018_au_re_input": 84,
        "gold_sakai_2025_table_s8_lineage": 35,
        "gold_sakai_2025_table_s14_cu_anchor": 36,
        "gold_sakai_2025_dewaele_2004_au_input": 37,
        "gold_sakai_2025_takemura_2008_au_input": 86,
        "gold_sakai_2025_hirose_2008_au_mgo_input": 21,
    }
    for identifier, count in expected.items():
        dataset = material.get_dataset(identifier)
        assert len(dataset) == count
        assert dataset.checksum_verified
        assert identifier in check["dataset_identifiers"]
    au = material.get_dataset("gold_sakai_2025_table_s4_au_w_volumes")
    assert np.all(np.isnan(au["volume_1_error_a3"]))
    assert np.all(np.isfinite(au["volume_2_error_a3"]))
    assert au["phase_1"].tolist() == ["fcc"] * 19
    assert au["phase_2"].tolist() == ["bcc"] * 19


def test_source_hashes_and_every_normalized_workbook_row_round_trip(tmp_path):
    manifest = json.loads(
        (DATA / "sakai-2025-au-source.json").read_text(encoding="utf-8")
    )
    for item in manifest["products"].values():
        path = ROOT / "peritheos/data" / item["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"]
    raw = json.loads(
        (DATA / "sakai-2025-au-workbook-cells.json").read_text(encoding="utf-8")
    )
    assert raw["workbook_sha256"] == (
        "c559b29fea9c793831e4bc6deef4d62712b5e78821fc420969f709cbd1060ff9"
    )
    rows = paired_rows(raw["sheets"])
    reconstructed = tmp_path / "paired.csv"
    write_csv(reconstructed, rows)
    assert (
        reconstructed.read_bytes()
        == (DATA / "sakai-2025-au-calibration-paired-volumes.csv").read_bytes()
    )
    assert len({(r["source_sheet"], r["source_row"]) for r in rows}) == 326
    # Independent source landmarks: original run order, literal OOXML tokens
    # and blank errors survive; none becomes a fabricated experimental sigma.
    assert rows[0]["volume_1_a3_source_token"] == "46.761648999999998"
    assert {r["run"] for r in rows if r["source_table"] == "S1"} == {"1a", "1b"}
    assert {r["run"] for r in rows if r["source_table"] == "S7"} == {"7a", "7b", "7c"}
    assert len([r for r in rows if r["source_table"] == "S4"]) == 19
    # Context tables stay archived without being passed off as Au ancestors.
    assert "Table S5, run#5" in raw["sheets"]
    assert "Table S6, run#6" in raw["sheets"]
    assert not {"S5", "S6"} & {r["source_table"] for r in rows}


def test_dewaele_printed_rows_and_conventional_volumes_preserve_both_peaks():
    rows = read_rows("gold-dewaele-2018-tables-4-5-au-re.csv")
    assert [sum(r["run"] == str(run) for r in rows) for run in (1, 2, 3)] == [
        36,
        31,
        17,
    ]
    assert rows[35]["au_a111_angstrom"] == "3.3019"
    assert rows[66]["au_a200_angstrom"] == "3.4422"
    assert (
        json.loads(rows[0]["source_tokens_json"])["reported_re_pressure_gpa"] == "10.0"
    )
    assert (
        json.loads(rows[36]["source_tokens_json"])["reported_re_pressure_gpa"] == "9.30"
    )
    for row in rows:
        assert float(row["au_volume_a3"]) == pytest.approx(
            float(row["au_a111_angstrom"]) ** 3
        )
        assert float(row["re_volume_a3"]) == pytest.approx(
            np.sqrt(3)
            * float(row["re_a_angstrom"]) ** 2
            * float(row["re_c_angstrom"])
            / 2
        )
        assert row["volume_error_a3"] == ""
    assert rows[75]["al_fcc_a111_angstrom"] == "3.2"
    assert rows[75]["al_hcp_a_angstrom"] == "2.265"
    assert rows[75]["al_hcp_c_angstrom"] == "3.671"


def test_ruby_coordinate_reconstruction_uses_reported_source_scale():
    mao = get_pressure_calibration("ruby_mao_1986")
    p = np.array(
        [
            float(r["ruby_pressure_classical_gpa"])
            for r in read_rows("gold-dewaele-2004-table1-compression.csv")
        ]
    )
    shift = mao.wavelength_ratio(p) - 1
    assert ruby2020_from_mao1986(p) == pytest.approx(
        1870 * shift * (1 + 5.63 * shift), abs=1e-10
    )
    # The two ruby-pressure columns are retained in the existing resource.
    assert set(read_rows("gold-takemura-2008-table3-compression.csv")[0]) == {
        "run",
        "ruby_pressure_zha_gpa",
        "ruby_pressure_dorogokupets_gpa",
        "lattice_111_angstrom",
    }


def test_reconstructed_pressure_states_and_objective_sensitivity_are_explicit():
    result = reproduce()
    saved = json.loads(
        (ROOT / "docs/data/sakai-2025-gold-reconstruction.json").read_text(
            encoding="utf-8"
        )
    )
    assert result["source_recovery"] == saved["source_recovery"]
    reductions = result["pressure_reductions"]
    assert len(reductions) == 235
    assert sum(r["included_in_diagnostic"] for r in reductions) == 210
    assert {r["source"] for r in reductions} == {
        "Sakai 2025 Table S4",
        "Hirose 2008 Table 1",
        "Dewaele 2018 Tables 4-5",
        "Dewaele 2004 Table I",
        "Takemura-Dewaele 2008 Table III",
    }
    excluded = [r for r in reductions if not r["included_in_diagnostic"]]
    assert len(excluded) == 25
    assert all(
        r["source"] == "Dewaele 2018 Tables 4-5" and r["calculated_pressure_gpa"] < 200
        for r in excluded
    )
    assert result["native_vs_independent_max_difference_gpa"] < 1e-10
    pressure_fit, volume_fit = result["conditional_au_fits"]
    assert pressure_fit["k0_prime"] == pytest.approx(5.7800965844, abs=1e-8)
    assert volume_fit["k0_prime"] == pytest.approx(5.7877715421, abs=1e-8)
    assert round(pressure_fit["k0_prime"], 3) == 5.780
    assert pressure_fit["conditional_numerical_parity_achieved"] is True
    assert volume_fit["conditional_numerical_parity_achieved"] is False
    assert result["conditional_numerical_parity_achieved"] is True
    assessment = result["reconstruction_assessment"]
    assert assessment == saved["reconstruction_assessment"]
    document = get_material_document("gold")
    record = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert (
        assessment
        == record["scientific_validation"]["primary_data_check"][
            "reconstruction_assessment"
        ]
    )
    manifest = json.loads(
        (DATA / "sakai-2025-au-source.json").read_text(encoding="utf-8")
    )
    assert manifest["reconstruction_assessment"] == assessment
    assert (
        manifest["reconstruction_report_sha256"]
        == hashlib.sha256(
            (ROOT / manifest["reconstruction_report"]).read_bytes()
        ).hexdigest()
    )
    for path, sha256 in manifest["recovery_code"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha256
    assert volume_fit["k0_prime"] - pressure_fit["k0_prime"] > 0.004
    assert all(f["optimizer_success"] for f in result["conditional_au_fits"])
    for key in (
        "author_fit_parity",
        "original_weights_recovered",
        "original_covariance_recovered",
        "full_sequential_calibration_reproduced",
        "uncertainty_parity_achieved",
    ):
        assert result[key] is False
    assert result["cu_anchor_conversion_diagnostic"]["rows"] == 35
    for live, archived in zip(
        result["conditional_au_fits"], saved["conditional_au_fits"]
    ):
        assert live["k0_prime"] == pytest.approx(archived["k0_prime"], abs=1e-10)
    # The first branch truly recalibrates W before reducing the Au states.
    branch = result["measurement_only_sequential_branch"][0]
    assert branch["w_fit"]["rows"] == 54
    assert branch["au_fit"]["rows"] == 19
    assert branch["au_fit"]["k0_prime"] == pytest.approx(5.7155171823, abs=1e-8)
