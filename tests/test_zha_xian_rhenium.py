import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_zha_2004_xian_2022_rhenium import (
    DATA,
    DATASET,
    OUTPUT,
    RECORD,
    ledger_outcome,
    pressure,
    reproduce,
    rows,
)

ROOT = Path(__file__).resolve().parents[1]


def test_reproduction_is_current_and_recovers_source_coefficients(assert_audit_close):
    result = reproduce()
    assert_audit_close(result, json.loads(OUTPUT.read_text()), rel=1e-10, abs=1e-12)
    fit = result["zha_isochore_refit"]
    assert fit["observations"] == 8
    assert fit["parameters"]["alpha_KT_ref"] == pytest.approx(0.00776, abs=0.000005)
    assert fit["parameters"]["dK_dT_V"] == pytest.approx(-0.00815, abs=0.000005)
    assert fit["published_pressure_rmse_gpa"] == pytest.approx(0.40849, abs=0.00001)
    assert all(
        r["max_abs_difference_gpa"] < 0.1
        for r in result["zha_table4_free_branch_checks"]
    )
    assert result["zha_table5"]["states"] == 147
    assert result["zha_table5"]["fixed_V0_branch_max_abs_residual_gpa"] < 0.02
    assert result["zha_table5"]["continuous_eq6_max_abs_difference_gpa"] > 0.7


def test_zha_public_evaluation_and_inversion_follow_equation6():
    record = get_eos_record(RECORD)
    v = np.array([29.482, 29.4087, 29.4087 / 1.2])
    t = np.array([1480.2, 3000.0, 3000.0])
    p = record.pressure(v, t)
    assert p == pytest.approx(pressure(v, t), rel=1e-12, abs=1e-10)
    assert record.volume(p, t) == pytest.approx(v, rel=1e-10)
    # Table IV's independently printed free-V0 3000 K fit gives 116.053 GPa.
    assert p[-1] == pytest.approx(116.053, abs=0.1)
    assert record.pressure(29.4087, 300) == pytest.approx(0, abs=1e-12)


def test_gold_discrepancy_is_preserved_and_not_explained_by_rounding():
    diagnostic = reproduce()["zha_gold_calibration_diagnostic"]
    observations = rows("table2-paired-pvt")
    volume = np.array([float(r["au_volume_a3"]) for r in observations])
    temperature = np.array([float(r["temperature_k"]) for r in observations])
    # Validate the actual public record against an independent source formula,
    # applying Zha's reference volume through equivalent compression ratios.
    au = get_eos_record("gold_anderson_1989_bm3_1")
    calculated = au.pressure(volume * 67.79 / 67.847, temperature)
    assert calculated == pytest.approx(diagnostic["calculated_pressure_gpa"], abs=1e-10)
    assert calculated[-1] == pytest.approx(6.97774, abs=0.00001)
    assert diagnostic["max_abs_volume_minus_printed_mean_a_cubed_a3"] < 0.00005
    assert diagnostic["max_abs_pressure_change_using_mean_of_four_a_gpa"] < 0.02
    assert diagnostic["rounding_corner_residual_interval_gpa"][-1][0] > 0.29
    alternatives = diagnostic["anderson_section3_exploratory_pairs"]
    assert alternatives[0]["max_abs_residual_gpa"] < 0.06
    assert all(r["max_abs_residual_gpa"] > 0.005 for r in alternatives)
    # No closer exploratory pair is silently substituted into the Re calibration.
    re_doc = get_material_document("rhenium")
    record = next(r for r in re_doc["eos_records"] if r["identifier"] == RECORD)
    calibration = record["pressure_calibration"]
    assert calibration["status"] == "partially_resolved"
    method = calibration["methods"][0]
    assert method["related_eos_record"] == "gold_anderson_1989_bm3_1"
    assert "reference_eos_record" not in method


def test_zha_primary_rows_markers_and_derived_branches_are_preserved():
    doc = get_material_document("rhenium")
    record = next(r for r in doc["eos_records"] if r["identifier"] == RECORD)
    assert record["fit_datasets"] == [DATASET]
    assert record["experimental_pressure_range_gpa"] == [6.41, 8.47]
    assert record["experimental_temperature_range_k"] == [1380.3, 1914.5]
    assert ledger_outcome(record)["status"] == "similar"
    selected = [
        d for d in doc["datasets"] if d["identifier"].startswith("rhenium_zha_2004_")
    ]
    assert sorted(d["row_count"] for d in selected) == [6, 8, 12, 147]
    for dataset in selected:
        path = DATA.parent / dataset["resource"]["path"]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )
        assert len(path.read_text().splitlines()) - 1 == dataset["row_count"]
    primary = rows("table2-paired-pvt")
    assert float(primary[1]["volume_a3"]) == 29.482
    assert float(primary[1]["au_volume_a3"]) == 67.8225
    assert float(primary[6]["au_a_error_angstrom"]) == 0.0038
    assert float(primary[-1]["temperature_k"]) == 1914.5
    material = Material.from_eosmat(doc)
    dataset = material.get_dataset(DATASET)
    assert len(dataset["temperature_k"]) == 8


def test_xian_unresolved_equations_are_not_promoted_or_sign_corrected():
    source = json.loads(
        (ROOT / "docs/data/xian-2022-rhenium-source-audit.json").read_text()
    )
    assert source["outcome"] == "deferred_incomplete_model"
    assert source["printed_surrogate"]["specific_volume_cm3_g_coefficients"] == [
        0.04789,
        -3.635e-7,
        -1.651e-11,
    ]
    checks = reproduce()["xian_literal_equations_37_40"]
    assert all(r["zero_pressure_expansivity_per_k"] < 0 for r in checks)
    assert checks[-1]["pressure_at_printed_cold_V0_gpa"] == pytest.approx(
        -7.52618, abs=0.00001
    )
    assert all(
        r["reference"].get("doi") != source["doi"]
        for r in get_material_document("rhenium")["eos_records"]
    )
