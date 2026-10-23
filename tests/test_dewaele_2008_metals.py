import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import (
    get_eos_record,
    get_material_document,
    resolve_dataset_observations,
    search_eos_records,
)
from scripts.reproduce_dewaele_2008_metals import (
    experiment_rows,
    reproduce,
    table4_specs,
    table6_specs,
    vinet,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
RESEARCH = ROOT / "peritheos/data/research/dewaele-2008-paw"


def test_source_table_rows_tokens_and_mixed_gauges_survive():
    with (DATA / "dewaele-2008-table2-complete.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert Counter(r["material"] for r in rows) == {
        "silver": 34,
        "molybdenum": 29,
        "cobalt_hcp": 25,
        "nickel": 42,
        "zinc_hcp": 43,
    }
    assert len(experiment_rows("silver")) == 31
    assert len(experiment_rows("molybdenum")) == 26
    assert len(experiment_rows("fe")) == 10
    he = [
        r
        for r in rows
        if r["material"] == "nickel" and r["pressure_from_helium_eos"] == "True"
    ]
    assert len(he) == 13
    assert he[0]["pressure_mao_token"] == "(63.2)"
    assert float(he[0]["ruby_pressure_dewaele_gpa"]) == 65
    with (DATA / "dewaele-2008-table3-fe.csv").open() as stream:
        fe = list(csv.DictReader(stream))
    assert len(fe) == 45
    tungsten = [r for r in fe if r["pressure_w_gpa"]]
    assert len(tungsten) == 12
    assert float(tungsten[0]["a_w_angstrom"]) == 3.0833


def test_all_table4_fits_recover_printed_95pct_intervals_with_source_constraints():
    outcomes = reproduce()["records"]
    assert len(list(table4_specs())) == 13
    for material, rid, parameters, _, fixed in table4_specs():
        row = outcomes[rid]
        assert row["status"] == "parity"
        assert row["source_error_confidence"] == 0.95
        assert row["original_author_workflow_parity"] is False
        assert row["fixed_parameters"] == fixed
        assert all(p["within_reported_95pct"] for p in row["parameters"])
        assert all(p["refit_error"] is None for p in row["parameters"])
        native = get_eos_record(rid)
        z = 4 if material in {"nickel", "silver"} else 2
        volume = np.array([0.9, 0.8]) * parameters[0]
        np.testing.assert_allclose(
            native.pressure(volume * z), vinet(volume, parameters), rtol=1e-12
        )
        np.testing.assert_allclose(
            native.volume(native.pressure(volume * z)), volume * z, rtol=1e-10
        )
        assert (
            native.pressure_calibration["methods"][0]["reference_calibration_record"]
            == "ruby_dewaele_2004"
        )
    ag = outcomes["silver_dewaele_2008_table4_vinet"]
    assert ag["published_rmse_gpa"] == pytest.approx(0.3580823795)
    assert ag["rmse_gpa"] == pytest.approx(0.3564005583)


def test_historical_ag_ni_identifiers_keep_numbers_and_correct_the_attribution():
    for material, coefficients in (
        ("silver", [68.28, 100.2, 5.7]),
        ("nickel", [43.816, 177.5, 4.83]),
    ):
        rid = f"{material}_dewaele_2008_vinet_1"
        raw = next(
            r
            for r in get_material_document(material)["eos_records"]
            if r["identifier"] == rid
        )
        assert [
            raw["eos"]["parameters"][k] for k in ("V0", "K0", "K0_prime")
        ] == coefficients
        assert raw["reference"]["year"] == 2019
        assert raw["pressure_scale_variant_of"] == f"{material}_dewaele_2019_mao_vinet"
        legacy = get_eos_record(rid)
        canonical = get_eos_record(f"{material}_dewaele_2019_mao_vinet")
        v = np.array([0.9, 0.8]) * coefficients[0]
        np.testing.assert_array_equal(legacy.pressure(v), canonical.pressure(v))
        assert rid not in {r.identifier for r in search_eos_records()}


def test_paw_states_and_source_figure_limits_remain_explicit():
    assert len(list(table6_specs())) == 54
    outcomes = reproduce()["records"]
    theory = [r for r in outcomes.values() if "_paw_" in r["record_identifier"]]
    assert sum(r["curve_validation"]["passed"] for r in theory) == 50
    assert sum(r["status"] == "parity" for r in theory) == 50
    failed = {
        r["record_identifier"] for r in theory if r["status"] == "parity_not_achieved"
    }
    assert failed == {
        f"iron_dewaele_2008_paw_lda_{state}_vinet" for state in ("static", "298k")
    }
    for material, rid, _, _, state in table6_specs():
        raw = next(
            r
            for r in json.loads((RESEARCH / f"{material}.eosmat").read_text())[
                "eos_records"
            ]
            if r["identifier"] == rid
        )
        assert raw["temperature_ref"] == (0 if state == "static" else 298)
        assert all(v is None for v in (p["error"] for p in raw["parameters"].values()))
        assert not raw.get(
            "fit_datasets"
        )  # graphical curves are not raw DFT fit inputs
        assessment = raw["evidence"]["reproduction"]["assessment"]
        assert assessment["scope"] == "source_equation"
        assert not assessment.get("parameters")
        if outcomes[rid]["curve_validation"]["passed"]:
            assert assessment["status"] == "parity"
            assert (
                assessment["basis_code"] == "source_curve_reading_and_print_precision"
            )
        check = raw["evidence"]["reproduction"]["source_equation_check"]
        assert check["independent_paw_refit_available"] is False
        assert check["passed"] == outcomes[rid]["curve_validation"]["passed"]
        if rid in failed:
            assert (
                raw["evidence"]["scientific_validation"]["status"] == "not_reproduced"
            )
            assert (
                raw["evidence"]["scientific_validation"]["usage_recommendation"]
                == "not_recommended"
            )
            assert "thin path" in raw["evidence"]["scientific_validation"]["note"]
        elif material == "fe" and "_lda_" in rid:
            assert raw["evidence"]["scientific_validation"]["status"] == "deferred"
        else:
            assert (
                raw["evidence"]["scientific_validation"]["status"]
                == "primary_source_validated"
            )


def test_ruby_2008_derivation_does_not_become_table4_input_scale():
    source = json.loads((DATA / "dewaele-2008-source.json").read_text())
    assert source["calibration_lineage"]["table4_input"] == "ruby_dewaele_2004"
    assert (
        source["calibration_lineage"]["new_calibration_output"] == "ruby_dewaele_2008"
    )
    assert len(source["calibration_lineage"]["ten_metal_comparison"]) == 10
    assert "Zn" not in source["calibration_lineage"]["ten_metal_comparison"]


def test_common_ledger_retains_the_dedicated_source_outcomes():
    from scripts.standardize_reconstruction_metadata import assessment_for_record
    from scripts.validate_primary_eos_refits import check_ledger_values

    common = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    research = json.loads((RESEARCH / "reproduction.json").read_text())
    rows = {r["record_identifier"]: r for r in common["records"] + research["records"]}
    for rid, expected in reproduce()["records"].items():
        subset = {key: rows[rid][key] for key in expected}
        check_ledger_values(expected, subset)
        path = (
            ROOT / f"peritheos/data/materials/{expected['material_identifier']}.eosmat"
        )
        if "_paw_" in rid:
            path = RESEARCH / f"{expected['material_identifier']}.eosmat"
        native = next(
            r
            for r in json.loads(path.read_text())["eos_records"]
            if r["identifier"] == rid
        )
        expected_assessment = assessment_for_record(native, rows[rid])
        if "_paw_" in rid:
            import copy

            from scripts.move_dewaele_2008_paw_to_research import relocate_record_links

            temporary = copy.deepcopy(native)
            temporary["evidence"]["reproduction"]["assessment"] = expected_assessment
            relocate_record_links(temporary)
            expected_assessment = temporary["evidence"]["reproduction"]["assessment"]
        assert native["evidence"]["reproduction"]["assessment"] == expected_assessment


@pytest.mark.parametrize("material", ["nickel", "silver"])
def test_historical_compatibility_views_use_the_published_mao_coordinates(material):
    document = get_material_document(material)
    dataset = f"{material}_dewaele_2008_table2_compression"
    compatibility = resolve_dataset_observations(
        document, f"{material}_dewaele_2008_vinet_1", dataset
    )
    parent = resolve_dataset_observations(
        document, f"{material}_dewaele_2019_mao_vinet", dataset
    )
    with (DATA / f"{material}-dewaele-2008-table2-compression.csv").open() as stream:
        source = list(csv.DictReader(stream))
    assert compatibility.pressure_gpa == pytest.approx(
        [float(row["ruby_pressure_mao_gpa"]) for row in source]
    )
    assert np.array_equal(compatibility.pressure_gpa, parent.pressure_gpa)
    assert np.array_equal(compatibility.volume_a3, parent.volume_a3)
