"""Case 3 native-execution evidence and conditional anchor-rounding diagnostics."""

import hashlib
import json
import re
from pathlib import Path, PurePosixPath

import numpy as np
import pytest

from scripts.run_rubidium_halides_farla_2025_case3_eosfit import (
    CASES,
    INPUTS,
    PARAMETERS,
    ROUNDING_ANCHOR_EXAMPLES,
    VARIANTS,
    commands,
)
from scripts.run_rubidium_halides_farla_2025_eosfit import covariance, prepare_input

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/rubidium-halides-farla-2025-case3-eosfit7c"
REPORT = json.loads((OUTPUT / "manifest.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize("name", CASES)
def test_native_case3_models_keep_four_anchors_fixed_and_preserve_evidence(
    variant, name
):
    folder = OUTPUT / variant / name
    case = REPORT["variants"][variant]["cases"][name]
    source = VARIANTS[variant]
    data, mapping = prepare_input(source["compound"], **CASES[name])
    assert (folder / "input.dat").read_text(encoding="utf-8") == data
    assert json.loads((folder / "row-map.json").read_text(encoding="utf-8")) == mapping
    archived = (folder / "run.mcr").read_text(encoding="utf-8")
    original_folder = PurePosixPath(
        archived.splitlines()[0].removeprefix("log ")
    ).parent
    assert commands(variant, original_folder, CASES[name]) == archived
    for filename, digest in case["files_sha256"].items():
        assert hashlib.sha256((folder / filename).read_bytes()).hexdigest() == digest
    assert case["converged"] and case["convergence_messages"] == 2
    assert set(case["refined_parameters_final_cycle"]) == set(PARAMETERS)
    for filename in ("initial.eos", "fitted.eos"):
        text = (folder / filename).read_text(encoding="utf-8")
        params = {
            int(k): float(v) for k, v in re.findall(r"Param\s*=\s*(\d+)\s+(\S+)", text)
        }
        actual = [params[1], params[2], params[3], params[11]]
        np.testing.assert_allclose(actual, case["fixed_parameters"], rtol=0, atol=1e-5)
        assert params[13] == 2 and params[14] == 0 and params[51] == 1
        assert "Thermal =  7," in text and "Order =  3" in text
    matrix = covariance(folder / "fitted.eos", indices=(17, 18))
    np.testing.assert_array_equal(matrix, case["saved_thermal_covariance"])
    assert np.linalg.eigvalsh(matrix).min() > 0
    errors = [case["refined_parameters_final_cycle"][p]["esd"] for p in PARAMETERS]
    np.testing.assert_allclose(np.sqrt(np.diag(matrix)), errors, atol=1e-5)
    np.testing.assert_array_equal(
        covariance(folder / "fitted.eos", indices=(0, 1, 2, 10)), np.zeros((4, 4))
    )
    if variant == "rbbr_kprime_4p98":
        assert case["fixed_parameters"][2] == 4.98


@pytest.mark.parametrize("variant", VARIANTS)
def test_nominal_case3_thermal_parameters_are_close_without_exact_rounding(variant):
    case = REPORT["variants"][variant]["cases"]["weighted_printed_theta"]
    comparison = case["published_comparison"]
    assert comparison["maximum_absolute_difference_over_published_esd"] < 0.09
    assert not comparison["both_parameters_round_to_publication"]
    assert case["fixed_parameters"] == VARIANTS[variant]["fixed_parameters"]


@pytest.mark.parametrize("variant", VARIANTS)
def test_sub_last_digit_anchor_examples_establish_conditional_rounding_consistency(
    variant,
):
    case = REPORT["variants"][variant]["cases"]["weighted_anchor_rounding_example"]
    anchors = ROUNDING_ANCHOR_EXAMPLES[variant]
    source = VARIANTS[variant]["fixed_parameters"][:3]
    assert all(round(a, 2) == b for a, b in zip(anchors, source))
    assert max(abs(np.array(anchors) - source)) < 0.005
    assert case["published_comparison"]["both_parameters_round_to_publication"]
    assert "not recovered author values" in case["fixed_anchor_qualification"]


@pytest.mark.parametrize("variant", VARIANTS)
def test_300k_rows_do_not_move_the_two_free_thermal_coefficients(variant):
    cases = REPORT["variants"][variant]["cases"]
    expected = cases["weighted_printed_theta"]["published_comparison"][
        "eosfit_parameters"
    ]
    for name in ["weighted_all_lvp_printed_theta", "weighted_lvp_only_printed_theta"]:
        np.testing.assert_allclose(
            cases[name]["published_comparison"]["eosfit_parameters"],
            expected,
            atol=1e-5,
            rtol=0,
        )


def test_source_provenance_keeps_case1_case3_and_author_recovery_distinct():
    assert REPORT["version"] == "7.60"
    assert (
        REPORT["source_facts_sha256"] == hashlib.sha256(INPUTS.read_bytes()).hexdigest()
    )
    assert (
        REPORT["executable_sha256"]
        == "6b48e46e5bb01c2dd13060bc29b4a9f7c67e70c7161866b77f8de3f238445e6d"
    )
    assert REPORT["rounding_examples_all_match_publication"]
    assert not REPORT["author_fit_parity"]
    assert not REPORT["original_author_input_files_recovered"]
    assert REPORT["author_contact"] == "none"
