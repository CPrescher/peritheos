"""Verify archived real-console fits; the external EosFit binary is not required."""

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import R

from scripts.investigate_rubidium_halides_farla_2025 import AUTHOR_THETA
from scripts.reproduce_rubidium_halides_farla_2025 import (
    DATA,
    PREFIX,
    SOURCE,
    read_rows,
    reference_pressure,
)
from scripts.run_rubidium_halides_farla_2025_eosfit import (
    CASES,
    PARAMETERS,
    covariance,
    macro,
    prepare_input,
)

OUTPUT = Path(__file__).resolve().parents[1] / "docs/data" / f"{PREFIX}-eosfit7c"
REPORT = json.loads((OUTPUT / "manifest.json").read_text())


@pytest.mark.parametrize("compound", SOURCE)
@pytest.mark.parametrize("name", CASES)
def test_native_inputs_row_selection_model_covariance_and_file_hashes(compound, name):
    folder = OUTPUT / compound / name
    case = REPORT["cases"][compound][name]
    options = CASES[name]
    data, mapping = prepare_input(compound, **options)
    assert data == (folder / "input.dat").read_text()
    assert mapping == json.loads((folder / "row-map.json").read_text())
    assert len(mapping) == case["observations"]
    # Archived macros name their original absolute checkout, including on clones.
    archived_macro = (folder / "run.mcr").read_text()
    original_folder = Path(archived_macro.splitlines()[0].removeprefix("log ")).parent
    assert macro(compound, original_folder, **options) == archived_macro
    for filename, digest in case["files_sha256"].items():
        assert hashlib.sha256((folder / filename).read_bytes()).hexdigest() == digest
    assert case["converged"]
    assert case["convergence_messages"] == 2
    assert set(case["refined_parameters_final_cycle"]) == set(PARAMETERS)
    model = (folder / "fitted.eos").read_text()
    assert "Model =  2," in model and "Order =  3" in model
    assert "Thermal =  7," in model and "Type = Volume" in model
    assert "Pscale =GPa" in model and "Vscale =cm^3/mol" in model
    values = {
        int(k): float(v) for k, v in re.findall(r"Param\s*=\s*(\d+)\s+(\S+)", model)
    }
    assert values[1] == SOURCE[compound][0]
    assert values[11] == (
        AUTHOR_THETA[compound] if options.get("precise_theta") else SOURCE[compound][1]
    )
    assert values[13] == 2 and values[14] == 0 and values[51] == 1
    matrix = covariance(folder / "fitted.eos")
    np.testing.assert_array_equal(matrix, case["saved_parameter_covariance"])
    np.testing.assert_allclose(matrix, matrix.T, atol=1e-8)
    assert np.linalg.eigvalsh(matrix).min() > 0
    errors = [case["refined_parameters_final_cycle"][p]["esd"] for p in PARAMETERS]
    np.testing.assert_allclose(np.sqrt(np.diag(matrix)), errors, atol=1e-5)
    final = (
        (folder / "stdout.txt")
        .read_text(errors="replace")
        .split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1]
    )
    weights = [
        float(fields[-1])
        for line in final.splitlines()
        if len(fields := line.split()) == 15
        and fields[0].isdigit()
        and fields[1] == "1"
    ]
    assert len(weights) == len(mapping)
    assert np.isfinite(weights).all() and min(weights) > 0
    if not options["weighted"]:
        assert set(weights) == {1}


@pytest.mark.parametrize("compound", SOURCE)
def test_real_console_confirms_irls_without_reproducing_all_published_digits(compound):
    case = REPORT["cases"][compound]["weighted_precise_theta"]
    np.testing.assert_array_less(
        abs(np.array(case["eosfit_minus_independent_irls"])),
        [0.0005, 0.0001, 0.0001, 0.0001],
    )
    comparison = case["published_comparison"]
    assert comparison["maximum_difference_over_published_esd"] < 1.65
    assert not comparison["all_coefficients_round_to_publication"]
    rows = [
        r
        for r in read_rows(compound, "observations")
        if r["volume_a3"]
        and (r["source_table"] == "S1" or float(r["temperature_k"]) > 300)
    ]
    v, t, p = [
        np.array([float(r[k]) for r in rows])
        for k in ("volume_a3", "temperature_k", "pressure_gpa")
    ]
    molar, _, coefficients = SOURCE[compound]
    theta = AUTHOR_THETA[compound]
    pressure = reference_pressure(v, t, molar, theta, coefficients)
    cold = reference_pressure(v, np.full_like(t, 300), molar, theta, coefficients)
    expected_residual = cold + (pressure - cold) * 8.314 / R - p
    native_residual = [r["residual_gpa"] for r in case["initial_source_replay"]["rows"]]
    np.testing.assert_allclose(native_residual, expected_residual, rtol=0, atol=6e-5)


def test_executable_identity_and_source_bytes_remain_explicit():
    assert REPORT["version"] == "7.60"
    assert (
        REPORT["executable_sha256"]
        == "6b48e46e5bb01c2dd13060bc29b4a9f7c67e70c7161866b77f8de3f238445e6d"
    )
    assert not REPORT["runtime_provenance"]["console_executable_modified"]
    assert not REPORT["original_author_input_files"]
    assert not REPORT["original_author_fit_parity"]
    assert REPORT["author_contact"] == "none"
    for compound, digest in REPORT["source_observation_hashes"].items():
        path = DATA / f"{PREFIX}-{compound}-observations.csv"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
