"""Verify the stated two-stage procedure and plotted numerical comparisons."""

import json

import numpy as np
import pytest

from peritheos.eos.rt import BM3
from scripts.brucite_replay_checks import check_report
from scripts.brucite_reproduce import observations, reference_volume, source
from scripts.brucite_staged_comparison import REPORT, reproduce


@pytest.fixture(scope="module")
def report():
    result = reproduce()
    check_report(result, json.loads(REPORT.read_text(encoding="utf-8")))
    return result


def test_ambient_stage_reproduces_the_printed_constant_law(report):
    ambient = report["stage1_ambient"]
    assert ambient["solver_success"]
    # Supported solver builds differ by a few 1e-9 Å³ in this diagnostic fit.
    # These absolute replay bounds do not change archived or published values.
    assert ambient["V300_a3"] == pytest.approx(40.9165851955, rel=0, abs=1e-8)
    assert ambient["alpha_per_k"] == pytest.approx(6.69185183775e-5, rel=0, abs=1e-13)
    assert ambient["source_rows"] == list(range(1, 12))
    assert not report["author_workflow_reproduced"]


@pytest.mark.parametrize(
    "index,temperature,count", [(0, 300, 13), (1, 473, 8), (2, 673, 8), (3, 873, 8)]
)
def test_fixed_reference_stage_and_native_residuals(report, index, temperature, count):
    result = report["isotherms"][index]
    saved = json.loads(REPORT.read_text(encoding="utf-8"))["isotherms"][index]
    rows = [
        row
        for row in observations("fukui_2003", "table4")
        if float(row["temperature_k"]) == temperature
    ]
    v = np.array([float(row["cell_volume_a3"]) for row in rows])
    p = np.array([float(row["pressure_gpa"]) for row in rows])
    assert result["solver_success"]
    assert result["observations"] == count == len(rows)
    assert result["source_rows"] == [int(row["source_row"]) for row in rows]
    ambient = report["stage1_ambient"]
    assert result["staged_refit"]["V0_a3"] == ambient["V300_a3"] * np.exp(
        ambient["alpha_per_k"] * (temperature - 300)
    )
    assert result["published"]["V0_a3"] == reference_volume(temperature)
    published = source()["sources"]["fukui_2003"]["parameters"]["isotherms"][index]
    for parameter in ("K0_gpa", "K0_prime"):
        assert result["published"][parameter] == published[parameter]
    for key in ["published", "staged_refit"]:
        curve = result[key]
        native = BM3(curve["V0_a3"], curve["K0_gpa"], curve["K0_prime"])
        np.testing.assert_allclose(
            native.pressure(v) - p, curve["pressure_residuals_gpa"], atol=1e-11
        )
        predicted_v = native.volume(p)
        np.testing.assert_allclose(
            predicted_v - v, curve["volume_residuals_a3"], atol=1e-9
        )
        assert np.sqrt(np.mean((predicted_v - v) ** 2)) == pytest.approx(
            curve["volume_rmse_a3"], abs=1e-10
        )
        assert curve["pressure_rmse_gpa"] == pytest.approx(
            saved[key]["pressure_rmse_gpa"], rel=0, abs=1e-9
        )


def test_873k_improvement_requires_no_change_to_the_expansion_law(report):
    hottest = report["isotherms"][-1]
    assert hottest["staged_refit"]["K0_gpa"] == pytest.approx(25.475924, abs=1e-5)
    assert hottest["staged_refit"]["pressure_rmse_gpa"] == pytest.approx(
        0.15803249, abs=1e-8
    )
    assert hottest["staged_refit"]["volume_rmse_a3"] == pytest.approx(
        0.05865854, abs=1e-8
    )
    assert hottest["published"]["volume_rmse_a3"] == pytest.approx(0.17072962, abs=1e-8)
