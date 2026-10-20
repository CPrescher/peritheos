"""Actual GUI artifacts, retained input rows and the bounded scientific finding."""

import hashlib
import json

import numpy as np
import pytest

from scripts.audit_fes_gui_replay import AUTHOR, FOLDER, model, reproduce


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_gui_uses_confirmed_selection_and_preserves_every_input_token():
    lines = (AUTHOR / "FeS6EOSfittot-SataOhfuji-300K.dat").read_bytes().splitlines()
    for name, hot, count in (("cold", False, 21), ("hot", True, 146)):
        expected = [
            r for r in lines[1:] if r.strip() and (float(r.split()[2]) > 300) == hot
        ]
        actual = (FOLDER / (name + ".dat")).read_bytes().splitlines()
        assert actual == [lines[0], *expected]
        assert len(expected) == count
    assert lines[7] == lines[8]  # Duplicate thermal source observation is preserved.


def test_gui_saved_fit_and_console_comparison_remain_qualified(report):
    weighted = report["cases"]["gamma-weighted"]
    gamma = weighted["refined_parameters"]["18"]
    assert gamma["value"] == pytest.approx(2.46316)
    assert gamma["error"] == pytest.approx(0.03949, abs=1e-5)
    staged = report["cases"]["staged-weighted"]
    assert staged["refined_parameters"]["18"]["value"] == pytest.approx(2.46420)
    assert staged["refit"]["thermal"]["max_abs_gpa"] > 5
    assert len(staged["refit"]["thermal"]["rows_exceeding_3_gpa"]) == 12
    assert not report["published_residual_envelope_reproduced"]
    assert report["gui_version"] == "20211014"
    assert "20210609" in report["limitation"]
    for c in report["cases"].values():
        assert c["gui_saved_native_max_difference_from_display_gpa"] < 0.02


def test_staged_thermal_initial_keeps_exact_saved_cold_tokens_and_full_mgd():
    _, cold, _ = model(FOLDER / "cold-weighted.eos")
    _, thermal, tref = model(FOLDER / "thermal-initial.eos")
    assert [thermal[i] for i in (1, 2, 3)] == [cold[i] for i in (1, 2, 3)]
    assert [thermal[i] for i in (11, 13, 14, 19)] == [417, 2, 0, 1]
    assert tref == 300
    # The actual saved GUI covariance retains one refined thermal parameter.
    report = json.loads((FOLDER / "manifest.json").read_text(encoding="utf-8"))
    covariance = np.asarray(
        report["cases"]["staged-weighted"]["conditional_covariance"]
    )
    assert covariance.shape == (1, 1)
    assert covariance[0, 0] > 0
    for name, checksum in report["artifact_sha256"].items():
        from scripts.audit_fes_gui_replay import ROOT

        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == checksum
