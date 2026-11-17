"""Source preservation and independent properties of the grossular audit."""

import hashlib
import json
from collections import Counter

import numpy as np
import pytest

from scripts.audit_garnet_angel_2022 import FOLDER, independent_pressure, source
from scripts.audit_grossular_angel_2022_observations import (
    DATA,
    GAS_CONSTANT,
    audit,
    observations,
)
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_subset import (
    arrays,
    cold,
    prepared_rows,
    pressure,
    solve,
    thermal,
)
from scripts.run_grossular_angel_2022_eosfit import check_output


def test_original_tokens_missing_errors_and_nominal_inventory_are_preserved():
    manifest = json.loads((DATA / "input-recovery.json").read_text())
    for item in manifest["sources"]:
        path = DATA / item["transcription"]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == item["transcription_sha256"]
        )
        rows = observations(item["transcription"])
        assert len(rows) == item["rows"]
        assert [int(r["source_row"]) for r in rows] == list(range(1, len(rows) + 1))
    assert sum(s["nominal_s1_used_count"] for s in manifest["sources"]) == 161
    assert sum(s["nominal_s1_used_count"] for s in manifest["unrecovered"]) == 59
    assert manifest["count_qualification"]
    milani = observations("milani-2017-table-3.csv")
    assert milani[-1]["coordinate_token"] == "7.674(19)"
    assert observations("milani-2017-table-2.csv")[3]["lattice_a_token"] == "11.841(13)"
    isaak = observations("isaak-1992-table-3b.csv")
    assert sum(bool(r["ks_gpa_reported_error"]) for r in isaak) == 3
    assert all(r["temperature_error_k"] == "" for r in isaak)
    expansion = observations("isaak-1992-table-2.csv")
    assert all("GR-2B" in r["observed_specimen"] for r in expansion)
    assert expansion[0]["observed_volume_strain_x1e2_token"] == ".013"
    assert expansion[0]["calculated_volume_strain_x1e2_token"] == ".014"
    assert all(
        r["temperature_error_k"] == "" for r in observations("kono-2010-table-1.csv")
    )


def test_grossular_source_and_equation_checks_remain_independent_of_subset_fit():
    source_path = FOLDER / "sources/grossular-official.eos"
    assert (
        hashlib.sha256(source_path.read_bytes()).hexdigest()
        == "279139b3835582e97e5419283271151a9b10491f7eef89e986ce4e84e9d3dca5"
    )
    receipt = json.loads((DATA / "eosfit-properties/checks.json").read_text())
    assert receipt["all_within_tolerances"]
    replay = check_output((DATA / "eosfit-properties/stdout.txt").read_text())
    assert len(replay["states"]) == len(receipt["states"]) == 9
    for actual, expected in zip(replay["states"], receipt["states"]):
        assert actual.keys() == expected.keys()
        for key in expected:
            if key in {"independent", "official_minus_independent"}:
                # Property evaluation has cross-library roundoff, including
                # cancellation in near-zero differences. These bounds remain
                # well below the receipt's declared property tolerances.
                assert actual[key] == pytest.approx(expected[key], rel=1e-10, abs=1e-12)
            else:
                assert actual[key] == expected[key]
    result = audit()
    assert result["source_nonzero_covariance_parameters"] == [2, 3, 11, 18, 53, 54]
    assert not result["source_coefficients_changed"]
    assert not result["author_workflow_reproduced"]
    assert not result["isaak_expansion_identity"]["included_in_endmember_refit"]
    assert (
        result["moduli"]["kono"]["calculated_modulus_multiplier_to_observed_scale"]
        == 0.99455
    )
    assert result["milani_compression"]["residual_summary_cm3_mol"]["count"] == 7


def test_fast_refit_equation_and_ks_inversion_match_independent_quadrature():
    p, _ = source("grossular")
    source_c = [p[2], p[3], p[11], p[18], 1, 1, 1]
    report = json.loads((DATA / "subset-refits.json").read_text())
    assert (
        hashlib.sha256((DATA / "refit-plan.json").read_bytes()).hexdigest()
        == report["plan_sha256"]
    )
    for name, digest in report["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    fitted_c = list(report["fits"][0]["parameters"].values())
    for c in (source_c, fitted_c):
        pp = p | dict(zip([2, 3, 11, 18], c[:4]))
        for t, pressure_gpa in [
            (20, 0.0001),
            (100, 0.0001),
            (298, 7.674),
            (1000, 0.0001),
            (1650, 16.91),
        ]:
            values = properties(pp, pressure_gpa, t, GAS_CONSTANT)
            v, ks = values["volume_cm3_mol"], values["ks_gpa"]
            pth, _ = thermal(t, c)
            np.testing.assert_allclose(
                cold(v, c)[0] + pth[0],
                independent_pressure(pp, v, t, GAS_CONSTANT),
                atol=1e-9,
            )
            # Both data types and a non-unit calculated-to-observed scale.
            rows = [
                dict(
                    P=pressure_gpa,
                    T=t,
                    Y=y * 1.003,
                    sY=0.01,
                    sT=None,
                    sP=None,
                    KS=k,
                    scale=4,
                )
                for y, k in [(v, False), (ks, True)]
            ]
            cc = np.array(c)
            cc[4] = 1.003
            np.testing.assert_allclose(
                pressure(cc, arrays(rows)), pressure_gpa, atol=1e-8
            )


def test_conditional_fit_selection_convergence_and_scope():
    rows = prepared_rows()
    assert Counter(r["dataset"] for r in rows) == dict(
        milani_vt=35, du_vt=13, milani_pv=7, kono_ks=36, isaak_ks=3, bass_ks=1
    )
    assert all(
        r["sT"] is None
        for r in rows
        if r["dataset"] in {"kono_ks", "isaak_ks", "bass_ks"}
    )
    report = json.loads((DATA / "subset-refits.json").read_text())
    assert report["channels"] == 95
    assert not report["author_workflow_reproduced"]
    assert not report["source_coefficients_changed"]
    a, b = report["fits"]
    assert a["numerically_converged"] and b["numerically_converged"]
    assert not a["at_parameter_bound"] and not b["at_parameter_bound"]
    np.testing.assert_allclose(
        list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-6
    )
    assert 4 < a["conditional_reduced_chi_square"] < 5
    p, _ = source("grossular")
    replay = solve(rows, [p[2], p[3], p[11], p[18], 1, 1, p[53]])
    np.testing.assert_allclose(
        list(replay["parameters"].values()), list(a["parameters"].values()), rtol=1e-6
    )
