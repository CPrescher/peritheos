"""Evidence levels, source tokens and numerical scope of the grossular recovery."""

import hashlib
import json
from collections import Counter

import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA, observations
from scripts.audit_grossular_angel_2022_recovery import audit, fingerprint


def test_primary_gwanmesia_tokens_units_and_missing_error_semantics():
    file = DATA / "gwanmesia-2014-table-1.csv"
    assert (
        fingerprint("\n".join(file.read_text().splitlines()[1:]) + "\n") == "c93742ca"
    )
    rows = observations(file.name)
    assert len(rows) == 52
    assert all(len(row) == 12 for row in rows)
    assert rows[0]["cell_volume_a3"] == "1656.20"
    assert rows[0]["tp_microsecond"] == "0.2333"
    assert rows[-1]["ks_gpa"] == "205.6"
    assert Counter(r["temperature_c"] for r in rows) == {
        "30": 10,
        **{str(t): 6 for t in range(127, 728, 100)},
    }
    result = audit()
    assert result["gwanmesia"]["rows"][0]["temperature_k"] == 303.15
    assert result["gwanmesia"]["rows"][-1]["temperature_k"] == 1000.15
    assert not result["gwanmesia"]["xrd_volume_constraints_added"]
    metadata = json.loads((DATA / "gwanmesia-2014-source.json").read_text())
    assert metadata["reported_error_estimates_percent"]["ks"] == 0.8
    assert (
        "not reconstructed per-row standard deviations"
        in metadata["error_qualification"]
    )
    # Rounding each of L, G and KS to 0.1 GPa allows this maximum discrepancy.
    assert result["gwanmesia"]["maximum_printed_modulus_identity_difference_gpa"] <= (
        0.05 + 4 * 0.05 / 3 + 0.05
    )
    assert 0.71 < result["gwanmesia"]["scaled_residual_gpa"]["rms"] < 0.72
    assert (
        result["gwanmesia"]["calculated_modulus_multiplier_to_observed_scale"]
        == 1.00551
    )


def test_figure_evidence_and_isaak_identity_are_qualified_at_plot_resolution():
    receipt = json.loads((DATA / "figure4-recovery.json").read_text())
    for item in receipt["groups"].values():
        assert (
            hashlib.sha256((DATA / item["transcription"]).read_bytes()).hexdigest()
            == item["sha256"]
        )
    selected = observations("figure4-thieblot-selected.csv")
    rejected = observations("figure4-thieblot-excluded-visible.csv")
    assert len(selected) == 7 and len(rejected) == 6
    assert all(float(r["temperature_k"]) < 1000 for r in selected)
    assert all(float(r["temperature_k"]) > 1000 for r in rejected)
    bounds = receipt["conservative_digitization_bounds"]
    assert bounds["temperature_k"] == 2 and bounds["molar_volume_cm3_mol"] == 0.003
    original = observations("isaak-1992-table-2.csv")
    plotted = observations("figure4-isaak.csv")
    assert len(original) == len(plotted) == 15
    for a, b in zip(original, plotted):
        assert abs(float(a["temperature_k"]) - float(b["temperature_k"])) < 2
        expected = 125.35 * (1 + 0.01 * float(a["observed_volume_strain_x1e2_token"]))
        assert abs(expected - float(b["molar_volume_cm3_mol"])) < 0.003
    assert all("GR-2B" in r["observed_specimen"] for r in original)
    result = audit()["isaak"]
    assert result["primary_reference_temperature_k"] == 293
    assert not result["included_in_endmember_refit"]
    # The calculated strain column also agrees: plot resolution cannot identify
    # a unique Table 2 column or the original fit-project contents.
    for kind in ("observed", "calculated"):
        assert (
            result["plot_comparisons"][kind]["fixed_reference_residual_cm3_mol"][
                "max_absolute"
            ]
            < 0.003
        )


def test_followup_preserves_historical_receipts_and_does_not_claim_full_fit():
    old = json.loads((DATA / "artifact-manifest.json").read_text())
    for item in old["files"]:
        assert (
            hashlib.sha256((DATA / item["file"]).read_bytes()).hexdigest()
            == item["sha256"]
        )
    saved = json.loads((DATA / "followup-recovery.json").read_text())
    result = audit()
    assert saved == result
    assert result["coverage"]["tabulated_or_table_derived_channels"] == 213
    assert result["coverage"]["additional_digitized_processed_channels"] == 7
    assert (
        result["coverage"]["original_thieblot_selected_table_channels_still_missing"]
        == 7
    )
    assert not any(
        result[k]
        for k in [
            "source_coefficients_changed",
            "author_workflow_reproduced",
            "new_fit_performed",
            "author_contacted",
        ]
    )
    for name, digest in result["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    gwan = result["gwanmesia"]
    residuals = [
        r["observed_ks_gpa"] - 1.00551 * r["calculated_ks_gpa"] for r in gwan["rows"]
    ]
    np.testing.assert_allclose(
        np.sqrt(np.mean(np.square(residuals))), gwan["scaled_residual_gpa"]["rms"]
    )
