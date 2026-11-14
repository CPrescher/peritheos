"""Audit the follow-up grossular recovery without fitting or changing the EOS.

Repository-only: python -m scripts.audit_grossular_angel_2022_recovery.
The 161-channel audit and 95-channel refit are retained as historical receipts.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import (
    DATA,
    GAS_CONSTANT,
    observations,
)
from scripts.audit_pyrope_angel_2022_observations import properties, summarize


def fingerprint(body):
    """FNV-1a matches the literal 624-cell body extracted from publisher HTML."""
    value = 2166136261
    for character in body:
        value = ((value ^ ord(character)) * 16777619) & 0xFFFFFFFF
    return f"{value:08x}"


def audit():
    p, _ = source("grossular")
    file = DATA / "gwanmesia-2014-table-1.csv"
    body = "\n".join(file.read_text().splitlines()[1:]) + "\n"
    if fingerprint(body) != "c93742ca":
        raise ValueError("Gwanmesia transport fingerprint differs from publisher table")
    rows = []
    for i, r in enumerate(observations(file.name), 1):
        pressure, t = float(r["pressure_gpa"]), float(r["temperature_c"]) + 273.15
        observed = float(r["ks_gpa"])
        calculated = properties(p, pressure, t, GAS_CONSTANT)["ks_gpa"]
        rows.append(
            {
                "source_row": i,
                "pressure_gpa": pressure,
                "temperature_k": t,
                "observed_ks_gpa": observed,
                "calculated_ks_gpa": calculated,
                "calculated_ks_on_observed_scale_gpa": p[54] * calculated,
                "observed_minus_scaled_calculated_gpa": observed - p[54] * calculated,
                "observed_minus_unscaled_calculated_gpa": observed - calculated,
                "reported_relative_ks_width_gpa": 0.008 * observed,
                "printed_modulus_identity_difference_gpa": observed
                - (
                    float(r["longitudinal_modulus_gpa"])
                    - 4 * float(r["shear_modulus_gpa"]) / 3
                ),
            }
        )
    figure = json.loads((DATA / "figure4-recovery.json").read_text())
    original = observations("isaak-1992-table-2.csv")
    plotted = observations("figure4-isaak.csv")
    comparisons = {}
    for kind in ("observed", "calculated"):
        ratios = np.array(
            [1 + 0.01 * float(r[f"{kind}_volume_strain_x1e2_token"]) for r in original]
        )
        volumes = np.array([float(r["molar_volume_cm3_mol"]) for r in plotted])
        comparisons[kind] = {
            "fixed_reference_volume_cm3_mol": 125.35,
            "fixed_reference_residual_cm3_mol": summarize(volumes - 125.35 * ratios),
        }
    tdiff = [
        float(a["temperature_k"]) - float(b["temperature_k"])
        for a, b in zip(plotted, original)
    ]
    result = {
        "audit_date": "2026-11-12",
        "scope": "Follow-up observation recovery and fixed-source comparison; no new fit",
        "historical_receipts": [
            "input-recovery.json",
            "observation-audit.json",
            "subset-refits.json",
            "verification.json",
        ],
        "coverage": {
            "nominal_s1_used_channels": 220,
            "tabulated_or_table_derived_channels": 213,
            "additional_digitized_processed_channels": 7,
            "original_thieblot_selected_table_channels_still_missing": 7,
            "qualification": "161 earlier upstream channels plus 52 Gwanmesia table rows plus seven Thiéblot plotted points cover the nominal S1 inventory at mixed evidence levels. This is not 220 recovered original author-fit residuals. Anchor counts and selection prose conflict; original normalization, errors and fit inputs remain incomplete.",
        },
        "gwanmesia": {
            "rows_recovered": len(rows),
            "calculated_modulus_multiplier_to_observed_scale": p[54],
            "scaled_residual_gpa": summarize(
                [r["observed_minus_scaled_calculated_gpa"] for r in rows]
            ),
            "unscaled_residual_gpa": summarize(
                [r["observed_minus_unscaled_calculated_gpa"] for r in rows]
            ),
            "maximum_printed_modulus_identity_difference_gpa": max(
                abs(r["printed_modulus_identity_difference_gpa"]) for r in rows
            ),
            "error_semantics": "Reported 0.8% KS error estimate; confidence and covariance unspecified. Not asserted to be independent one-sigma uncertainty. No temperature error supplied. Pressure precision 0.7%; not propagated here. No formal chi-square or p-value computed.",
            "xrd_volume_constraints_added": False,
            "rows": rows,
        },
        "isaak": {
            "primary_reference_temperature_k": 293,
            "reference_evidence": "Isaak et al. (1992), p.109 explicitly defines the strain reference Tref=293 K; this resolves the earlier audit's reference-temperature assumption.",
            "measured_dilatometry_specimen": "GR-2B: 76.1 mol% grossular, 22.1 mol% andradite",
            "s1_composition_belongs_to": "GGD elasticity specimen, about 96.53 mol% grossular",
            "plot_comparisons": comparisons,
            "plot_minus_table_temperature_k": summarize(tdiff),
            "finding": "All 15 Figure 4 points match 125.35 times the GR-2B Table 2 measured volume-strain ratios within the declared digitization bounds. Thus the plotted series is consistent with GR-2B strain scaled directly to the adopted volume, despite S1's GGD composition label. The nearly identical calculated GR-2B strain column cannot be uniquely distinguished at plot resolution. No resolved compositional correction or 293-to-298 K re-referencing is evident from this agreement.",
            "limitation": "This identifies the plotted series, not the exact fit-project inputs. The endmember expansivity estimate is calculated, not a measured pure-grossular volume series. These mixed-composition observations are not added to the diagnostic endmember refit.",
            "included_in_endmember_refit": False,
        },
        "figure4_recovery": {
            "receipt": "figure4-recovery.json",
            "selected_thieblot_points": figure["groups"]["thieblot_selected"]["rows"],
            "excluded_visible_thieblot_points": figure["groups"][
                "thieblot_excluded_visible"
            ]["rows"],
            "selection_conflict": "Main article p.8 says Skinner, Thieblot and Rodehorst (2004) expansion are excluded. S1 and EOS comments retain seven Thieblot points below 1000 K; Figure 4 shows seven filled and six open magenta markers. The exact author mask remains unconfirmed.",
        },
        "source_coefficients_changed": False,
        "author_workflow_reproduced": False,
        "new_fit_performed": False,
        "author_contacted": False,
        "input_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [
                file.name,
                "gwanmesia-2014-source.json",
                "figure4-recovery.json",
                "figure4-isaak.csv",
                "figure4-thieblot-selected.csv",
                "figure4-thieblot-excluded-visible.csv",
                "isaak-1992-table-2.csv",
            ]
        },
    }
    return result


if __name__ == "__main__":
    (DATA / "followup-recovery.json").write_text(
        json.dumps(audit(), indent=2, allow_nan=False) + "\n"
    )
