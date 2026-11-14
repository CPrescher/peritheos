"""Audit recovered grossular tables against the unchanged Angel (2022) EOS.

Repository-only: python -m scripts.audit_grossular_angel_2022_observations.
Profiled volume scales are conditional curve comparisons, not the author's
normalization or a refit of the EOS. Missing errors remain missing in CSVs.
"""

from __future__ import annotations

import csv
import json

import numpy as np
from scipy.constants import Avogadro

from scripts.audit_garnet_angel_2022 import FOLDER, source
from scripts.audit_pyrope_angel_2022_observations import properties, summarize
from scripts.refit_pyrope_angel_2022_expanded import recalibrate

DATA = FOLDER / "grossular"
GAS_CONSTANT = float(np.float32(8.314))


def observations(filename):
    with (DATA / filename).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def volume_rows():
    """Original units converted exactly; no measured value is replaced by EOS."""
    groups = {}
    for name, filename in [
        ("milani_vt", "milani-2017-table-2.csv"),
        ("du_vt", "du-2015-preprint-table-2.csv"),
        ("bosenick_vt", "bosenick-geiger-1997-table-2.csv"),
        ("rodehorst_vt", "rodehorst-2002-table-4a.csv"),
    ]:
        rows = []
        for r in observations(filename):
            t = float(r["temperature_k"])
            if name == "rodehorst_vt" and t == 298:
                continue
            if name == "bosenick_vt":
                v, sv = (
                    float(r[k])
                    for k in ["molar_volume_cm3_mol", "molar_volume_sd_cm3_mol"]
                )
            elif name == "rodehorst_vt":
                a, sa = (
                    float(r["lattice_a_angstrom"]),
                    float(r["lattice_a_sd_angstrom"]),
                )
                v, sv = a**3 * Avogadro / 8e24, 3 * a**2 * sa * Avogadro / 8e24
            else:
                v, sv = (
                    float(r[k]) * Avogadro / 8e24
                    for k in ["cell_volume_a3", "cell_volume_sd_a3"]
                )
            rows.append({"row": int(r["source_row"]), "T": t, "V": v, "sV": sv})
        groups[name] = rows
    return groups


def audit():
    p, covariance = source("grossular")
    result = {
        "audit_date": "2026-11-12",
        "scope": "partial_original_observation_audit",
        "source_coefficients_changed": False,
        "author_workflow_reproduced": False,
        "nominal_s1_used_count": 220,
        "nominal_channels_recovered": 161,
        "nominal_channels_unrecovered": 59,
        "count_qualification": "S1 used counts sum to 220, but its normalization-anchor footnote conflicts with some unchanged original/used counts. Recovered counts identify upstream channels, not 161 confirmed author-fit residuals. Isaak expansion specimen identity is unresolved.",
        "gas_constant_j_mol_k": GAS_CONSTANT,
        "source_nonzero_covariance_parameters": [
            int(i + 1) for i in np.flatnonzero(np.diag(covariance))
        ],
        "volume_shape_comparisons": {},
    }
    for name, rows in volume_rows().items():
        g = np.array(
            [
                properties(p, 0.0001, r["T"], GAS_CONSTANT)["volume_cm3_mol"] / p[1]
                for r in rows
            ]
        )
        v, sv = (np.array([r[k] for r in rows]) for k in ["V", "sV"])
        scale = float(np.sum(g * v / sv**2) / np.sum(g**2 / sv**2))
        residual = v - scale * g
        result["volume_shape_comparisons"][name] = {
            "protocol": "Fixed source EOS; one constant molar-volume scale profiled with diagonal reported volume errors only. Temperature errors and shared calibration/anchor covariance are omitted. This is a conditional shape check, not the original normalization or a formal goodness-of-fit test.",
            "profiled_volume_at_298k_cm3_mol": scale,
            "residual_summary_cm3_mol": summarize(residual),
            "volume_error_only_standardized_rms": float(
                np.sqrt(np.mean((residual / sv) ** 2))
            ),
            "rows": [
                {
                    **r,
                    "calculated_volume_cm3_mol": float(scale * gi),
                    "observed_minus_calculated_cm3_mol": float(ri),
                }
                for r, gi, ri in zip(rows, g, residual)
            ],
        }
    pv = observations("milani-2017-table-3.csv")
    anchor = float(pv[0]["cell_volume_a3"])
    # The room-temperature table supplies no numerical T: 298 K is explicit
    # diagnostic adoption, not a recovered experimental thermometer reading.
    pv_rows = []
    for r in pv[1:]:
        original_p, original_sp = float(r["pressure_gpa"]), float(r["pressure_sd"])
        pressure, sp, quartz_v = recalibrate(original_p, original_sp)
        calc = properties(p, pressure, 298, GAS_CONSTANT)["volume_cm3_mol"]
        observed = p[1] * float(r["cell_volume_a3"]) / anchor
        pv_rows.append(
            {
                "row": int(r["source_row"]),
                "original_pressure_gpa": original_p,
                "recalibrated_pressure_gpa": pressure,
                "recalibrated_pressure_error_gpa": sp,
                "calculated_equivalent_quartz_volume_a3": quartz_v,
                "normalized_molar_volume_cm3_mol": observed,
                "calculated_molar_volume_cm3_mol": calc,
                "observed_minus_calculated_cm3_mol": observed - calc,
            }
        )
    result["milani_compression"] = {
        "protocol": "Seven loaded rows; measured ambient row retained as normalization anchor. V=125.35*Vobs/Vambient. Equivalent quartz volumes inferred from the rounded original pressures using the author 1997 quartz EOS, then pressures calculated with Scheidl 2016 BM4. Original paired quartz readings are missing. T=298 K adopted; shared anchor covariance not used for a significance test.",
        "source_s1_discrepancy": "S1 says Milani 2015, 0-8.4 GPa; Milani 2017 deposit Table 3 has eight rows ending at 7.674 GPa, consistent with the 2017 citation in the article/EOS file. No pyrope table substituted.",
        "residual_summary_cm3_mol": summarize(
            [r["observed_minus_calculated_cm3_mol"] for r in pv_rows]
        ),
        "rows": pv_rows,
    }
    moduli = {}
    for name, filename, scale in [
        ("kono", "kono-2010-table-1.csv", p[53]),
        ("isaak_ggd", "isaak-1992-table-3b.csv", 1),
        ("bass", "bass-1989-ambient.csv", 1),
        ("conrad", "conrad-1999-table-6.csv", 1),
    ]:
        rows = []
        for r in observations(filename):
            pressure = float(r.get("pressure_gpa") or 0.0001)
            temperature = float(r.get("temperature_k") or 298)
            if name == "conrad":
                k = (float(r["c11_gpa"]) + 2 * float(r["c12_gpa"])) / 3
                error = None  # C11/C12 covariance was not supplied.
            else:
                k = float(r["ks_gpa"])
                error = r.get("ks_gpa_reported_error", r.get("ks_gpa_error"))
                error = float(error) if error else None
            calc = properties(p, pressure, temperature, GAS_CONSTANT)["ks_gpa"]
            rows.append(
                {
                    "row": int(r["source_row"]),
                    "pressure_gpa": pressure,
                    "temperature_k": temperature,
                    "temperature_adopted": not bool(r.get("temperature_k")),
                    "observed_ks_gpa": k,
                    "reported_ks_error_gpa": error,
                    "calculated_ks_gpa": calc,
                    "calculated_ks_on_observed_scale_gpa": scale * calc,
                    "observed_minus_calculated_gpa": k - scale * calc,
                    "selected": r.get("selection_angel_2022", "selected") == "selected",
                }
            )
        selected = [r for r in rows if r["selected"]]
        moduli[name] = {
            "calculated_modulus_multiplier_to_observed_scale": scale,
            "residual_summary_gpa": summarize(
                [r["observed_minus_calculated_gpa"] for r in selected]
            ),
            "qualification": "Fixed published coefficients and source modulus scale. Ambient T absent from Bass/Conrad tables is adopted as 298 K. No missing modulus error or covariance is synthesized. Conrad KS is derived from measured C11/C12, not an independent bulk-modulus measurement.",
            "rows": rows,
        }
    result["moduli"] = moduli
    rows = []
    for r in observations("isaak-1992-table-2.csv"):
        t = float(r["temperature_k"])
        observed = 1 + 0.01 * float(r["observed_volume_strain_x1e2_token"])
        calc = (
            properties(p, 0.0001, t, GAS_CONSTANT)["volume_cm3_mol"]
            / properties(p, 0.0001, 293, GAS_CONSTANT)["volume_cm3_mol"]
        )
        rows.append(
            {
                "row": int(r["source_row"]),
                "temperature_k": t,
                "gr2b_tabulated_one_plus_volume_strain": observed,
                "source_eos_volume_ratio_to_adopted_293k": calc,
                "difference": observed - calc,
            }
        )
    result["isaak_expansion_identity"] = {
        "reference_temperature_assumption": "293 K adopted for this illustrative ratio comparison from the stated start of the dilatometer runs. The tabulated strain reference is described as ambient, not independently established here as exactly 293 K.",
        "finding": "Primary Table 2 observed strain is GR-2B (76.1% grossular, 22.1% andradite), not the 96.53% grossular GGD specimen listed for this expansion set in S1. Endmember alpha in Table 2 is a calculated compositional estimate, not measured endmember volumes. The original author input is needed to identify which column/reduction entered the 15 expansion channels.",
        "included_in_endmember_refit": False,
        "rows": rows,
    }
    return result


if __name__ == "__main__":
    (DATA / "observation-audit.json").write_text(
        json.dumps(audit(), indent=2, allow_nan=False) + "\n"
    )
