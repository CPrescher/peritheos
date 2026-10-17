"""Read-only catalog screening; evaluate existing coefficients, never fit them."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from peritheos import get_eos_record, list_materials

ROOT = Path.cwd()
DATA = ROOT / "peritheos/data"
OUT = ROOT / "output/eos-recommendation-audit-2026-10-17"
paths = sorted((DATA / "materials").glob("*.eosmat"))
before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
documents = [json.loads(p.read_text()) for p in paths]
raw = {r["identifier"]: (d, r) for d in documents for r in d["eos_records"]}
defaults = {m.default_record().identifier for m in list_materials()}
ledger_path = ROOT / "docs/data/primary-eos-refits.json"
ledger = json.loads(ledger_path.read_text())
outcomes = {r["record_identifier"]: r for r in ledger["records"]}
assert set(outcomes) == set(raw)

new = {
    "ca_perovskite_fu_2023_bm3_mgd_refit": (
        "docs/literature-reproductions/fu-2023-eos.md",
        "Completed source/normalization audit: published physical-n=5 model gives 5.917 GPa RMS on 140 Sun pressure observations; the high-temperature discrepancy remains unresolved. This historical _refit identifier denotes Fu's published fit, not the separate Peritheos refit.",
    ),
    "palladium_baty_2024_bm3_1": (
        "docs/literature-reproductions.md#palladium-baty-2024",
        "Experimental BM3: all 78 source rows and equation checked; published curve RMS 1.2731 GPa versus diagnostic 0.7043 GPa. Completed audit excludes transcription, volume basis and EOS formalism; stated fit is not recovered. This does not affect Baty's separate 0 K DFT record.",
    ),
}
for kind in ["type1", "type2", "type4_1", "type4_2", "type4_3"]:
    new[f"fe09ni01_hcp_sakai_2014_{kind}_thermal"] = (
        "docs/literature-reproductions/iron-fei-input-papers.md",
        "Direct replay on all 32 P7 hot rows gives 7.085-7.214 GPa RMS, versus about 1.1 GPa reported in Table 7. Table 8 extrapolation checkpoints agree but do not resolve the hot-data discrepancy. Hold recommendation of the current source implementation until reconciled; no source error or replacement is inferred.",
    )

artifact_groups = [
    (
        [k for k in raw if "jacobsen_2002" in k],
        "Wrong observable: generic audit selected quartz calibrant volume instead of sample lattice/volume. Correct sample-column replay is recorded below.",
        "peritheos/data/datasets/ferropericlase-fe27-jacobsen-2002-table5-compression.csv",
    ),
    (
        ["fe093o_b1_jacobsen_2005_bm3_1", "mg073fe027o_jacobsen_2005_bm3_1"],
        "Wrong material/phase observable: both generic diagnostics selected fe093o_rhombohedral_volume_a3, not the respective cubic sample volume and source fit mask.",
        "peritheos/data/datasets/jacobsen-2005-table1-compression.csv",
    ),
    (
        [
            "mg0215fe0762vac0023o_finkelstein_2017_bm3_helium_1",
            "mg0215fe0762vac0023o_finkelstein_2017_bm3_neon_cubic_2",
        ],
        "Wrong observable and selection: hexagonal_c_over_a was treated as volume and all 29 rows pooled. Correct cubic-column/source-mask replay agrees with the dedicated audit.",
        "docs/literature-reproductions/finkelstein-2017-mw78.md",
    ),
    (
        [k for k in raw if "dobrosavljevic_2019" in k],
        "Both generic fits pool all 28 source rows across phase branches, ignoring explicit used_in_b1_fit/used_in_rhombohedral_fit flags. A phase-specific replay is required before judging coefficient parity.",
        "peritheos/data/datasets/mgfe94o-dobrosavljevic-2019-table-s1-compression.csv",
    ),
    (
        ["e_feooh_hc_low_spin_thompson_2017_bm3_1"],
        "Fit-scope error: all 17 rows include HOC/transition observations. The published HC branch's 12 rows at >=30 GPa give 0.1082 GPa source-curve RMS, not the generic all-row 1.5641 GPa.",
        "peritheos/data/datasets/epsilon-feooh-thompson-2017-supplement-table-s1.csv",
    ),
    (
        ["mg092fe008sio3_bridgmanite_gong_2004_adiabatic_bm3"],
        "Wrong thermodynamic path: the generic audit compares a shock-reduced reference isentrope directly with Hugoniot P-density observations without thermal reduction.",
        "peritheos/data/materials/mg092fe008sio3_bridgmanite.eosmat",
    ),
    (
        ["iron_dewaele_2019_mao_vinet", "iron_dewaele_2019_dor_vinet"],
        "Generic audit pools all 63 upstream rows and does not perform the source ruby-scale re-reduction. The existing dedicated 53-row hcp/two-scale audit reports both fits reproduced.",
        "docs/literature-reproductions/dewaele-2019-static-dac-metals.md",
    ),
    (
        [
            "mg0490fe0483ti0027o_solomatova_2016_hs_bm3_reference_1",
            "mg0490fe0483ti0027o_solomatova_2016_ls_bm3_reference_2",
        ],
        "Wrong model comparison: individual HS/LS reference branches were fitted separately to the mixed-spin dataset. The dedicated coupled Fp48 calculation is classified similar (0.1198 A^3 volume RMS).",
        "docs/literature-reproductions/solomatova-2016-ferropericlase.md",
    ),
]
artifacts = {
    k: (evidence, reason) for ids, reason, evidence in artifact_groups for k in ids
}
assert len(artifacts) == 15

notes = {
    "dorogokupets_2015": (
        "Joint thermoelastic fit is incompletely reconstructed; a 298 K P-V slice alone cannot establish failure of the full source objective.",
        "docs/literature-reproductions/dorogokupets-2015-mgsio3-refit.md",
    ),
    "b4c_somayazulu_2023": (
        "Published heated curve has chi-square/N=0.961. The q-only shift concerns weak identifiability and unspecified objective; it does not establish failed reproduction of the published curve.",
        "docs/literature-reproductions.md#c01-boron-carbide-somayazulu-et-al-2023",
    ),
    "coesite_i_iii_bykova": (
        "Recovered 24-row numerical subset omits source-study fit observations; exact objective/weights are unavailable. Partial refit is not a complete source test.",
        "docs/data/primary-eos-refits.json",
    ),
    "coo_clendenen": (
        "Only a smoothed nine-point table survives, not the original regression observations. Strong coefficient tradeoff and missing objective limit the inference; investigate the curve discrepancy within source tolerances before a new designation.",
        "docs/literature-reproductions.md#coo-clendenen-1966",
    ),
    "diamond_": (
        "Downstream DFT pressure/energy checkpoints or digitized validation markers do not recover the upstream cold-energy/phonon fit objective; individual coefficient non-parity is insufficient.",
        "docs/data/primary-eos-refits.json",
    ),
    "sakai_2014": (
        "Static diagnostic conditions on the published preliminary g-G V0; source preliminary regression/weights/covariance are not recovered. Do not equate this with the separately confirmed five thermal-curve discrepancies.",
        "docs/literature-reproductions/iron-fei-input-papers.md",
    ),
    "nickel_oxide_b1_campbell": (
        "NiO's current-table diagnostic does not by itself establish failure of the complete coupled/source fit. Review separately from the already deferred fcc-Fe/FeO records.",
        "docs/literature-reproductions/campbell-2009-buffers.md",
    ),
    "dorfman_2012": (
        "Available-row pair inventory is not confirmed as the author inventory; fixed-metal results lie within the authors' approximately 2% realistic uncertainty. Published table pressures reproduce within rounding. Exact free-fit coefficient parity remains unresolved, not a blanket six-record rejection.",
        "docs/literature-reproductions/dorfman-2012-cocompression.md",
    ),
    "lead_hcp_dewaele": (
        "Generic fit uses a partial digitized upstream subset without the complete source row selection/calibration handling; dedicated audit identifies missing older Pb inputs.",
        "docs/literature-reproductions/dewaele-2019-static-dac-metals.md",
    ),
    "iron_yamazaki": (
        "Published source checkpoints pass and pressure RMS changes only 0.661/0.690 to 0.622/0.620 GPa; differing conditional coefficients with unknown weights/covariance alone are insufficient for rejection.",
        "docs/literature-reproductions/iron-fei-input-papers.md",
    ),
    "mg075fe025o_solomatova": (
        "Coupled crossover diagnostic remains non-parity; standalone branch records are components of that coupled source representation. Review the coupled model and source priors before deciding, rather than reject a reference branch on a standalone mixed-spin fit.",
        "docs/literature-reproductions/solomatova-2016-ferropericlase.md",
    ),
    "mg080fe020o_fei": (
        "Low-spin validation differs, but exact source weights and row mask are unavailable and the selection endpoint is approximate. Existing consult-audit qualification should remain; no automatic new designation from the coefficient flag alone.",
        "docs/literature-reproductions/fei-2007-ferropericlase.md",
    ),
    "nacl_b2_fei": (
        "Published calculated curves agree within source line width. q is highly sensitive to digitization, with no measurement-error model; its refit failure alone does not justify rejecting the source curve.",
        "docs/literature-reproductions/fei-2007-nacl-b2.md",
    ),
    "platinum_fei_2004": (
        "Diagnostic fixes source choices that are not fully specified and does not reproduce the complete upstream joint calibration objective. Do not reject solely on its conditional coefficient comparison.",
        "docs/data/primary-eos-refits.json",
    ),
    "platinum_fei_2007": (
        "Thermal curve spacing agrees within source readout precision, whereas the diagnostic coefficients depart from the plotted branches. Exact fitted/fixed list and weights remain unknown. Existing investigation is paused; no new fitting was performed.",
        "docs/literature-reproductions/fei-2007-platinum.md",
    ),
    "silicon_vii": (
        "High-pressure-only data do not identify an unconstrained ambient EOS; refit hits a bound. Recover source constraints/selection before interpreting coefficient non-parity.",
        "docs/data/primary-eos-refits.json",
    ),
    "sno2_": (
        "Both aliases use the same high-pressure-only dataset and share extrapolation/non-identifiability. Bound-hitting refit is not two independent source failures.",
        "docs/data/primary-eos-refits.json",
    ),
}

screened = []
for key, row in outcomes.items():
    if row["status"] != "parity_not_achieved":
        continue
    doc, meta = raw[key]
    if key in new:
        decision, (evidence, reason) = "recommend_not_recommended", new[key]
    elif (
        meta["scientific_validation"]["status"] != "primary_source_validated"
        or key == "goethite_gleason_2008_bm3_1"
    ):
        decision = "already_flagged_or_held"
        evidence, reason = (
            "docs/data/primary-eos-refits.json",
            "Existing not_reproduced/deferred status or explicit not-recommended designation. Preserve the restriction; standardize recommendation metadata where absent.",
        )
    elif key in artifacts:
        decision, (evidence, reason) = "audit_comparison_defect", artifacts[key]
    elif key == "wadsleyite_katsura_2009_bm3_1":
        decision = "hold_current_implementation_pending_mapping_review"
        evidence = "scripts/validate_primary_eos_refits.py"
        reason = "Direct replay on all 85 normalized-volume rows confirms 1.8944 GPa pressure RMS and 3.4638 GPa maximum. Catalog default. Existing audit still identifies unresolved calibration and possible MGD reference/Debye convention mismatch. Do not recommend the current implementation pending that check; do not attribute the mismatch to a source error."
    else:
        matches = [
            (reason, evidence)
            for token, (reason, evidence) in notes.items()
            if token in key
        ]
        assert len(matches) == 1, (key, matches)
        reason, evidence = matches[0]
        decision = "insufficient_for_automatic_downgrade"
    screened.append(
        {
            "identifier": key,
            "material": doc["identifier"],
            "is_material_default": key in defaults,
            "scientific_validation": meta["scientific_validation"]["status"],
            "current_usage_recommendation": meta["scientific_validation"].get(
                "usage_recommendation"
            ),
            "ledger_status": row["status"],
            "decision": decision,
            "reason": reason,
            "evidence": evidence,
            "ledger_pressure_rmse_gpa": row.get("published_rmse_gpa"),
            "ledger_diagnostic_columns": row.get("columns"),
            "ledger_diagnostic_selection": row.get("selection"),
        }
    )
assert len(screened) == 75
assert Counter(r["decision"] for r in screened) == {
    "recommend_not_recommended": 7,
    "already_flagged_or_held": 6,
    "audit_comparison_defect": 15,
    "hold_current_implementation_pending_mapping_review": 1,
    "insufficient_for_automatic_downgrade": 46,
}

already = []
for key, (doc, meta) in raw.items():
    status = meta["scientific_validation"]["status"]
    if status != "primary_source_validated" or key == "goethite_gleason_2008_bm3_1":
        already.append(
            {
                "identifier": key,
                "material": doc["identifier"],
                "status": status,
                "current_usage_recommendation": meta["scientific_validation"].get(
                    "usage_recommendation"
                ),
                "is_material_default": key in defaults,
            }
        )
assert len(already) == 11

checks = []


def check(key, column, mask=None, cube=False, scale=1):
    doc, meta = raw[key]
    dataset = next(
        d
        for d in doc["datasets"]
        if d["identifier"] in meta.get("fit_datasets", [])
        or ("katsura_2009" in d["identifier"] and "wadsleyite" in key)
    )
    rows = list(csv.DictReader((DATA / dataset["resource"]["path"]).open()))
    selected = [r for r in rows if r.get(column) and (mask is None or mask(r))]
    model = get_eos_record(key)
    residuals = []
    for row in selected:
        volume = float(row[column]) ** (3 if cube else 1) * scale
        temp = (
            float(row["temperature_k"])
            if model.is_thermal
            else model.reference_temperature
        )
        residuals.append(
            float(model.pressure(volume, temp)) - float(row["pressure_gpa"])
        )
    checks.append(
        {
            "identifier": key,
            "observations": len(selected),
            "volume_column": column,
            "cube_lattice_parameter": cube,
            "volume_scale": scale,
            "pressure_rmse_gpa": float(np.sqrt(np.mean(np.array(residuals) ** 2))),
            "maximum_absolute_pressure_residual_gpa": max(abs(x) for x in residuals),
        }
    )


for prefix in ["ferropericlase_fe27", "magnesiowustite_fe56", "magnesiowustite_fe75"]:
    check(
        next(k for k in raw if k.startswith(prefix) and "jacobsen_2002" in k),
        "lattice_a_angstrom",
        cube=True,
    )
for key, prefix in [
    ("mg073fe027o_jacobsen_2005_bm3_1", "mg073fe027o"),
    ("fe093o_b1_jacobsen_2005_bm3_1", "fe093o_b1"),
]:
    check(
        key,
        prefix + "_derived_volume_a3",
        lambda r, prefix=prefix: r[prefix + "_fit_included"].lower() in ["1", "true"],
    )
for suffix, flag in [
    ("helium_1", "used_in_helium_cubic_fit"),
    ("neon_cubic_2", "used_in_neon_cubic_fit"),
]:
    check(
        "mg0215fe0762vac0023o_finkelstein_2017_bm3_" + suffix,
        "cubic_volume_a3",
        lambda r, flag=flag: r[flag].lower() in ["1", "true"],
    )
for prefix, flag, scale in [
    ("mgfe94o_b1", "used_in_b1_fit", 8),
    ("mgfe94o_rhombohedral", "used_in_rhombohedral_fit", 2),
]:
    check(
        next(k for k in raw if k.startswith(prefix)),
        "atomic_volume_angstrom3",
        lambda r, flag=flag: r[flag].lower() in ["1", "true"],
        scale=scale,
    )
check(
    "e_feooh_hc_low_spin_thompson_2017_bm3_1",
    "unit_cell_volume_a3",
    lambda r: float(r["pressure_gpa"]) >= 30,
)
check(
    "wadsleyite_katsura_2009_bm3_1",
    "normalized_volume",
    scale=raw["wadsleyite_katsura_2009_bm3_1"][1]["eos"]["parameters"]["V0"],
)

doc, _ = raw["fe09ni01_hcp_sakai_2014_type1_thermal"]
dataset = next(d for d in doc["datasets"] if d["identifier"] == "sakai_2014_table2_pvt")
hot = [
    r
    for r in csv.DictReader((DATA / dataset["resource"]["path"]).open())
    if r.get("pressure_p7_gpa") and float(r["temperature_k"]) > 300
]
assert len(hot) == 32 and {r["material"] for r in hot} == {"fe09ni01_hcp"}
for key in new:
    if "sakai_2014" not in key:
        continue
    model = get_eos_record(key)
    residuals = [
        float(
            model.pressure(
                float(r["volume_a3_conventional_cell"]), float(r["temperature_k"])
            )
        )
        - float(r["pressure_p7_gpa"])
        for r in hot
    ]
    checks.append(
        {
            "identifier": key,
            "observations": 32,
            "pressure_column": "pressure_p7_gpa",
            "pressure_rmse_gpa": float(np.sqrt(np.mean(np.array(residuals) ** 2))),
            "maximum_absolute_pressure_residual_gpa": max(abs(x) for x in residuals),
        }
    )

report = {
    "audit_date": "2026-10-17",
    "scope": "Catalog-wide screening against existing curation evidence; no new fitting, catalog edits, source correction or author contact. This is not an independent revalidation of every publication.",
    "material_documents": len(documents),
    "stored_records": len(raw),
    "validation_counts": dict(
        Counter(r["scientific_validation"]["status"] for _, r in raw.values())
    ),
    "ledger_sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
    "parity_failure_screen": sorted(screened, key=lambda r: r["identifier"]),
    "screen_decision_counts": dict(Counter(r["decision"] for r in screened)),
    "already_flagged_or_held": sorted(already, key=lambda r: r["identifier"]),
    "direct_evaluations_no_refitting": checks,
    "additional_qualification_review": [
        {
            "identifier": "argon_hcp_wittlinger_1997_bm2_1",
            "evidence": "docs/literature-reproductions/argon-fcc.md#wittlinger-1997-published-fit-not-reproduced",
            "reason": "Material default, primary_source_validated with reproduction_status=not_reproduced. Source-data/figure-only limitations need explicit recommendation semantics; missing full raw data alone is not a demonstrated source-curve failure.",
        },
        {
            "identifier": "argon_fcc_errandonea_2006_bm3",
            "evidence": "docs/literature-reproductions/argon-errandonea-2006.md",
            "reason": "Primary_source_validated with reproduction_status=not_reproduced; eight recovered markers and bound-dependent refit cannot identify the full author fit. Qualify evidence without equating missing complete inputs with a demonstrated numerical contradiction.",
        },
    ],
    "catalog_fingerprint": hashlib.sha256(
        json.dumps(before, sort_keys=True).encode()
    ).hexdigest(),
}
assert before == {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "audit.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
print(
    json.dumps(
        {
            "stored_records": len(raw),
            "screen": report["screen_decision_counts"],
            "already_flagged_or_held": len(already),
            "direct_evaluations": len(checks),
            "catalog_unchanged": True,
        }
    )
)
