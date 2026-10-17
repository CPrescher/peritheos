"""Apply the reviewed October 17 recommendation holds without changing source values."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
AUDIT = "docs/research/eos-recommendation-dispositions-2026-10-17.md"
NOT_RECOMMENDED = {
    "ca_perovskite_fu_2023_bm3_mgd_refit": "The published thermal parameterization retains an unresolved high-temperature discrepancy against the reconstructed source observations; source normalization checks do not reconcile it.",
    "palladium_baty_2024_bm3_1": "The verified experimental Table S1 observations and source equation do not recover the published experimental BM3 fit under the completed audit checks.",
    **{
        f"fe09ni01_hcp_sakai_2014_{kind}_thermal": "The 32 P7 hot observations give about 7.1 GPa pressure RMS with the published coefficients, versus about 1.1 GPa reported in Table 7; passing Table 8 checkpoints does not resolve that discrepancy."
        for kind in ("type1", "type2", "type4_1", "type4_2", "type4_3")
    },
}
WADSLEYITE = "wadsleyite_katsura_2009_bm3_1"


def apply() -> None:
    found = set()
    for path in sorted(MATERIALS.glob("*.eosmat")):
        document = json.loads(path.read_text())
        changed = False
        for record in document["eos_records"]:
            identifier = record["identifier"]
            if identifier not in NOT_RECOMMENDED and identifier != WADSLEYITE:
                continue
            found.add(identifier)
            validation = record["scientific_validation"]
            if identifier == WADSLEYITE:
                record.pop("catalog_access", None)
                validation.update(
                    status="deferred",
                    usage_recommendation="on_hold_not_recommended_pending_source_mapping_review",
                    recommendation_audit=AUDIT,
                    recommendation_date="2026-10-17",
                    unresolved=[
                        "Pressure calibration and MGD reference-state/Debye convention mapping remain unresolved.",
                        "Direct replay of all 85 source rows gives 1.89443 GPa pressure RMS and 3.46380 GPa maximum residual.",
                    ],
                )
                decision = "Current wadsleyite default on hold and not recommended pending source/model mapping review. This does not establish an error in the publication; no replacement default is selected."
                document["recommendation_hold"] = {
                    "record": identifier,
                    "date": "2026-10-17",
                    "status": "on_hold",
                    "reason": decision,
                    "audit": AUDIT,
                }
            else:
                record["catalog_access"] = "explicit_selection"
                record["default"] = False
                validation.update(
                    status="not_reproduced",
                    reproduction_status="not_reproduced",
                    usage_recommendation="not_recommended_for_quantitative_use",
                    recommendation_date="2026-10-17",
                    recommendation_audit=AUDIT,
                )
                decision = (
                    "Not recommended for scientific predictions or pressure calibration. Retained only for explicit diagnostic inspection. "
                    + NOT_RECOMMENDED[identifier]
                )
            validation["recommendation_reason"] = decision
            # Keep the original source audit and its date distinct from this
            # recommendation decision; no coefficient or observation is edited.
            changed = True
        if changed:
            path.write_text(json.dumps(document, indent=1, ensure_ascii=False) + "\n")
    expected = set(NOT_RECOMMENDED) | {WADSLEYITE}
    if found != expected:
        raise ValueError(f"Missing reviewed records: {expected - found}")


if __name__ == "__main__":
    apply()
