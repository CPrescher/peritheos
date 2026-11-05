"""Repository-only compatibility pins for reviewed source-metadata changes.

The original release and cohort fingerprints remain immutable. The supplemental
file records both parent revisions and exact merged fingerprints; it does not
replace source coefficients, resource checksums, or historical audit inputs.
The separate licensing ledger chains exact license-only updates to those
existing predecessors; tests recover the previous complete declarations.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROOF = ROOT / "docs/data/merged-source-bindings.json"
LICENSE_PROOF = ROOT / "docs/data/license-metadata-source-bindings.json"


@lru_cache(maxsize=1)
def merge_bindings():
    return json.loads(PROOF.read_text(encoding="utf-8"))["bindings"]


@lru_cache(maxsize=1)
def license_metadata_bindings():
    return json.loads(LICENSE_PROOF.read_text(encoding="utf-8"))["bindings"]


def accepts_merge_binding(
    context: str, identifier: str, expected: str, actual: str
) -> bool:
    """Accept exact reviewed pairs, retaining the historical/predecessor chain."""
    entry = merge_bindings().get(context, {}).get(identifier, {})
    if entry and (
        entry.get("historical_sha256") == expected
        and entry.get("merged_sha256") == actual
    ):
        return True
    review = license_metadata_bindings().get(context, {}).get(identifier, {})
    if bool(review) and (
        review.get("historical_sha256") == expected
        and review.get("licensed_sha256") == actual
        and (not entry or entry.get("historical_sha256") == expected)
        and review.get("previous_reviewed_sha256")
        == (entry["merged_sha256"] if entry else expected)
    ):
        return True
    # A separately authorized structural substitution is not a license-only
    # edit. Keep its exact pair in its own receipt and the original pins intact.
    replacement = json.loads(
        (ROOT / "docs/data/brucite-structure-replacement.json").read_bytes()
    )
    structure = replacement["bindings"].get(context, {}).get(identifier, {})
    return bool(structure) and (
        structure["historical_sha256"] == expected
        and structure["replacement_sha256"] == actual
    )
