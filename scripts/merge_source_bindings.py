"""Repository-only compatibility pins for the reviewed origin/main integration.

The original release and cohort fingerprints remain immutable. The supplemental
file records both parent revisions and exact merged fingerprints; it does not
replace source coefficients, resource checksums, or historical audit inputs.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROOF = ROOT / "docs/data/merged-source-bindings.json"


@lru_cache(maxsize=1)
def merge_bindings():
    return json.loads(PROOF.read_text(encoding="utf-8"))["bindings"]


def accepts_merge_binding(
    context: str, identifier: str, expected: str, actual: str
) -> bool:
    """Accept only the explicitly reviewed original-to-merged fingerprint pair."""
    entry = merge_bindings().get(context, {}).get(identifier, {})
    return (
        entry.get("historical_sha256") == expected
        and entry.get("merged_sha256") == actual
    )
