"""The domain review freezes scientific values independently of metadata pins."""

import json
from pathlib import Path

from scripts.scientific_source_snapshot import scientific_sha256

ROOT = Path(__file__).resolve().parents[1]


def test_reviewed_coefficients_units_and_source_arrays_remain_exact():
    review = json.loads(
        (ROOT / "docs/data/pr27-scientific-domain-review.json").read_text(
            encoding="utf-8"
        )
    )
    assert scientific_sha256() == review["catalog_scientific_sha256"]
