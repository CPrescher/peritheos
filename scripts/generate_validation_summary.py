"""Synchronize the validation guide's numerical totals with the refit ledger.

Run ``python -m scripts.generate_validation_summary`` after updating the ledger;
``--check`` reports drift without writing or rerunning scientific fits.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def render_summary(document: str, summary: dict[str, int]) -> str:
    """Update only the five campaign counts, preserving the scientific narrative."""
    replacements = (
        (
            r"(fit for every record with sufficient direct observations and documents all )\d+(\s+records)",
            summary["total"],
        ),
        (r"(The current campaign finds )\d+( parity matches)", summary["parity"]),
        (
            r"(parity matches and )\d+( additional numerically similar results)",
            summary["similar"],
        ),
        (
            r"(\[)\d+( direct refits\]\(primary-eos-refits\.md#parity-not-achieved\))",
            summary["parity_not_achieved"],
        ),
        (r"(while )\d+( records cannot be)", summary["not_refittable"]),
    )
    for pattern, count in replacements:
        document, matches = re.subn(
            pattern, lambda m: m[1] + str(count) + m[2], document
        )
        if matches != 1:
            raise ValueError(f"Expected one validation summary field: {pattern}")
    return document


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    path = ROOT / "docs/validation.md"
    original = path.read_text(encoding="utf-8")
    summary = json.loads(
        (ROOT / "docs/data/primary-eos-refits.json").read_text(encoding="utf-8")
    )["summary"]
    expected = render_summary(original, summary)
    if args.check:
        if original != expected:
            print(
                "Validation totals are stale; run python -m scripts.generate_validation_summary"
            )
            return 1
    elif original != expected:
        if path.read_text(encoding="utf-8") != original:
            raise RuntimeError("Validation guide changed while updating totals")
        path.write_text(expected, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
