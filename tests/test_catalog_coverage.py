"""Keep the published coverage overview synchronized with accepted material files."""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_generated_catalog_coverage_is_current_and_all_records_have_ledger_entries():
    import json

    subprocess.run(
        [sys.executable, "-m", "scripts.generate_catalog_coverage", "--check"],
        cwd=ROOT,
        check=True,
    )
    coverage = json.loads(
        (ROOT / "docs/data/catalog-coverage.json").read_text(encoding="utf-8")
    )
    assert coverage["ledgers"]["missing_from_refit_ledger"] == []
    assert coverage["ledgers"]["missing_from_source_audit"] == []
    assert coverage["ledgers"]["absent_from_catalog_refit_entries"] == []
    assert coverage["ledgers"]["absent_from_catalog_source_audit_entries"] == []
