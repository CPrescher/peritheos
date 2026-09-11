import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = ROOT / "peritheos" / "data" / "primary-source-audit.json"
NONPRODUCTION_PATH = ROOT / "docs" / "data" / "nonproduction-paper-investigations.json"
LEDGER_PATH = ROOT / "docs" / "paper-investigation-ledger.md"


def test_paper_investigation_ledger_is_complete_and_current():
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))["records"]
    nonproduction = json.loads(NONPRODUCTION_PATH.read_text(encoding="utf-8"))["papers"]
    paper_keys = {
        ("doi", row["doi"].lower())
        if row.get("doi")
        else ("url", row["primary_source_check"]["access_url"])
        for row in audit
    }

    assert len(paper_keys) == 235
    assert len(nonproduction) >= 41
    assert {row["outcome"] for row in nonproduction} == {
        "withheld_unreproduced",
        "deferred_incomplete_model",
        "direct_refit_unavailable",
    }
    subprocess.run(
        [sys.executable, "scripts/generate_paper_investigation_ledger.py", "--check"],
        cwd=ROOT,
        check=True,
    )
    ledger = LEDGER_PATH.read_text(encoding="utf-8")
    assert f"**{len(paper_keys) + len(nonproduction)} primary papers**" in ledger
    assert "| Reproduced |" in ledger
    assert "| Coefficient parity not achieved |" in ledger
    assert "| Direct refit unavailable |" in ledger
    assert "| Source reconstruction | 1 |" in ledger
    assert "Katsura et al. (2004)" in ledger
    assert "Wang et al. (2026)" in ledger
    assert (
        "Zhang et al. (2025)</a>" not in ledger
        and "3 final-input parity, upstream reduction partial" in ledger
    )


def test_source_reconstruction_does_not_count_as_an_independent_paper_refit():
    from scripts.generate_paper_investigation_ledger import classify

    assert classify(Counter(source_reconstruction=2)) == "source_reconstruction"
    assert (
        classify(Counter(parity=1, source_reconstruction=2)) == "partial_reproduction"
    )
    assert classify(Counter(not_refittable=1)) == "direct_refit_unavailable"
