import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = ROOT / "peritheos" / "data" / "primary-source-audit.json"
NONPRODUCTION_PATH = (
    ROOT / "docs" / "data" / "nonproduction-paper-investigations.json"
)
LEDGER_PATH = ROOT / "docs" / "paper-investigation-ledger.md"


def test_paper_investigation_ledger_is_complete_and_current():
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))["records"]
    nonproduction = json.loads(NONPRODUCTION_PATH.read_text(encoding="utf-8"))[
        "papers"
    ]
    paper_keys = {
        ("doi", row["doi"].lower())
        if row.get("doi")
        else ("url", row["primary_source_check"]["access_url"])
        for row in audit
    }

    assert len(paper_keys) == 144
    assert len(nonproduction) == 2
    assert {row["outcome"] for row in nonproduction} == {
        "withheld_unreproduced",
        "deferred_incomplete_model",
    }
    subprocess.run(
        [sys.executable, "scripts/generate_paper_investigation_ledger.py", "--check"],
        cwd=ROOT,
        check=True,
    )
    ledger = LEDGER_PATH.read_text(encoding="utf-8")
    assert "**146 primary papers**" in ledger
    assert "| Reproduced | 106 |" in ledger
    assert "| Coefficient parity not achieved | 9 |" in ledger
    assert "| Direct refit unavailable | 23 |" in ledger
    assert "Katsura et al. (2004)" in ledger
    assert "Wang et al. (2026)" in ledger
