"""Research commands must round-trip UTF-8 evidence on an ANSI Windows locale."""

import copy
import importlib
import json
import sys
from pathlib import Path

import pytest

from peritheos import load_eosmat


@pytest.fixture
def ansi_path_defaults(monkeypatch):
    read_text, write_text = Path.read_text, Path.write_text

    def read(path, encoding=None, errors=None, **kwargs):
        return read_text(path, encoding=encoding or "cp1252", errors=errors, **kwargs)

    def write(path, data, encoding=None, errors=None, **kwargs):
        return write_text(
            path, data, encoding=encoding or "cp1252", errors=errors, **kwargs
        )

    monkeypatch.setattr(Path, "read_text", read)
    monkeypatch.setattr(Path, "write_text", write)


def test_integration_preserves_utf8_catalog_and_source_values(
    tmp_path, monkeypatch, ansi_path_defaults
):
    from scripts import review_similar_parameter_scope as audit

    original_root = audit.ROOT
    material = "argon_fcc"
    identifier = "argon_fcc_chen_2010_bm3_digitized_refit"
    source = original_root / f"peritheos/data/materials/{material}.eosmat"
    before_bytes = source.read_bytes()
    assert "—" in before_bytes.decode("utf-8")
    destination = tmp_path / source.relative_to(original_root)
    destination.parent.mkdir(parents=True)
    destination.write_bytes(before_bytes)
    before = load_eosmat(destination)
    ledger = json.loads(
        (original_root / "docs/data/primary-eos-refits.json").read_text(
            encoding="utf-8"
        )
    )
    ledger["records"] = [
        r for r in ledger["records"] if r["record_identifier"] == identifier
    ]
    ledger_path = tmp_path / "docs/data/primary-eos-refits.json"
    ledger_path.parent.mkdir(parents=True)
    ledger_path.write_bytes(json.dumps(ledger).encode("utf-8"))
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    monkeypatch.setattr(
        audit, "saved_review", lambda: {identifier: {"status": "parity"}}
    )
    # Isolate the I/O path: calculations and source verification have their own
    # independent tests; use the real native document and assessment formatter.
    monkeypatch.setattr(
        audit, "augment_outcome", lambda record, row: copy.deepcopy(row)
    )
    monkeypatch.setattr(audit, "verify_sources", lambda: None)
    audit.integrate()
    after = load_eosmat(destination)
    for document in (before, after):
        for record in document["eos_records"]:
            record["evidence"]["reproduction"].pop("assessment", None)
    assert before == after
    assert "—" in destination.read_bytes().decode("utf-8")
    json.loads(ledger_path.read_bytes().decode("utf-8"))
    (tmp_path / "docs/primary-eos-refits.md").read_bytes().decode("utf-8")
    assert source.read_bytes() == before_bytes


def test_scope_report_write_then_check_uses_utf8(
    tmp_path, monkeypatch, ansi_path_defaults, capsys
):
    from scripts import review_similar_parameter_scope as audit

    report = {"records": [], "note": "Measured α — source unchanged"}
    rendered = "# Source review\n\nMeasured α — source unchanged.\n"
    monkeypatch.setattr(audit, "OUTPUT", tmp_path / "report.json")
    monkeypatch.setattr(audit, "REPORT", tmp_path / "report.md")
    monkeypatch.setattr(audit, "calculate", lambda: report)
    monkeypatch.setattr(audit, "render_report", lambda result: rendered)

    def check_report(left, right):
        assert left == right

    monkeypatch.setattr(audit, "check_report", check_report)
    monkeypatch.setattr(sys, "argv", ["review"])
    audit.main()
    assert audit.REPORT.read_bytes() == rendered.encode("utf-8")
    assert json.loads(audit.OUTPUT.read_bytes().decode("utf-8")) == report
    monkeypatch.setattr(sys, "argv", ["review", "--check"])
    audit.main()
    capsys.readouterr()


@pytest.mark.parametrize(
    "module_name, outputs",
    [
        ("scripts.plot_sata_2010_confidence", ("REPORT", "ELLIPSES", "BANDS")),
        (
            "scripts.refit_sata_2010_recalculated_pressure",
            ("REPORT", "INPUTS", "RESIDUALS"),
        ),
    ],
)
def test_diagnostic_write_then_check_uses_utf8(
    module_name, outputs, tmp_path, monkeypatch, ansi_path_defaults, capsys
):
    audit = importlib.import_module(module_name)
    report = {"inputs": [{"label": "α —"}], "residuals": [{"label": "β —"}]}
    for name in outputs:
        monkeypatch.setattr(audit, name, tmp_path / (name.lower() + ".txt"))
    if outputs[1] == "ELLIPSES":
        monkeypatch.setattr(
            audit, "build", lambda: (report, report["inputs"], report["residuals"])
        )
    else:
        monkeypatch.setattr(audit, "reproduce", lambda: report)
    monkeypatch.setattr(sys, "argv", ["audit", "--write"])
    audit.main()
    for name in outputs[1:]:
        text = getattr(audit, name).read_bytes().decode("utf-8")
        assert "—" in text and "\r\n" not in text
    monkeypatch.setattr(sys, "argv", ["audit", "--check"])
    audit.main()
    capsys.readouterr()
