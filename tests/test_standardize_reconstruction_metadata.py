"""The reconstruction-metadata drift check on bundled (record + sidecar) files."""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

from scripts import standardize_reconstruction_metadata as tool

ROOT = Path(__file__).resolve().parents[1]
MATERIAL = "siderite_fe095mn005"
RECORD = "siderite_fe095mn005_litasov_2013_bm3"


def _run(monkeypatch, *arguments) -> int:
    monkeypatch.setattr(sys, "argv", ["standardize", *map(str, arguments)])
    return tool.main()


def _copy(tmp_path: Path) -> tuple[Path, Path]:
    record = tmp_path / "materials" / f"{MATERIAL}.eosmat"
    sidecar = tmp_path / "audits" / f"{MATERIAL}.audit.json"
    for target, folder in ((record, "materials"), (sidecar, "audits")):
        target.parent.mkdir()
        shutil.copy(ROOT / "peritheos/data" / folder / target.name, target)
    return record, sidecar


def test_bundled_catalog_has_no_drift(monkeypatch, capsys):
    assert _run(monkeypatch, "--check") == 0
    assert "0 changed" in capsys.readouterr().out


def test_a_missing_sidecar_assessment_is_reported_and_never_written(
    tmp_path, monkeypatch, capsys
):
    record, sidecar = _copy(tmp_path)
    assert _run(monkeypatch, "--check", record) == 0
    document = json.loads(sidecar.read_bytes())
    del document["records"][RECORD]["evidence"]["reproduction"]["assessment"]
    sidecar.write_text(json.dumps(document), encoding="utf-8")
    before = record.read_bytes(), sidecar.read_bytes()
    capsys.readouterr()
    assert _run(monkeypatch, "--check", record) == 1
    assert "Needs reconstruction metadata" in capsys.readouterr().out
    with pytest.raises(SystemExit, match="apply_catalog_revision"):
        _run(monkeypatch, record)
    assert (record.read_bytes(), sidecar.read_bytes()) == before
