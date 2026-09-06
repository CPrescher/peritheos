"""Regression checks for the final citation/compilation source audit."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDITED_DOIS = {
    "10.1029/gl017i008p01153",
    "10.1029/jb095ib12p19311",
    "10.1029/92gl02960",
    "10.1029/93gl01265",
    "10.1126/science.259.5091.66",
    "10.1029/94eo01093",
    "10.1029/96gl03951",
    "10.2138/am-2000-1013",
    "10.1046/j.1365-246x.2001.00437.x",
    "10.1103/physrevb.65.104114",
    "10.1002/pssb.200302047",
    "10.1063/1.1780510",
    "10.1088/0953-8984/16/30/006",
    "10.1103/physrevlett.93.215502",
    "10.1073/pnas.0608609104",
    "10.1142/s0217979208038910",
    "10.12693/aphyspola.115.709",
    "10.2138/am.2010.3368",
    "10.1142/s0217984912501461",
    "10.1088/1742-6596/653/1/012095",
    "10.1063/1.4967779",
    "10.1139/cjp-2019-0326",
    "10.1186/s40645-020-00333-3",
    "10.1016/j.fmre.2021.12.013",
    "10.69626/sea.2024.0152",
    "10.1038/s43247-025-02383-1",
    "10.15407/mfint.47.06.0601",
    "10.3724/j.issn.1007-2802.20240151",
}


def _load(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def test_every_candidate_is_explicitly_in_the_audit() -> None:
    candidates = _load("docs/data/litcurate-eos-candidates.json")["records"]
    selected = [
        row
        for row in candidates
        if row["publication"]["doi"].casefold() in AUDITED_DOIS
    ]
    audit = (
        ROOT / "docs/literature-reproductions/"
        "litcurate-source-exhaustion-citation-audit.md"
    ).read_text(encoding="utf-8")

    assert len(selected) == 122
    assert {row["publication"]["doi"].casefold() for row in selected} == AUDITED_DOIS
    assert all(f"`{row['identifier']}`" in audit for row in selected)


def test_every_extracting_paper_is_in_nonproduction_registry() -> None:
    papers = _load("docs/data/nonproduction-paper-investigations.json")["papers"]
    registered = {paper["doi"].casefold() for paper in papers if paper.get("doi")}

    assert AUDITED_DOIS <= registered


def test_every_litcurate_publication_has_a_paper_level_disposition() -> None:
    candidates = _load("docs/data/litcurate-eos-candidates.json")["records"]
    source_audits = _load("peritheos/data/primary-source-audit.json")["records"]
    nonproduction = _load("docs/data/nonproduction-paper-investigations.json")["papers"]

    candidate_dois = {row["publication"]["doi"].casefold() for row in candidates}
    production_dois = {row["doi"].casefold() for row in source_audits if row.get("doi")}
    nonproduction_dois = {
        row["doi"].casefold() for row in nonproduction if row.get("doi")
    }

    assert candidate_dois <= production_dois | nonproduction_dois
