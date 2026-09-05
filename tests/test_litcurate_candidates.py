import importlib.util
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "import_litcurate_eos_candidates.py"
LEDGER_PATH = ROOT / "docs" / "data" / "litcurate-eos-candidates.json"

SPEC = importlib.util.spec_from_file_location("litcurate_candidates", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
litcurate_candidates = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(litcurate_candidates)


def test_published_litcurate_candidate_ledger_is_internally_consistent():
    ledger = json.loads(LEDGER_PATH.read_text(encoding="utf-8"))

    assert ledger["format"] == "peritheos.literature-eos-candidates"
    assert ledger["format_version"] == 1
    assert ledger["source"]["doi"] == "10.5281/zenodo.22118629"
    assert ledger["source"]["license"] == "CC-BY-4.0"
    assert len(ledger["source"]["input_sha256"]) == 64
    assert ledger["summary"] == {
        "records": 1309,
        "unique_publication_dois": 204,
        "source_reported_publication_dois": 168,
        "source_reported_publication_dois_not_bundled": 160,
        "by_origin": {
            "Citation-reported": 727,
            "Source-reported": 578,
            "Unspecified": 4,
        },
        "by_review_bucket": {
            "already_bundled_source": 39,
            "citation_trace": 727,
            "equation_audit": 115,
            "manual_triage": 4,
            "model_work": 23,
            "not_an_eos_fit": 17,
            "primary_source_audit": 110,
            "review_first": 274,
        },
    }

    records = ledger["records"]
    identifiers = [record["identifier"] for record in records]
    assert len(identifiers) == len(set(identifiers)) == 1309
    assert all(
        record["publication"]["doi"] != "10.2138/am-2024-9562"
        for record in records
    )
    assert all("evidence" not in record["litcurate"] for record in records)

    possible_duplicate_dois = {
        record["publication"]["doi"]
        for record in records
        if "possible_duplicate_publication_identifier"
        in record["classification"]["flags"]
    }
    assert {
        "10.2138/am-2018-6694",
        "10.2138/am-2019-6694",
    } <= possible_duplicate_dois


def test_review_classification_keeps_discovery_separate_from_execution(tmp_path):
    source = tmp_path / "candidates.csv"
    source.write_text("source fixture\n", encoding="utf-8")
    rows = [
        {
            "doi": "https://doi.org/10.1000/new",
            "title": "New EOS",
            "year": 2026,
            "phase": "Example",
            "composition": "X",
            "eos_model": "third-order Birch-Murnaghan",
            "origin": "Source-reported",
            "V0": 10,
            "K0": 100,
            "Kp": 4,
            "evidence": "Table 1",
        },
        {
            "doi": "10.1000/cited",
            "title": "Compilation",
            "year": 2025,
            "phase": "Example",
            "eos_model": "third-order Birch-Murnaghan",
            "origin": "Citation-reported",
            "V0": 10,
            "K0": 100,
            "Kp": 4,
            "evidence": "Table 2",
        },
    ]

    ledger = litcurate_candidates.build_ledger(
        rows,
        source_path=source,
        material_dois=set(),
        candidate_backlog_dois={"10.1000/new"},
    )

    first, second = ledger["records"]
    assert first["classification"]["review_bucket"] == "review_first"
    assert first["classification"]["suggested_peritheos_mapping"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
    }
    assert first["classification"]["documented_backlog_doi"] is True
    assert second["classification"]["review_bucket"] == "citation_trace"
    assert "eos_records" not in ledger
    assert "scientific_validation" not in first


def test_native_litcurate_json_export_is_flattened_and_origin_is_normalized(tmp_path):
    source = tmp_path / "database.json"
    source.write_text(
        json.dumps(
            {
                "papers": [
                    {
                        "paper_id": "10.1000_example",
                        "source": {
                            "payload": {
                                "doi": "10.1000/example",
                                "title": "Example EOS",
                                "year": 2026,
                            }
                        },
                        "peritheos_eos_candidate": {
                            "payload": {
                                "eos_entries": [
                                    {
                                        "material": {
                                            "phase": "Example",
                                            "composition": "X",
                                        },
                                        "equation": {
                                            "reported_name": "Vinet",
                                            "family": "Vinet",
                                        },
                                        "origin": "source_reported",
                                        "parameters": [
                                            {
                                                "name": "V0",
                                                "reported_value": 10,
                                                "unit": "angstrom^3",
                                                "determination": "fitted",
                                            },
                                            {
                                                "name": "K0",
                                                "reported_value": 100,
                                                "unit": "GPa",
                                                "determination": "fitted",
                                            },
                                            {
                                                "name": "K0'",
                                                "reported_value": 4,
                                                "unit": None,
                                                "determination": "fitted",
                                            },
                                        ],
                                        "evidence_locations": [
                                            {
                                                "field": "parameters",
                                                "locator": "Table 1",
                                                "excerpt": "V0=10, K0=100, K0'=4",
                                            }
                                        ],
                                    }
                                ]
                            }
                        },
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    rows = litcurate_candidates.load_json(source)
    ledger = litcurate_candidates.build_ledger(
        rows,
        source_path=source,
        material_dois=set(),
        candidate_backlog_dois=set(),
    )

    assert ledger["records"][0]["reported_eos"]["origin"] == "Source-reported"
    assert ledger["records"][0]["classification"]["review_bucket"] == "review_first"
    assert len(ledger["records"][0]["reported_eos"]["all_parameters"]) == 3
    assert (
        ledger["records"][0]["litcurate"]["candidate_context"]["evidence_locations"][0][
            "locator"
        ]
        == "Table 1"
    )


def test_standard_library_xlsx_reader_handles_sparse_inline_strings(tmp_path):
    workbook = """<?xml version="1.0" encoding="UTF-8"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
 <sheets><sheet name="database" sheetId="1" r:id="rId1"/></sheets>
</workbook>"""
    relationships = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
 <Relationship Id="rId1"
  Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"
  Target="worksheets/sheet1.xml"/>
</Relationships>"""
    worksheet = """<?xml version="1.0" encoding="UTF-8"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
 <sheetData>
  <row r="1">
   <c r="A1" t="inlineStr"><is><t>doi</t></is></c>
   <c r="C1" t="inlineStr"><is><t>year</t></is></c>
  </row>
  <row r="2">
   <c r="A2" t="inlineStr"><is><t>10.1000/example</t></is></c>
   <c r="C2"><v>2026</v></c>
  </row>
 </sheetData>
</worksheet>"""
    path = tmp_path / "fixture.xlsx"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", relationships)
        archive.writestr("xl/worksheets/sheet1.xml", worksheet)

    assert litcurate_candidates.load_xlsx(path) == [
        {"doi": "10.1000/example", "year": 2026}
    ]
