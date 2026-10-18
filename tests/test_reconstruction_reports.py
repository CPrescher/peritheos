"""Report links resolve exactly, and figures never manufacture scientific data."""

import copy
from pathlib import Path
from urllib.parse import urlsplit

import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    convert_eosmat,
    eosmat_schema,
    load_eosmat,
    validate_eosmat_document,
)
from peritheos.errors import EosmatError
from scripts.generate_reconstruction_reports import (
    comparison_svg,
    current_figure,
    explicit_reports,
    generated_outputs,
    illustrated_mentions,
    report_page,
)
from scripts.reconstruction_links import DOCS_URL, report_links
from scripts.standardize_reconstruction_metadata import standardize_document
from tests.test_eosmat_reconstruction import sample

ROOT = Path(__file__).resolve().parents[1]


def linked_sample():
    document, row = sample()
    row["material_identifier"] = document["identifier"]
    return standardize_document(document, {row["record_identifier"]: row}), row


def test_navigation_upgrade_preserves_every_existing_field_and_refuses_different_findings():
    document, row = linked_sample()
    before = copy.deepcopy(document)
    assessment = document["eos_records"][0]["evidence"]["reproduction"]["assessment"]
    assessment.pop("report_url")
    assessment.pop("report_path")
    upgraded = standardize_document(document, {row["record_identifier"]: row})
    assert upgraded == before
    assert "report_url" not in assessment
    assert standardize_document(upgraded, {row["record_identifier"]: row}) == upgraded
    assessment["report_url"] = "https://example.org/authored-report"
    with pytest.raises(ValueError, match="Existing assessment differs"):
        standardize_document(document, {row["record_identifier"]: row})


def test_linked_report_keeps_source_values_missing_errors_and_scope_clear(tmp_path):
    document, _ = linked_sample()
    original = copy.deepcopy(document)
    page = report_page(document, root=tmp_path)
    assert "Parity not achieved for this comparison" in page
    assert "160" in page and "180.12345678901235" in page
    assert "Unknown" not in page or "not been inferred" in page
    assert "Missing errors mean unknown, not zero" in page
    assert "not an observation or residual plot" in page
    assert "Fixed/adopted values are constraints" in page
    assert "Synthetic conditional diagnostic" in page
    assert document == original
    legacy = convert_eosmat(document, version=3)
    assert (
        legacy["eos_records"][0]["reproduction"]["assessment"]
        == document["eos_records"][0]["evidence"]["reproduction"]["assessment"]
    )
    assessment = document["eos_records"][0]["evidence"]["reproduction"]["assessment"]
    assessment.pop("parameters")
    page = report_page(document, root=tmp_path)
    assert "Why no coefficient graph is shown" in page
    assert "No numeric source/reconstruction pair" in page
    assert ".svg)" not in page
    assert comparison_svg(assessment) is None


def test_figure_is_reproducible_fingerprinted_and_does_not_change_input():
    document, _ = linked_sample()
    assessment = document["eos_records"][0]["evidence"]["reproduction"]["assessment"]
    original = copy.deepcopy(assessment)
    svg = comparison_svg(assessment)
    assert svg and svg == comparison_svg(assessment)
    assert current_figure(svg, assessment)
    assert not current_figure(svg.replace(b"Source", b"Edited", 1), assessment)
    assert not current_figure(b"<svg/>", assessment)
    assessment["parameters"][0]["refit"] += 1
    assert not current_figure(svg, assessment)
    assessment["parameters"][0]["refit"] -= 1
    assert assessment == original


@pytest.mark.parametrize("field", ["report_url", "report_path"])
@pytest.mark.parametrize("value", [None, "", [], 123])
def test_invalid_navigation_fields_are_rejected_by_schema_and_python(field, value):
    document, _ = linked_sample()
    document["eos_records"][0]["evidence"]["reproduction"]["assessment"][field] = value
    assert not Draft202012Validator(eosmat_schema(4)).is_valid(document)
    with pytest.raises(EosmatError):
        validate_eosmat_document(document)


def test_supporting_reports_use_only_explicit_existing_local_references(tmp_path):
    directory = tmp_path / "docs/literature-reproductions"
    directory.mkdir(parents=True)
    (directory / "existing.md").write_text("# Explicit report\n", encoding="utf-8")
    record = {
        "notes": "See docs/literature-reproductions/existing.md and literature-reproductions/missing.md"
    }
    assert explicit_reports(record, tmp_path) == [
        "literature-reproductions/existing.md"
    ]


def test_illustrated_reports_match_exact_identifiers_without_guessing_from_a_material_name(
    tmp_path,
):
    directory = tmp_path / "docs/literature-reproductions"
    directory.mkdir(parents=True)
    (directory / "exact.md").write_text(
        "# Existing plot\n`example_bm3`\n![Recorded plot](plot.png)", encoding="utf-8"
    )
    (directory / "different.md").write_text(
        "# Other EOS\n`example_bm3_thermal`\n![Other plot](other.png)", encoding="utf-8"
    )
    assert illustrated_mentions("example_bm3", tmp_path) == [
        "literature-reproductions/exact.md"
    ]
    assert illustrated_mentions("example", tmp_path) == []


def test_every_bundled_link_opens_its_record_and_each_graph_matches_its_input():
    records = charts = 0
    outputs = generated_outputs(plots=False)
    for relative, payload in outputs.items():
        assert (ROOT / "docs/reconstruction-reports" / relative).read_bytes() == payload
    for path in sorted((ROOT / "peritheos/data/materials").glob("*.eosmat")):
        document = load_eosmat(path)
        for record in document["eos_records"]:
            assessment = record["evidence"]["reproduction"]["assessment"]
            expected = report_links(document["identifier"], record["identifier"])
            assert all(assessment[k] == v for k, v in expected.items())
            assert assessment["report_url"].startswith(DOCS_URL)
            local, fragment = assessment["report_path"].split("#", 1)
            public = urlsplit(assessment["report_url"])
            assert public.fragment == fragment
            assert public.path.endswith("/" + document["identifier"] + "/")
            assert f'<a id="{fragment}"></a>' in (ROOT / local).read_text(
                encoding="utf-8"
            )
            figure = (
                ROOT
                / "docs/reconstruction-reports/figures"
                / (record["identifier"] + ".svg")
            )
            page = (ROOT / local).read_text(encoding="utf-8")
            if figure.exists():
                assert current_figure(figure.read_bytes(), assessment)
                assert f"figures/{record['identifier']}.svg" in page
                charts += 1
            records += 1
    assert records == 663 and charts > 0
