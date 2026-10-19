"""Report links resolve exactly, and figures never manufacture scientific data."""

import copy
from pathlib import Path
from urllib.parse import urlsplit
from xml.etree import ElementTree

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    convert_eosmat,
    eosmat_schema,
    get_material_document,
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
from scripts.reconstruction_figures import current_curve, curve_svg
from scripts.reconstruction_links import DOCS_URL, report_links
from scripts.reconstruction_plot_data import curve_comparison
from scripts.report_numbers import value_error
from scripts.report_plot_themes import dark_path, dark_svg, themed_images
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
    assert "160" in page and "180.1" in page and "2.5" in page
    assert "180.12345678901235" not in page
    assert "Unknown" not in page or "not been inferred" in page
    assert "Missing errors mean unknown, not zero" in page
    assert "not an observation or residual plot" in page
    assert "Fixed/adopted values are constraints" in page
    assert "Synthetic conditional diagnostic" in page
    assert "Curve/residual comparison unavailable" in page
    assert "Recorded pressure RMS error" in page
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


def test_dark_variant_changes_only_paint_and_preserves_plot_geometry_and_labels():
    document, _ = linked_sample()
    assessment = document["eos_records"][0]["evidence"]["reproduction"]["assessment"]
    light = comparison_svg(assessment)
    dark = dark_svg(light)
    assert dark == dark_svg(light)
    assert b"fill: #ffffff" not in dark
    assert b"fill: #1e2129" in dark
    assert b"text { fill: #e6e8ef; }" in dark
    assert b"#82baff" in dark and b"#ffad73" in dark
    assert b"peritheos-dark source=" in dark
    assert not current_figure(dark, assessment)

    def scientific_content(payload):
        svg = ElementTree.fromstring(payload)
        return [
            (
                element.tag,
                {
                    k: v
                    for k, v in element.attrib.items()
                    if k not in {"style", "fill", "stroke"}
                },
                element.text,
            )
            for element in svg.iter()
            if element.tag.rsplit("}", 1)[-1] != "style"
        ]

    assert scientific_content(light) == scientific_content(dark)
    assert themed_images("Comparison", "figures/example.svg") == [
        "![Comparison](figures/example.svg#only-light)",
        "![Comparison](figures/example-dark.svg#only-dark)",
    ]


@pytest.mark.parametrize(
    "value,error,expected",
    [
        (25.447259, 0.091348, ("25.45", "0.09")),
        (24.9234, 0.1125, ("24.92", "0.11")),
        (180.123456, 2.5, ("180.1", "2.5")),
        (160.12345, 3.18, ("160", "3")),
        (-4.237, 0.29, ("-4.24", "0.29")),
        (12.3456, 0.099, ("12.35", "0.10")),
        (0.00000123456, 0.000000121, ("12.3e-7", "1.2e-7")),
        (123.456789, None, ("123.457", "Not reported")),
        (123.456789, 0, ("123.457", "0 (reported)")),
    ],
)
def test_display_rounding_matches_uncertainty_place(value, error, expected):
    assert value_error(value, error) == expected


def test_curve_residuals_replay_recorded_shen_inputs_without_optimizing(monkeypatch):
    def no_fit(*args, **kwargs):
        raise AssertionError("Report generation must not fit")

    monkeypatch.setattr("scripts.validate_primary_eos_refits.fit_rt_eos", no_fit)
    document = get_material_document("nacl_b2")
    record = next(
        r
        for r in document["eos_records"]
        if r["identifier"] == "nacl_b2_shen_2026_vinet_2"
    )
    original = copy.deepcopy(document)
    comparison, reason = curve_comparison(document, record)
    assert comparison is not None, reason
    assert len(comparison.series.pressure) == 125
    assert comparison.source_rmse == pytest.approx(0.6223809764, abs=1e-9)
    assert comparison.refit_rmse == pytest.approx(0.2932287344, abs=1e-9)
    assert np.all(comparison.source_pressure - comparison.series.pressure < 0)
    assert document == original
    svg = curve_svg(comparison)
    assert svg == curve_svg(comparison)
    assert current_curve(svg, comparison)
    assert b"Observed / calibrated pressure (GPa)" in svg
    assert b"0.622 GPa" in svg and b"0.293 GPa" in svg
    assert not current_curve(svg.replace(b"Original", b"Edited", 1), comparison)
    comparison.series.pressure[0] += 0.1
    assert not current_curve(svg, comparison)


def test_fei_staged_rows_do_not_duplicate_cold_observations():
    document = get_material_document("nacl_b2")
    record = next(
        r for r in document["eos_records"] if r["identifier"] == "nacl_b2_fei_2007"
    )
    comparison, reason = curve_comparison(document, record)
    assert comparison is not None, reason
    assert len(comparison.series.pressure) == 24
    assert list(np.unique(comparison.series.temperature, return_counts=True)[1]) == [
        12,
        12,
    ]
    assert comparison.series.pressure_sigma is None
    assert comparison.series.volume_sigma is None
    record["reproduction"]["assessment"]["fit"]["observations"] = 36
    comparison, reason = curve_comparison(document, record)
    assert comparison is None and "row count differs" in reason


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
    records = charts = curves = 0
    outputs = generated_outputs(plots=False)
    index = outputs[Path("index.md")].decode()
    assert "471 of 471" in index
    assert "3 explicitly partial replays" in index
    assert "None among assessments with stored coefficient pairs" in index
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
                comparison, reason = curve_comparison(document, record)
                assert comparison is not None, (record["identifier"], reason)
                curve = figure.with_stem(figure.stem + "-curves")
                assert current_curve(curve.read_bytes(), comparison)
                curves += 1
            records += 1
    assert records == 663 and charts == curves == 471
    for light in (ROOT / "docs/reconstruction-reports/figures").glob("*.svg"):
        if not light.stem.endswith("-dark"):
            dark = light.with_name(Path(dark_path(light.name)).name)
            assert dark.read_bytes() == dark_svg(light.read_bytes())
