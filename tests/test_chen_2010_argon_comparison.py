"""Preserve source/refit identities and prevent accidental scientific promotion."""

import json

import numpy as np
import pytest

from peritheos import get_material, get_material_document
from peritheos.catalog import get_eos_record
from peritheos.materials import Material
from scripts.compare_chen_2010_argon import ORIGINAL, REFIT, REPORT, compare
from scripts.generate_paper_investigation_ledger import build_papers
from scripts.reproduce_chen_2010_argon import MASS_PER_CELL, ROOT


def test_original_record_is_explicit_diagnostic_and_preserves_local_constraints():
    doc = get_material_document("argon_fcc")
    raw = next(r for r in doc["eos_records"] if r["identifier"] == ORIGINAL)
    assert raw["record_kind"] == "diagnostic"
    assert raw["catalog_access"] == "explicit_selection"
    assert raw["scientific_validation"]["status"] == "not_reproduced"
    assert raw["reproduction"]["do_not_use"] is True
    assert raw["fit_datasets"] == []
    assert not raw["default"]
    assert "DO NOT USE" in raw["label"]
    assert get_material("argon_fcc").default_record().identifier != ORIGINAL
    with pytest.raises(ValueError, match="not_reproduced"):
        Material.from_eosmat(doc, record_identifiers=[ORIGINAL])
    record = get_eos_record(ORIGINAL)
    rho = 2.18
    step = 1e-4

    def p(density):
        return record.pressure(MASS_PER_CELL / density, 290)

    first = (p(rho + step) - p(rho - step)) / (2 * step)
    second = (p(rho + step) - 2 * p(rho) + p(rho - step)) / step**2
    assert p(rho) == pytest.approx(2, abs=1e-10)
    assert rho * first == pytest.approx(15.1, rel=1e-6)
    assert 1 + rho * second / first == pytest.approx(5.4, rel=1e-5)


def test_both_records_link_explanation_and_preserve_distinct_status():
    records = get_material_document("argon_fcc")["eos_records"]
    for identifier in [ORIGINAL, REFIT]:
        record = next(r for r in records if r["identifier"] == identifier)
        details = record["reproduction"]
        assert (ROOT / details["documentation"]).is_file()
        assert (ROOT / details["comparison_figure"]).is_file()
        assert set(details["related_record_identifiers"]) == {ORIGINAL, REFIT}
    refit = next(r for r in records if r["identifier"] == REFIT)
    assert refit["record_kind"] == "refit"
    assert (
        refit["scientific_validation"]["original_publication_reproduction_status"]
        == "not_reproduced"
    )
    assert refit["parameter_covariance"] is None


def test_comparison_uses_fixed_pressure_and_reproduces_saved_report():
    actual = compare()
    expected = json.loads(REPORT.read_text())
    assert actual["pressure_range_gpa"] == expected["pressure_range_gpa"]
    for name, values in actual["plot"]["density_g_cm3"].items():
        np.testing.assert_allclose(
            values, expected["plot"]["density_g_cm3"][name], atol=1e-10
        )
    point = next(r for r in actual["checkpoints"] if r["pressure_gpa"] == 20)
    density = point["density_g_cm3"]
    for name, difference in point["volume_difference_from_refit_percent"].items():
        assert difference == pytest.approx(100 * (density["refit"] / density[name] - 1))
    assert point["volume_difference_from_refit_percent"][
        "reported_constants"
    ] == pytest.approx(11.14877, abs=1e-4)
    assert point["volume_difference_from_refit_percent"][
        "published_curve"
    ] == pytest.approx(-0.28447, abs=1e-4)


def test_paper_outcome_never_promotes_independent_refit_to_original_reproduction():
    papers, _ = build_papers()
    chen = next(
        p
        for p in papers
        if str(p.get("doi", "")).lower() == "10.1103/physrevb.81.144110"
    )
    assert chen["outcome"] == "unreproduced_original_with_independent_refit"
    assert {r["record_identifier"] for r in chen["records"]} == {ORIGINAL, REFIT}
    ledger = (ROOT / "docs/paper-investigation-ledger.md").read_text()
    assert "data/chen-2010-argon-comparison.png" in ledger
    assert "DO NOT USE" in ledger
    assert "no error in Chen's" in ledger
