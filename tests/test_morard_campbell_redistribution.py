"""Trace original observations and keep new calculations independently named."""

import csv
import hashlib
import io
import json
from decimal import Decimal

import pytest

from peritheos import Dataset, get_material_document
from peritheos_build import ROOT, digest
from scripts.generate_campbell_2009_polynomial_tables import (
    DATA,
    ORIGINAL,
    OUTPUT,
    SOURCE,
    compare_original,
    evaluate,
    generate,
)


def test_morard_all_cold_rows_match_identified_sources_without_rewriting_errors():
    trace = json.loads((DATA / "morard-cold-source-lineage.json").read_bytes())
    for name, checksum in trace["bindings"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == checksum
    original = list(
        csv.DictReader(
            io.StringIO(
                (DATA / "fes-morard-2026-author-literature-cold.csv").read_text()
            )
        )
    )
    assert len(original) == len(trace["rows"]) == 21
    previous = json.loads(
        (ROOT / trace["sources"]["ohfuji_2007"]["prior_receipt"]).read_bytes()
    )
    prior = {
        r["author_literature_id"]: r
        for r in previous["cold_vi_matches_to_author_literature"]
    }
    differences = []
    for row, source in zip(original, trace["rows"]):
        assert row["observation_id"] == source["observation_id"]
        published_volume = Decimal(
            source["source_cell_volume_token_angstrom3"].split("(")[0]
        )
        assert Decimal(row["molar_volume_cm3_mol"]) == published_volume * Decimal(
            "0.15055"
        )
        differences.append(
            Decimal(row["pressure_gpa"])
            - Decimal(source["source_pressure_token_gpa"].split("(")[0])
        )
        assert float(row["pressure_uncertainty_gpa"]) == pytest.approx(
            0.01 * float(row["pressure_gpa"])
        )
        assert float(row["molar_volume_uncertainty_cm3_mol"]) == pytest.approx(
            0.005 * float(row["molar_volume_cm3_mol"])
        )
        if source["source"] == "ohfuji_2007":
            assert (
                source["source_cell_volume_token_angstrom3"]
                == prior[row["observation_id"]]["volume_cell_token"]
            )
    assert differences.count(Decimal(0)) == 20
    assert differences[11] == Decimal("0.1")
    # The newly identified source includes a non-monotone measured volume pair;
    # it is retained, not smoothed to make an EOS fit look better.
    assert trace["rows"][15]["source_cell_volume_token_angstrom3"] == "78.59(9)"
    assert trace["rows"][16]["source_cell_volume_token_angstrom3"] == "78.63(9)"


def test_campbell_generator_uses_coefficients_without_original_output_access(
    tmp_path, monkeypatch
):
    import scripts.generate_campbell_2009_polynomial_tables as generator

    monkeypatch.setattr(generator, "ORIGINAL", tmp_path / "absent-original.csv")
    assert generate() == OUTPUT.read_bytes()
    configuration = json.loads(SOURCE.read_bytes())
    iw = configuration["coefficients"]["Table S4"]
    nno = configuration["coefficients"]["Table S5"]
    # Separate hand arithmetic at a nonzero P, including the P**4 NNO term.
    assert evaluate(10, 1500, iw) == pytest.approx(
        6.54106
        + 0.0012324 * 10
        + (-28163.6 + 546.32 * 10 - 1.13412 * 100 + 0.0019274 * 1000) / 1500
    )
    assert evaluate(10, 1500, nno) == pytest.approx(
        8.699
        + 0.01642 * 10
        - 0.0002755 * 100
        + 2.683e-6 * 1000
        - 1.015e-8 * 10000
        + (-24205 + 444.73 * 10 - 0.59288 * 100 + 0.0015292 * 1000) / 1500
    )
    for p, t in [(-1, 1500), (101, 1500), (10, 999), (10, 2601), (float("nan"), 1500)]:
        with pytest.raises(ValueError):
            evaluate(p, t, iw)


def test_campbell_approximations_match_all_original_keys_but_are_not_originals():
    assert (
        hashlib.sha256(ORIGINAL.read_bytes()).hexdigest()
        == "95f0bfe5c2541885fc80bee1b3759b4ab2f6f23a9e9e13251b9dc0f81f1a2699"
    )
    comparison = compare_original(generate())
    assert {v["rows"] for v in comparison.values()} == {1717}
    iw, nno = comparison["Table S4"], comparison["Table S5"]
    assert iw["rmse_log10_fugacity"] == pytest.approx(0.01757083636)
    assert iw["maximum_absolute_difference_log10_fugacity"] == pytest.approx(
        0.05136112623
    )
    assert iw["maximum_absolute_difference_log10_fugacity"] > 0.046
    assert nno["maximum_absolute_difference_log10_fugacity"] < 0.054
    assert generate() != ORIGINAL.read_bytes()


@pytest.mark.parametrize("material", ["fe_fcc", "feo", "nickel", "nickel_oxide_b1"])
def test_new_calculated_declaration_preserves_the_complete_previous_material(material):
    proof = json.loads(
        (ROOT / "docs/data/morard-campbell-redistribution-followup.json").read_bytes()
    )
    document = json.loads(
        (ROOT / "peritheos/data/materials" / f"{material}.eosmat").read_bytes()
    )
    added = document["datasets"].pop()
    assert added["identifier"] == proof["prior_documents"][material]["added_identifier"]
    assert digest(document) == proof["prior_documents"][material]["semantic_sha256"]
    assert added["kind"] == "model_output"
    assert added["used_by_eos_records"] == []
    assert added["provenance"]["source_parameters"] == json.loads(SOURCE.read_bytes())
    assert len(Dataset.from_mapping(added)) == 3434
    loaded = get_material_document(material)
    assert any(d["identifier"] == added["identifier"] for d in loaded["datasets"])
