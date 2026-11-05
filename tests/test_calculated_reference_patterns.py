"""New powder arrays must be reproducible calculations with scoped source terms."""

import copy
import hashlib
import json
import math

import pytest

from peritheos import get_material_document, load_eosmat
from peritheos_build import ROOT, digest, stage_core
from scripts.merge_source_bindings import accepts_merge_binding
from scripts.scientific_source_snapshot import (
    assert_snapshot_lineage,
    scientific_sha256,
    scientific_snapshot,
)

PROOF = json.loads(
    (ROOT / "docs/data/calculated-pattern-replacement.json").read_bytes()
)
INPUTS = json.loads(
    (ROOT / "peritheos/data/calculated-pattern-inputs.json").read_bytes()
)


def test_replacement_preserves_every_non_authorized_field_and_scientific_history():
    assert scientific_sha256() == PROOF["after_catalog_scientific_sha256"]
    assert_snapshot_lineage(
        scientific_snapshot(), PROOF["before_catalog_scientific_sha256"]
    )
    assert (
        hashlib.sha256(
            (ROOT / "peritheos/data/calculated-pattern-inputs.json").read_bytes()
        ).hexdigest()
        == PROOF["inputs_sha256"]
    )
    assert len(PROOF["material_bindings"]) == 10
    assert PROOF["original_peak_count"] == 452
    assert PROOF["calculated_peak_count"] == 1154
    for material, binding in PROOF["material_bindings"].items():
        payload = (ROOT / binding["path"]).read_bytes()
        document = json.loads(payload)
        assert hashlib.sha256(payload).hexdigest() == binding["after_sha256"]
        assert (
            digest(
                {
                    k: v
                    for k, v in document.items()
                    if k not in binding["changed_root_fields"]
                }
            )
            == binding["unchanged_projection_sha256"]
        )
        assert digest(document["peaks"]) == binding["calculated_peaks_sha256"]
        assert digest(document["peaks"]) != binding["original_peaks_sha256"]
        assert document["source"]["pattern_kind"] == "independently_calculated_xray"
        assert "not measured" in document["notes"]
        assert len(document["peaks"]) == binding["calculated_peak_count"]


@pytest.mark.parametrize(
    "mutation", ["eos", "unrelated_peaks", "new_pattern", "proxy_sites"]
)
def test_historical_projection_cannot_accept_other_scientific_changes(mutation):
    snapshot = scientific_snapshot()
    documents = {d["identifier"]: d for d in snapshot["documents"]}
    if mutation == "eos":
        documents["feh2"]["eos_records"][0]["parameters"]["K0"]["value"] += 1
    elif mutation == "unrelated_peaks":
        documents["gold"]["peaks"] = [[1, 0, 0, 1, 100]]
    elif mutation == "new_pattern":
        documents["feh2"]["peaks"][0][-1] += 1
    else:
        documents["fe3s"]["atom_sites"][0]["x"] += 0.01
    with pytest.raises(AssertionError):
        assert_snapshot_lineage(snapshot, PROOF["before_catalog_scientific_sha256"])


@pytest.mark.parametrize("material", list(INPUTS["materials"]))
def test_full_recalculation_contents_spacings_absences_and_normalization(material):
    gemmi = pytest.importorskip(
        "gemmi", reason="maintenance reproduction: --with gemmi==0.7.4"
    )
    from scripts.calculate_reference_patterns import calculate

    entry = INPUTS["materials"][material]
    peaks, contents = calculate(entry)
    assert contents == pytest.approx(
        PROOF["material_bindings"][material]["cell_contents"]
    )
    assert peaks == get_material_document(material)["peaks"]
    assert max(p[4] for p in peaks) == 100
    assert all(0 < p[4] <= 100 and 0.8 <= p[3] <= 20 for p in peaks)
    assert len({p[3] for p in peaks}) == len(peaks)
    assert [p[3] for p in peaks] == sorted((p[3] for p in peaks), reverse=True)
    a, b, c = (entry["lattice"][k] for k in ("a", "b", "c"))
    operations = gemmi.SpaceGroup(entry["space_group"]).operations()
    for h, k, ell, d, _ in peaks:
        assert not operations.is_systematically_absent([h, k, ell])
        inverse_d2 = (
            4 / 3 * (h * h + h * k + k * k) / a**2 + ell**2 / c**2
            if entry["lattice"]["gamma"] == 120
            else h * h / a**2 + k * k / b**2 + ell**2 / c**2
        )
        assert d == pytest.approx(1 / math.sqrt(inverse_d2), abs=5.1e-11)


def test_rocksalt_intensity_ratio_matches_independent_two_sublattice_formula():
    gemmi = pytest.importorskip("gemmi")
    doc = get_material_document("mgfe60o")
    a = doc["lattice"]["a"]
    peaks = {tuple(row[:3]): row[4] for row in doc["peaks"]}

    def line(norm2, sign, multiplicity):
        stol2 = norm2 / (4 * a * a)
        f = {e: gemmi.Element(e).it92.calculate_sf(stol2) for e in ("Mg", "Fe", "O")}
        amplitude = 4 * (0.4 * f["Mg"] + 0.6 * f["Fe"] + sign * f["O"])
        return multiplicity * amplitude**2

    # Our deterministic representative for the cubic 200 family is 002.
    expected = line(3, -1, 8) / line(4, 1, 6)
    assert peaks[(1, 1, 1)] / peaks[(0, 0, 2)] == pytest.approx(expected, rel=2e-8)


def test_input_grants_and_proxy_limits_are_not_a_blanket_article_license():
    models = INPUTS["materials"]
    assert sum("cod_id" in e for e in models.values()) == 6
    for material in ("feh2", "sno2_cubic_27gpa", "sno2_pa_3_at_48gpa"):
        assert "No separate source grant" in models[material]["source_license"]
        assert (
            "source publications are separate"
            in models[material]["input_format_license"]
        )
    assert models["sno2_cubic_27gpa"]["source_doi"] == "10.1107/S1600577519003904"
    assert "not a refined Fe3S" in models["fe3s"]["qualification"]
    assert "3.7 wt% Fe" in models["phase_d"]["qualification"]
    assert "0.88 Mg + 0.12 Fe" in models["perovskite_orthorhombic"]["qualification"]


def test_core_ships_reproduction_inputs_and_scoped_notice_without_observations(
    tmp_path,
):
    core = stage_core(tmp_path / "core") / "peritheos/data"
    assert (core / "calculated-pattern-inputs.json").read_bytes() == (
        ROOT / "peritheos/data/calculated-pattern-inputs.json"
    ).read_bytes()
    assert not list(core.rglob("*.csv"))
    for material, entry in INPUTS["materials"].items():
        if "resource" in entry:
            assert (
                hashlib.sha256((core / entry["resource"]).read_bytes()).hexdigest()
                == entry["sha256"]
            )
        document = json.loads((core / f"materials/{material}.eosmat").read_bytes())
        assert (
            document["source"]["pattern_input_source_terms"] == entry["source_license"]
        )
    assert (core / "datasets/calculated-patterns.LICENSE.md").is_file()


def test_frozen_source_exceptions_accept_only_exact_reviewed_pairs():
    for context, entries in PROOF["bindings"].items():
        for name, binding in entries.items():
            path = (
                name
                if name.endswith(".eosmat")
                else f"peritheos/data/materials/{name}.eosmat"
            )
            document = copy.deepcopy(load_eosmat(ROOT / path))
            for record in document["eos_records"]:
                record.get("evidence", {}).get("reproduction", {}).pop(
                    "assessment", None
                )
            actual = hashlib.sha256(
                json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
            assert actual == binding["replacement_sha256"]
            old = binding["historical_sha256"]
            assert accepts_merge_binding(context, name, old, actual)
            assert not accepts_merge_binding(context, name, "0" * 64, actual)
            assert not accepts_merge_binding(context, name, old, "0" * 64)
            assert not accepts_merge_binding("unreviewed", name, old, actual)
