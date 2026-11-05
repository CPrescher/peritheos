"""COD reacquisition must preserve values, proxy limits and old source evidence."""

import copy
import hashlib
import json
import math
import re
from collections import Counter
from fractions import Fraction

from peritheos import Material, get_material_document, load_eosmat
from peritheos.eosmat import validate_eosmat_document
from peritheos_build import ROOT, digest, stage_core
from scripts.merge_source_bindings import accepts_merge_binding
from scripts.scientific_source_snapshot import (
    assert_snapshot_lineage,
    scientific_snapshot,
)

PROOF_PATH = ROOT / "docs/data/post-perovskite-cod-replacement.json"
PROOF = json.loads(PROOF_PATH.read_bytes())
PURE = "mgsio3_post_perovskite"
PROXY = "mg090fe009al0005ca0005sio3_post_perovskite"


def _cif_rows(text, last_header):
    """Read the specified, checksum-bound simple loop; not a general CIF parser."""
    lines = text.splitlines()
    start = lines.index(last_header) + 1
    end = next(
        (i for i in range(start, len(lines)) if lines[i].startswith(("loop_", "_"))),
        len(lines),
    )
    return [line.split() for line in lines[start:end] if line.strip()]


def _source():
    payload = (ROOT / "peritheos/data" / PROOF["resource"]["path"]).read_bytes()
    assert hashlib.sha256(payload).hexdigest() == PROOF["resource"]["sha256"]
    text = payload.decode()
    atoms = {
        row[0]: tuple(Fraction(x) for x in row[1:])
        for row in _cif_rows(text, "_atom_site_fract_z")
    }
    lengths = [
        float(re.search(rf"(?m)^_cell_length_{axis}\s+(\S+)", text)[1])
        for axis in ("a", "b", "c")
    ]
    return text, atoms, lengths


def test_cod_cell_symmetry_and_full_site_orbits_establish_mgsio3_z4():
    text, atoms, lengths = _source()
    assert lengths == [2.456, 8.042, 6.093]
    assert "_space_group_IT_number           63" in text
    assert "_cell_formula_units_Z            4" in text
    assert "_diffrn_ambient_pressure         1.21e+08" in text  # kPa = 121 GPa
    assert "_diffrn_ambient_temperature      300" in text
    assert "ab initio calculations using observed cell" in text
    operations = [
        r[0].split(",") for r in _cif_rows(text, "_space_group_symop_operation_xyz")
    ]
    assert len(operations) == 16

    def apply(expression, xyz):
        result = Fraction(0)
        for term in re.findall(r"[+-]?[^+-]+", expression):
            sign = -1 if term.startswith("-") else 1
            token = term.lstrip("+-")
            result += sign * (xyz[token] if token in xyz else Fraction(token))
        return result % 1

    contents = Counter()
    for label, position in atoms.items():
        xyz = dict(zip(("x", "y", "z"), position))
        orbit = {tuple(apply(e, xyz) for e in op) for op in operations}
        contents[re.sub(r"\d+$", "", label)] += len(orbit)
    assert contents == {"Mg": 4, "Si": 4, "O": 12}


def test_both_cells_and_sites_regenerate_from_cod_without_numeric_change():
    _, atoms, lengths = _source()
    factor = (162.2 / math.prod(lengths)) ** (1 / 3)
    pure_axes = [round(x * factor, 10) for x in lengths]
    proxy_factor = (164.9 / math.prod(pure_axes)) ** (1 / 3)
    proxy_axes = [round(x * proxy_factor, 10) for x in pure_axes]
    for material, axes in ((PURE, pure_axes), (PROXY, proxy_axes)):
        document = get_material_document(material)
        validate_eosmat_document(document)
        assert [document["lattice"][k] for k in ("a", "b", "c")] == axes
        assert document["formula_units_per_cell"] == 4
        assert document["space_group_number"] == 63
        assert not document["peaks"]
        for site, label in zip(document["atom_sites"], ("Si", "Mg", "O1", "O2")):
            assert tuple(site[k] for k in ("x", "y", "z")) == tuple(
                float(x) for x in atoms[label]
            )
            assert site["occupancy"] == 1
        assert document["source"]["structure_cod_id"] == "9009217"
        assert document["source"]["structure_reference_pressure_gpa"] == 121
        assert (
            "not experimentally refined"
            in document["source"]["structure_determination_method"]
        )
        assert "ICSD 158959" not in json.dumps(document.get("source", {}))
        assert Material.from_eosmat(document).formula == document["formula"]
    proxy = get_material_document(PROXY)
    assert proxy["source"]["topology_proxy"] is True
    assert proxy["source"]["structure_topology_proxy"] is True
    assert "no Fe/Al/Ca site occupancies" in proxy["source"]["topology_proxy_scope"]


def test_only_root_source_and_notes_change_with_exact_frozen_source_pairs():
    assert_snapshot_lineage(scientific_snapshot(), PROOF["catalog_scientific_sha256"])
    for material, binding in PROOF["material_bindings"].items():
        path = ROOT / "peritheos/data/materials" / f"{material}.eosmat"
        payload = path.read_bytes()
        assert hashlib.sha256(payload).hexdigest() == binding["after_sha256"]
        raw = json.loads(payload)
        projection = {k: v for k, v in raw.items() if k not in ("notes", "source")}
        assert digest(projection) == binding["unchanged_except_notes_and_source_sha256"]
        assert set(PROOF["changed_root_fields"][material]) == {"notes", "source"}
        normalized = copy.deepcopy(load_eosmat(path))
        for record in normalized["eos_records"]:
            record.get("evidence", {}).get("reproduction", {}).pop("assessment", None)
        actual = hashlib.sha256(
            json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        name = path.relative_to(ROOT).as_posix()
        before = PROOF["bindings"]["release_document"][name]["historical_sha256"]
        assert accepts_merge_binding("release_document", name, before, actual)
        assert not accepts_merge_binding("release_document", name, "0" * 64, actual)
        assert not accepts_merge_binding("release_document", name, before, "0" * 64)
        assert not accepts_merge_binding("unreviewed", name, before, actual)


def test_runtime_exports_keep_cod_resource_and_explicit_proxy_qualifications(tmp_path):
    core = stage_core(tmp_path / "core") / "peritheos/data"
    relative = PROOF["resource"]["path"]
    assert (core / relative).read_bytes() == (
        ROOT / "peritheos/data" / relative
    ).read_bytes()
    for material in (PURE, PROXY):
        document = json.loads((core / "materials" / f"{material}.eosmat").read_bytes())
        source = document["source"]
        assert source["structure_dataset"] == relative
        assert source["structure_dataset_sha256"] == PROOF["resource"]["sha256"]
        assert "CC0" in source["structure_license"]
        assert source["structure_topology_proxy"] is (material == PROXY)
        if material == PROXY:
            assert "no Fe/Al/Ca" in source["structure_topology_proxy_scope"]


def test_current_review_has_five_cod_components_and_no_icsd_case():
    review = json.loads(
        (ROOT / "peritheos/data/eos-structure-rights-review.json").read_bytes()
    )
    assert review["counts"]["assessments"]["cod_source_and_numerical_changes"] == 5
    assert "icsd_acquisition_review_required" not in review["counts"]["assessments"]
    assert sum(e["assessment"].endswith("required") for e in review["entries"]) == 7
    assert (
        review["post_perovskite_cod_replacement"]["sha256"]
        == hashlib.sha256(PROOF_PATH.read_bytes()).hexdigest()
    )
