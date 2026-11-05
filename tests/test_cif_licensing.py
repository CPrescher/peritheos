"""CIF source grants must be file-scoped and leave structural evidence intact."""

import csv
import hashlib
import json
from pathlib import Path

import pytest

from scripts.inventory_distribution_licenses import match_scope

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"


def review():
    return json.loads((DATA / "cif-license-review.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    "filename,expected",
    [
        (
            "vanadium-james-1960-9012770.cif",
            "3f81e29961c91b096d22f3dc3313b08dd87e56bd1449e0b7a24eb2dc01e41548",
        ),
        (
            "superhydrous-phase-b-koch-muller-2005-mkm0105.cif",
            "9118860e2e0e65a8a8fac0743d2ad154d5611d80f40b53e2444715dfad12adf3",
        ),
        (
            "brucite-zigan-1967-cod-1000054.cif",
            "1defd6b791e80c0ce967ea13a4457559e46d23028fb6f9796c5215871542eba1",
        ),
    ],
)
def test_review_preserves_exact_original_cif_bytes(filename, expected):
    row = next(r for r in review()["files"] if r["file"] == filename)
    assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == expected
    assert row["local_sha256"] == expected
    if "source_versions" in row:
        assert len(row["source_versions"]) == 2
        exact = [r for r in row["source_versions"] if r["exact_byte_match"]]
        assert len(exact) == 1
        assert exact[0]["sha256"] == expected
        payload = b"".join(
            line
            for line in (DATA / filename).read_bytes().splitlines(keepends=True)
            if not line.startswith((b"#$Date:", b"#$URL:"))
        )
        assert (
            hashlib.sha256(payload).hexdigest()
            == row["cif_payload_sha256_ignoring_only_svn_date_url"]
        )


def test_current_cod_review_is_scoped_and_historical_inventory_is_preserved():
    scopes = json.loads(
        (ROOT / "docs/data/license-inventory-scopes.json").read_bytes()
    )["scopes"]
    for name in (
        "vanadium-james-1960-9012770.cif",
        "superhydrous-phase-b-koch-muller-2005-mkm0105.cif",
    ):
        path = f"peritheos/data/datasets/{name}"
        assert (
            match_scope(path, scopes, False)["status"]
            == "limited_or_unspecified_source"
        )
        current = match_scope(path, scopes, True)
        assert current["status"] == "open_source_declared"
        assert current["declared_license"] == "CC0-1.0 (COD data declaration)"
    he = match_scope(
        "peritheos/data/datasets/brucite-he-2023-table43.cif", scopes, True
    )
    assert he["status"] == "limited_or_unspecified_source"
    assert he["id"] != "cod-cif-current-review"


def test_superseded_he_structure_is_not_a_distribution_input():
    removed = review()["removed_files"][0]
    assert removed["file"] == "brucite-he-2023-table43.cif"
    assert not (DATA / removed["file"]).exists()
    assert removed["distribution_status"] == "removed_after_0.17.4"
    assert removed["replacement"] == "brucite-zigan-1967-cod-1000054.cif"
    source = json.loads((DATA / "brucite-source.json").read_bytes())
    assert "structure_zigan_1967_archive" not in source["sources"]
    assert source["sources"]["structure"]["cod_id"] == "1000054"
    # Source/CIF change must remove copied thesis refinement, not just its filename.
    assert source["sources"]["structure"]["refinement"] == {"R_factor_all": 0.08}
    for path in [
        ROOT / "peritheos/data/materials/brucite.eosmat",
        *(ROOT / "crates/peritheos/tests/data").glob("brucite*.eosmat"),
        *(ROOT / "output/material-candidates/brucite").glob("*.eosmat"),
    ]:
        document = json.loads(path.read_bytes())
        assert document["source"]["structure_cod_id"] == "1000054", path
        assert document["atom_sites"][2]["z"] == 0.4303, path
        assert "structure_atomic_displacements_a2" not in document["source"], path


def test_inventory_retains_historical_provenance_without_thesis_refinement():
    with (ROOT / "docs/data/license-inventory/catalog-components.csv").open(
        encoding="utf-8", newline=""
    ) as f:
        rows = [
            r
            for r in csv.DictReader(f)
            if r["material"] == "brucite" and r["kind"] == "structure_and_diffraction"
        ]
    historical = next(r for r in rows if r["state"] == "release")
    source = json.loads(historical["source"])["source"]
    assert source["structure_reference"]["doi"] == "10.15083/0002012205"
    assert "structure_refinement" not in source
    assert "structure_atomic_displacements_a2" not in source
    current = next(r for r in rows if r["state"] == "local")
    assert current["status"] == "open_source_declared"
    assert current["scope"] == "cod-cif-current-review"
