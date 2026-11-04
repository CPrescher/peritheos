"""CIF source grants must be file-scoped and leave structural evidence intact."""

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
        (
            "brucite-he-2023-table43.cif",
            "9983ea18cf364ef6b7378b44b0b0a20e28b29c4d08b6c76399e9d79605b915ff",
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
