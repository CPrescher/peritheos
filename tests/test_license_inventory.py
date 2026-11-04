"""Rights scopes and archive coverage must not silently become blanket grants."""

import csv
import io
import json
import tarfile
from zipfile import ZipFile

import pytest

from scripts.inventory_distribution_licenses import (
    ROOT,
    assertions,
    classify_file,
    declaration_status,
    dependency_lock_fingerprint,
    digest,
    load_release,
    match_scope,
    source_inventory,
)


def test_release_inventory_allows_workspace_bumps_but_rejects_dependency_drift():
    original = b"""version = 4
[[package]]
name = "peritheos"
version = "0.17.3"
[[package]]
name = "peritheos-python"
version = "0.17.3"
[[package]]
name = "rand"
version = "0.8.8"
source = "registry+https://github.com/rust-lang/crates.io-index"
checksum = "original-checksum"
"""
    bumped = original.replace(b"0.17.3", b"0.17.4")
    assert dependency_lock_fingerprint(original) == dependency_lock_fingerprint(bumped)
    assert dependency_lock_fingerprint(original) != dependency_lock_fingerprint(
        bumped.replace(b"0.8.8", b"0.10.3")
    )
    assert dependency_lock_fingerprint(original) != dependency_lock_fingerprint(
        bumped.replace(b"original-checksum", b"changed-checksum")
    )


@pytest.mark.parametrize(
    "declaration,expected",
    [
        ({"license": "CC-BY-4.0"}, "open_license_declared"),
        ({"license": "CC BY-NC 4.0"}, "restricted_source"),
        (
            {
                "license": "CC0-1.0",
                "license_scope": "Contributor transcription only; source rights unspecified",
            },
            "limited_or_unspecified_source",
        ),
        ({"license": "CC0-1.0"}, "cc0_scope_needs_review"),
        (
            {"notes": "CC BY 4.0 article cited for provenance"},
            "source_terms_not_recorded",
        ),
    ],
)
def test_article_mentions_and_contributor_dedications_are_not_source_grants(
    declaration, expected
):
    assert declaration_status(assertions(declaration)) == expected


def fixture_archives(folder):
    data = b"pressure,volume\n1,2\n"
    entries = {
        "peritheos/data/datasets/example.csv": data,
        "peritheos/data/README.md": b"Source terms remain unspecified.\n",
    }
    path = folder / "peritheos-0.17.3.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name, content in entries.items():
            member = tarfile.TarInfo("peritheos-0.17.3/" + name)
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))
    paths = [(path, "sdist")]
    for label, crlf in (("unix", False), ("windows", True)):
        path = folder / f"peritheos-0.17.3-{label}.whl"
        with ZipFile(path, "w") as archive:
            for name, content in entries.items():
                if crlf and name.endswith("README.md"):
                    content = content.replace(b"\n", b"\r\n")
                archive.writestr(name, content)
            archive.writestr("peritheos/_rust.so", label.encode())
        paths.append((path, "bdist_wheel"))
    metadata = {
        "info": {"version": "0.17.3"},
        "urls": [
            {
                "filename": p.name,
                "packagetype": kind,
                "url": "https://example.invalid/" + p.name,
                "digests": {"sha256": digest(p.read_bytes())},
            }
            for p, kind in paths
        ],
    }
    (folder / "pypi.json").write_text(json.dumps(metadata), encoding="utf-8")
    return paths


def test_inventory_preserves_byte_variants_and_all_archive_members(tmp_path):
    fixture_archives(tmp_path)
    _, archives, _, variants, differences = load_release(tmp_path)
    for archive in archives:
        assert (
            sum(archive["id"] in row["archives"] for row in variants.values())
            == archive["files"]
        )
    data = [r for (p, _), r in variants.items() if p.endswith("example.csv")]
    assert len(data) == 1 and len(data[0]["archives"]) == 3
    readmes = [r for (p, _), r in variants.items() if p.endswith("README.md")]
    assert len(readmes) == 2
    difference = next(iter(differences.values()))
    assert difference["byte_differences"] == ["peritheos/data/README.md"]
    assert difference["differences_after_crlf_normalization"] == []


def test_changed_archive_fails_before_it_can_be_inventoried(tmp_path):
    paths = fixture_archives(tmp_path)
    paths[0][0].write_bytes(b"changed archive")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_release(tmp_path)


def test_local_author_grants_are_not_backdated_to_published_release():
    policy = json.loads(
        (ROOT / "docs/data/license-inventory-scopes.json").read_text(encoding="utf-8")
    )
    path = "peritheos/data/datasets/iron-carbide-fe7c3-prescher-2015-author-input.txt"
    assert match_scope(path, policy["scopes"], False) is None
    assert (
        match_scope(path, policy["scopes"], True)["status"]
        == "author_permission_declared"
    )
    path = "peritheos/data/datasets/morard_2026_author_sources/fes-sata-2010-table1.csv"
    scope = match_scope(path, policy["scopes"], True)
    assert scope is None or scope["status"] not in (
        "author_permission_declared",
        "open_source_declared",
    )
    prefix = "peritheos/data/datasets/ma_2024_sources/"
    assert (
        match_scope(prefix + "supporting-information.docx", policy["scopes"], True)[
            "status"
        ]
        == "restricted_source"
    )
    assert (
        match_scope(prefix + "data-tables.xlsx", policy["scopes"], True)["status"]
        == "open_source_declared"
    )


def test_project_code_license_does_not_cover_unknown_original_author_code():
    status = classify_file(
        "peritheos/data/datasets/unknown-author.py", b"", [], {}, True
    )[0]
    assert status == "source_terms_not_recorded"
    assert (
        classify_file("peritheos/materials.py", b"", [], {}, True)[0] == "project_mit"
    )


def test_source_index_accepts_literal_and_structured_citations():
    rows = [
        {
            "state": "local",
            "source": json.dumps(reference),
            "material": "example",
            "kind": "dataset",
            "status": "source_terms_not_recorded",
            "resource": "",
            "evidence": "example#/datasets/0",
        }
        for reference in ("Author2010", {"doi": "10.example/data"}, {})
    ]
    assert len(source_inventory(rows)) == 3


def test_committed_inventory_covers_release_members_and_inline_components():
    folder = ROOT / "docs/data/license-inventory"
    summary = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
    with (folder / "files.csv").open(encoding="utf-8", newline="") as stream:
        files = list(csv.DictReader(stream))
    for archive in summary["archives"]:
        assert (
            sum(archive["id"] in r["archive_ids"].split("|") for r in files)
            == archive["files"]
        )
    with (folder / "catalog-components.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        components = list(csv.DictReader(stream))
    for state in ("release", "local"):
        assert (
            sum(r["state"] == state and r["kind"] == "eos" for r in components) == 685
        )
        assert (
            sum(
                r["state"] == state and r["kind"] == "structure_and_diffraction"
                for r in components
            )
            == 234
        )
        assert (
            sum(r["state"] == state and r["kind"] == "dataset" for r in components)
            == 449
        )
    assert summary["published_wheels_checked"] == 35
    assert not summary["all_wheel_data_bytes_identical"]
    assert summary["all_wheel_data_identical_after_crlf_normalization"]
