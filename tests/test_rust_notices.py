"""Dependency notices retain source attribution and exclude dev-only claims."""

import hashlib
import io
import json
import tarfile
from pathlib import Path

import pytest

from scripts.generate_rust_notices import (
    OUTPUT,
    ROOT,
    crate_notices,
    dependency_roles,
    lock_checksums,
)


@pytest.mark.parametrize("ending", [b"\n", b"\r\n"])
def test_locked_checksums_preserve_the_audited_registry_set_with_either_ending(ending):
    lock = (ROOT / "Cargo.lock").read_bytes().replace(b"\r\n", b"\n")
    checksums = lock_checksums(lock.replace(b"\n", ending))
    manifest = json.loads((OUTPUT / "manifest.json").read_bytes())
    assert len(manifest["packages"]) == 47
    for package in manifest["packages"]:
        assert (
            checksums[(package["crate"], package["version"])]
            == package["archive_sha256"]
        )


def test_committed_notices_are_complete_and_hash_bound():
    manifest = json.loads((OUTPUT / "manifest.json").read_bytes())
    expected = {
        "NOTICE.md",
        "manifest.json",
        "stdlib.json",
        manifest["stdlib"]["notice"],
    }
    for package in manifest["packages"]:
        assert package["notice_files"]
        for notice in package["notice_files"]:
            expected.add(notice["path"])
            assert (
                hashlib.sha256((OUTPUT / notice["path"]).read_bytes()).hexdigest()
                == notice["sha256"]
            )
    assert (
        hashlib.sha256((OUTPUT / manifest["stdlib"]["notice"]).read_bytes()).hexdigest()
        == manifest["stdlib"]["sha256"]
    )
    assert {
        p.relative_to(OUTPUT).as_posix() for p in OUTPUT.rglob("*") if p.is_file()
    } == expected
    # These obligations are easy to miss when copying only root LICENSE files.
    assert (OUTPUT / "licenses/rand-0.8.8/COPYRIGHT").is_file()
    assert (
        "Alexander Huszagh"
        in (OUTPUT / "source-notices/serde_json-1.0.151.txt").read_text()
    )
    assert (OUTPUT / "licenses/unicode-ident-1.0.24/LICENSE-UNICODE").is_file()


def test_graph_roles_preserve_build_tools_macros_and_shared_dependencies():
    packages = [
        {
            "id": key,
            "name": "peritheos-python" if key == "root" else key,
            "source": None if key == "root" else "registry",
            "targets": [{"kind": ["proc-macro" if key == "macro" else "lib"]}],
        }
        for key in ("root", "lib", "build", "macro", "shared", "dev")
    ]

    def dep(key, kind=None):
        return {"pkg": key, "dep_kinds": [{"kind": kind}]}

    metadata = {
        "packages": packages,
        "resolve": {
            "nodes": [
                {
                    "id": "root",
                    "deps": [
                        dep("lib"),
                        dep("build", "build"),
                        dep("macro"),
                        dep("dev", "dev"),
                    ],
                },
                {"id": "lib", "deps": [dep("shared")]},
                {"id": "build", "deps": [dep("shared")]},
                {"id": "macro", "deps": []},
                {"id": "shared", "deps": []},
                {"id": "dev", "deps": []},
            ]
        },
    }
    assert dependency_roles(metadata) == {
        "lib": ["target-library"],
        "build": ["host-tool"],
        "macro": ["host-tool"],
        "shared": ["host-tool", "target-library"],
    }


def fixture_crate(tmp_path):
    base = tmp_path / "registry/src/index/toy-1.0"
    files = {
        "LICENSE-MIT": b"Permission and copyright\n",
        "COPYRIGHT": b"Additional author attribution\n",
        "src/lib.rs": b"// Copyright A different author\n// Licensed under MIT\n\npub fn value() {}\n",
    }
    archive = tmp_path / "registry/cache/index/toy-1.0.crate"
    archive.parent.mkdir(parents=True)
    with tarfile.open(archive, "w:gz") as stream:
        for name, content in files.items():
            path = base / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
            info = tarfile.TarInfo("toy-1.0/" + name)
            info.size = len(content)
            stream.addfile(info, io.BytesIO(content))
    package = {
        "name": "toy",
        "version": "1.0",
        "manifest_path": str(base / "Cargo.toml"),
    }
    return package, archive, hashlib.sha256(archive.read_bytes()).hexdigest()


def test_license_and_additional_source_attributions_are_copied(tmp_path):
    package, _, checksum = fixture_crate(tmp_path)
    result = crate_notices(package, checksum)
    assert result["licenses/toy-1.0/COPYRIGHT"] == b"Additional author attribution\n"
    assert b"Copyright A different author" in result["source-notices/toy-1.0.txt"]


@pytest.mark.parametrize("tamper", ["archive", "source"])
def test_changed_archive_or_license_is_rejected(tmp_path, tamper):
    package, archive, checksum = fixture_crate(tmp_path)
    if tamper == "archive":
        archive.write_bytes(archive.read_bytes() + b"altered")
    else:
        (Path(package["manifest_path"]).parent / "LICENSE-MIT").write_bytes(b"altered")
    with pytest.raises(ValueError, match="checksum"):
        crate_notices(package, checksum)
