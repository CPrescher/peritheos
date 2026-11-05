#!/usr/bin/env python3
"""Repository maintenance: inventory immutable release archives and local rights.

This records evidence and unresolved questions. It neither grants permissions
nor selects a redistributable subset. No scientific source is rewritten.
"""

from __future__ import annotations

import argparse
import ast
import csv
import fnmatch
import hashlib
import io
import json
import re
import subprocess
import tarfile
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "docs/data/license-inventory-scopes.json"
OUTPUT = ROOT / "docs/data/license-inventory"
OPEN_LICENSES = {"CC-BY-4.0", "CC-BY-SA-4.0", "MIT", "Apache-2.0"}
RIGHTS_KEYS = {
    "license",
    "licence",
    "source_license",
    "data_license",
    "source_data_license",
    "source_rights",
    "contributor_rights",
    "license_scope",
    "copyright_note",
    "rights",
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def dependency_lock_fingerprint(content):
    """Allow workspace version bumps without overlooking dependency changes."""
    blocks = content.replace(b"\r\n", b"\n").split(b"[[package]]")
    for index, block in enumerate(blocks):
        if re.search(
            rb'(?m)^name = "(?:peritheos|peritheos-python)"$', block
        ) and not re.search(rb"(?m)^source = ", block):
            blocks[index] = re.sub(
                rb'(?m)^version = "[^"\n]+"$',
                b'version = "<workspace>"',
                block,
                count=1,
            )
    return digest(b"[[package]]".join(blocks))


def compact(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def assertions(value, pointer=""):
    """Retain scoped statements verbatim; a CC token alone is not a grant."""
    result = []
    if isinstance(value, dict):
        for key, item in value.items():
            location = pointer + "/" + key.replace("~", "~0").replace("/", "~1")
            if key in RIGHTS_KEYS:
                result.append({"pointer": location, "statement": item})
            elif isinstance(item, (dict, list)):
                result.extend(assertions(item, location))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            result.extend(assertions(item, pointer + "/" + str(index)))
    return result


def declaration_status(items):
    """Conservative triage of declarations, never a legal clearance verdict."""
    text = compact(items)
    if (
        "CC BY-NC" in text
        or "CC-BY-NC" in text
        or "CC BY-ND" in text
        or "CC-BY-ND" in text
    ):
        return "restricted_source"
    # Check limitations before recognizing a bare license identifier.
    if any(
        word in text.lower()
        for word in (
            "contributor",
            "unspecified",
            "not state",
            "not identified",
            "no separate",
            "no explicit",
        )
    ):
        return "limited_or_unspecified_source"
    for item in items:
        token = (
            str(item["statement"])
            .replace("CC BY-SA ", "CC-BY-SA-")
            .replace("CC BY ", "CC-BY-")
        )
        if token in OPEN_LICENSES:
            return "open_license_declared"
        if token == "CC0-1.0":
            return "cc0_scope_needs_review"
    return "source_terms_not_recorded" if not items else "limited_or_unspecified_source"


def load_release(directory):
    metadata = json.loads((directory / "pypi.json").read_text(encoding="utf-8"))
    archives, payloads, variants = [], {}, {}
    wheel_data = None
    normalized_wheel_data = None
    wheel_differences = {}
    wheels = 0
    for index, info in enumerate(sorted(metadata["urls"], key=lambda x: x["filename"])):
        path = directory / info["filename"]
        actual = digest(path.read_bytes())
        if actual != info["digests"]["sha256"]:
            raise ValueError(f"Archive checksum mismatch: {path.name}")
        archive_id = "a" + str(index + 1).zfill(2)
        kind = info["packagetype"]
        archives.append(
            {
                "id": archive_id,
                "filename": path.name,
                "kind": kind,
                "sha256": actual,
                "url": info["url"],
                "size": path.stat().st_size,
            }
        )
        if kind == "sdist":
            with tarfile.open(path) as archive:
                entries = [
                    (m.name.partition("/")[2], archive.extractfile(m).read())
                    for m in archive.getmembers()
                    if m.isfile()
                ]
        else:
            with zipfile.ZipFile(path) as archive:
                entries = [
                    (n, archive.read(n))
                    for n in archive.namelist()
                    if not n.endswith("/")
                ]
            data = {n: digest(b) for n, b in entries if n.startswith("peritheos/data/")}
            normalized = {
                n: digest(b.replace(b"\r\n", b"\n"))
                for n, b in entries
                if n.startswith("peritheos/data/")
            }
            if wheel_data is None:
                wheel_data = data
                normalized_wheel_data = normalized
            elif data != wheel_data:
                wheel_differences[path.name] = {
                    "byte_differences": sorted(
                        n
                        for n in data.keys() | wheel_data.keys()
                        if data.get(n) != wheel_data.get(n)
                    ),
                    "differences_after_crlf_normalization": sorted(
                        n
                        for n in normalized.keys() | normalized_wheel_data.keys()
                        if normalized.get(n) != normalized_wheel_data.get(n)
                    ),
                }
            wheels += 1
        archives[-1]["files"] = len(entries)
        for name, content in entries:
            checksum = digest(content)
            key = (name, checksum)
            row = variants.setdefault(
                key,
                {
                    "path": name,
                    "sha256": checksum,
                    "bytes": len(content),
                    "archives": [],
                },
            )
            row["archives"].append(archive_id)
            # Keep one copy of shared source payloads, not every native binary.
            if not name.endswith((".so", ".pyd", ".dylib")):
                payloads.setdefault(name, content)
    if wheels != sum(x["packagetype"] == "bdist_wheel" for x in metadata["urls"]):
        raise ValueError("Incomplete wheel inventory")
    return metadata, archives, payloads, variants, wheel_differences


def current_payloads(released):
    result = {
        name: (ROOT / name).read_bytes() if (ROOT / name).is_file() else content
        for name, content in released.items()
    }
    paths = subprocess.check_output(
        [
            "git",
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "--",
            "peritheos/data",
        ],
        cwd=ROOT,
        text=True,
    ).splitlines()
    for name in paths:
        path = ROOT / name
        if path.is_file():
            result[name] = path.read_bytes()
    return result


def match_scope(path, scopes, current):
    for scope in scopes:
        if scope.get("local_only") and not current:
            continue
        if any(fnmatch.fnmatchcase(path, pattern) for pattern in scope["files"]):
            return scope
    return None


def components(payloads, scopes, current):
    rows = []
    for name, content in sorted(payloads.items()):
        if not name.startswith("peritheos/data/materials/") or not name.endswith(
            ".eosmat"
        ):
            continue
        document = json.loads(content)
        material = document["identifier"]
        for kind, key in (("dataset", "datasets"), ("eos", "eos_records")):
            for index, item in enumerate(document.get(key, [])):
                location = f"{name}#/{key}/{index}"
                rights = assertions(item)
                resource = item.get("resource", {}).get("path", "")
                resource = "peritheos/data/" + resource if resource else ""
                scope = match_scope(resource, scopes, current) if resource else None
                status = scope["status"] if scope else declaration_status(rights)
                rows.append(
                    {
                        "kind": kind,
                        "material": material,
                        "identifier": item["identifier"],
                        "resource": resource,
                        "evidence": location,
                        "status": status,
                        "scope": scope["id"] if scope else "",
                        "content_kind": item.get("kind", kind),
                        "scope_note": scope["scope"]
                        if scope
                        else "Existing scoped statements are retained verbatim. A source license does not automatically restrict uncopyrightable facts or contributor calculations, or grant rights in third-party content.",
                        "related_rights_evidence": "",
                        "source": compact(item.get("reference", {})),
                        "rights_statements": compact(rights),
                        "required_action": scope["action"] if scope else action(status),
                    }
                )
        # A material file is not just its datasets: do not omit structures/peaks.
        structural_source = {
            k: document[k] for k in ("reference", "source", "notes") if k in document
        }
        if material == "brucite" and not current:
            # Inventory the historical thesis source without redistributing its
            # superseded refinement through this generated report. Exact source
            # archives/hashes and the complete original inventory are preserved.
            structural_source["source"] = {
                k: v
                for k, v in document.get("source", {}).items()
                if k
                in (
                    "structure_reference",
                    "structure_access_url",
                    "structure_source_kind",
                    "structure_dataset",
                    "structure_dataset_sha256",
                    "source_manifest",
                )
            }
        rows.append(
            {
                "kind": "structure_and_diffraction",
                "material": material,
                "identifier": material,
                "resource": "",
                "evidence": name + "#/",
                "status": "source_terms_not_recorded",
                "scope": "",
                "content_kind": "structural facts and/or imported diffraction data",
                "scope_note": "No blanket license is inferred from the material's EOS or dataset licenses.",
                "related_rights_evidence": "",
                "source": compact(structural_source),
                "rights_statements": compact(
                    {k: document[k] for k in RIGHTS_KEYS if k in document}
                ),
                "required_action": "Assess structural facts, source compilation rights and any imported peak patterns separately; a dataset license is not a blanket grant for the material.",
            }
        )
        if material == "brucite":
            row = rows[-1]
            if not current:
                row["scope_note"] += (
                    " Superseded He refinement values are omitted from this"
                    " inventory; original values/inventory preserved privately."
                )
            elif document.get("source", {}).get(
                "structure_cod_id"
            ) == "1000054" and not document.get("peaks"):
                resource = "peritheos/data/" + document["source"]["structure_dataset"]
                scope = match_scope(resource, scopes, current)
                if scope and scope["id"] == "cod-cif-current-review":
                    row.update(
                        resource=resource,
                        status=scope["status"],
                        scope=scope["id"],
                        scope_note="CC0 COD 1000054 structural reference only; no imported peak pattern. Original paper not directly inspected; EOS/data source terms remain separate.",
                        rights_statements=compact(
                            {
                                "structure_license": document["source"][
                                    "structure_license"
                                ]
                            }
                        ),
                        required_action=scope["action"],
                    )
    by_doi = defaultdict(list)
    for row in rows:
        ref = json.loads(row["source"])
        if row["kind"] == "dataset" and isinstance(ref, dict) and ref.get("doi"):
            by_doi[ref["doi"]].append(
                {
                    "evidence": row["evidence"],
                    "status": row["status"],
                    "scope": row["scope"],
                }
            )
    for row in rows:
        ref = json.loads(row["source"])
        if row["kind"] == "eos" and isinstance(ref, dict):
            row["related_rights_evidence"] = compact(by_doi.get(ref.get("doi"), []))
    return rows


def action(status):
    if status in (
        "open_license_declared",
        "open_source_declared",
        "third_party_code_open",
    ):
        return "Preserve attribution, full license and changes; verify the declared grant's scope and third-party exclusions before downstream clearance."
    if status == "restricted_source":
        return "Obtain a suitable separate grant or assess factual-data applicability; otherwise exclude from a Debian-main source archive, retaining private evidence."
    if status == "mixed_author_permission":
        return "Author-held rights are CC BY 4.0; assess embedded literature inputs separately. Retain the maintainer's permission basis and attribution."
    if status == "author_permission_declared":
        return "Preserve author attribution, CC BY 4.0 and change notices; retain the dated maintainer declaration and third-party exclusions."
    if status == "project_mit":
        return "Retain the project MIT notice for contributor-written code/text; preserve source attribution and assess any embedded third-party source content separately."
    if status in ("license_notice", "build_metadata"):
        return "Retain the notice or metadata; its presence does not license unrelated resources."
    if status == "native_binary":
        return "Audit the normal/build Rust dependency closure and required third-party notices; the project MIT notice alone is not a dependency license inventory."
    return "Find an applicable source grant or obtain a documented rights/factual-data assessment; otherwise exclude from a cleared source archive. Contributor-only CC0 does not clear third-party rights."


def classify_file(path, content, scopes, resource_rows, current):
    if path.startswith("peritheos/data/third_party/rust/"):
        return (
            "license_notice",
            "rust-notice-bundle",
            "peritheos/data/third_party/rust/manifest.json",
            action("license_notice"),
        )
    if (
        path.endswith(".LICENSE.md")
        or path.endswith("/LICENSE.md")
        or path in ("LICENSE", "peritheos/data/LICENSING.md")
        or path.startswith("LICENSES/")
        or ".dist-info/licenses/" in path
    ):
        return "license_notice", "", path, action("license_notice")
    scope = match_scope(path, scopes, current)
    if scope:
        return (
            scope["status"],
            scope["id"],
            " | ".join(scope["evidence"]),
            scope["action"],
        )
    if path in resource_rows:
        items = resource_rows[path]
        statuses = {r["status"] for r in items}
        status = (
            next(iter(statuses))
            if len(statuses) == 1
            else "conflicting_or_mixed_declarations"
        )
        return status, "", " | ".join(r["evidence"] for r in items), action(status)
    if path.endswith((".so", ".pyd", ".dylib")):
        return (
            "native_binary",
            "",
            "Cargo.lock; cargo-components.csv",
            action("native_binary"),
        )
    if ".dist-info/" in path or path == "PKG-INFO":
        return "build_metadata", "", path, action("build_metadata")
    project_code = (
        path.startswith("scripts/")
        or path.startswith("peritheos/")
        and not path.startswith("peritheos/data/")
        or path.startswith("crates/")
        and "/tests/data/" not in path
        or path
        in (
            "Cargo.toml",
            "Cargo.lock",
            "pyproject.toml",
            "mkdocs.yml",
            "docs/javascripts/mathjax.js",
            "peritheos/data/__init__.py",
        )
        or path.startswith("peritheos/data/")
        and path.endswith(".schema.json")
    )
    if project_code:
        return (
            "project_mit",
            "",
            "LICENSE; peritheos/data/LICENSING.md",
            action("project_mit"),
        )
    if path.startswith(
        (
            "peritheos/data/",
            "docs/data/",
            "docs/reconstruction-reports/",
            "docs/material-summaries/",
            "docs/literature-reproductions/",
            "docs/research/",
            "output/",
            "crates/peritheos/tests/data/",
        )
    ):
        status = (
            "mixed_scientific_content"
            if path.endswith((".eosmat", ".md", ".json"))
            or path.startswith(
                (
                    "docs/reconstruction-reports/",
                    "docs/material-summaries/",
                    "docs/literature-reproductions/",
                    "docs/data/",
                )
            )
            and path.endswith((".svg", ".png", ".pdf"))
            else "source_terms_not_recorded"
        )
        return status, "", path + "; peritheos/data/LICENSING.md", action(status)
    if path.startswith(("docs/notebooks/", "docs/images/")):
        return "mixed_scientific_content", "", path, action("mixed_scientific_content")
    if path.endswith((".md", ".cff")):
        return "project_mit", "", "LICENSE; " + path, action("project_mit")
    return "source_terms_not_recorded", "", path, action("source_terms_not_recorded")


def cargo_inventory(metadata_path):
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    packages = {p["id"]: p for p in metadata["packages"]}
    nodes = {n["id"]: n for n in metadata["resolve"]["nodes"]}
    roots = [p["id"] for p in metadata["packages"] if p["name"] == "peritheos-python"]
    seen, pending = set(), list(roots)
    while pending:
        key = pending.pop()
        if key in seen:
            continue
        seen.add(key)
        for dep in nodes[key]["deps"]:
            if any(k["kind"] in (None, "build") for k in dep["dep_kinds"]):
                pending.append(dep["pkg"])
    rows = []
    for key in sorted(seen):
        p = packages[key]
        if p["source"] is None:
            continue
        base = Path(p["manifest_path"]).parent
        texts = []
        for file in sorted(base.iterdir()):
            if file.is_file() and file.name.upper().startswith(
                ("LICENSE", "COPYING", "NOTICE", "COPYRIGHT")
            ):
                texts.append({"file": file.name, "sha256": digest(file.read_bytes())})
        rows.append(
            {
                "crate": p["name"],
                "version": p["version"],
                "license_expression": p.get("license") or "unspecified",
                "source": f"https://crates.io/crates/{p['name']}/{p['version']}",
                "license_texts": compact(texts),
                "scope": "All-target normal/build dependency closure; includes build tools, excludes dev-only dependencies. Not a per-wheel linked-code determination.",
                "required_action": "Preserve the packaged Rust notice bundle and standard-library copyright document. Target roles and archive checksums are recorded in peritheos/data/third_party/rust/manifest.json; re-audit dependency, compiler or native-payload changes. Crate source code is not vendored in the source archive.",
            }
        )
    return rows


def csv_bytes(rows, fields):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def source_inventory(rows):
    groups = {}
    for row in rows:
        reference = json.loads(row["source"])
        key = (
            (reference.get("doi") or reference.get("title"))
            if isinstance(reference, dict)
            else reference
        )
        key = key or "unresolved-source:" + row["material"]
        entry = groups.setdefault(
            (row["state"], key),
            {
                "state": row["state"],
                "source": key,
                "citation": compact(reference),
                "component_count": 0,
                "kinds": set(),
                "materials": set(),
                "statuses": set(),
                "resources": set(),
                "evidence": set(),
            },
        )
        entry["component_count"] += 1
        for key, value in (
            ("kinds", row["kind"]),
            ("materials", row["material"]),
            ("statuses", row["status"]),
            ("resources", row["resource"]),
            ("evidence", row["evidence"]),
        ):
            if value:
                entry[key].add(value)
    result = []
    for _, row in sorted(groups.items()):
        result.append(
            {
                key: " | ".join(sorted(value)) if isinstance(value, set) else value
                for key, value in row.items()
            }
        )
    return result


def dataset_source_report(rows):
    """Readable publication index; exact file membership stays in the CSVs."""
    sources = source_inventory(
        [r for r in rows if r["state"] == "local" and r["kind"] == "dataset"]
    )
    text = [
        "# Dataset source-rights index",
        "",
        "Generated from all 449 catalog dataset entries. Counts are entries, not independent experiments or unique physical files. Source grouping does not extend a dataset grant to an entire publication. See [the inventory review](../../research/license-inventory.md) for status meanings and limitations.",
        "",
        "| Source | Dataset entries | Materials | Recorded statuses |",
        "| --- | ---: | --- | --- |",
    ]
    for row in sources:
        citation = json.loads(row["citation"])
        label = row["source"]
        if isinstance(citation, dict):
            authors = citation.get("authors", [])
            first = authors[0] if isinstance(authors, list) and authors else ""
            label = f"{first} {citation.get('year', '')}: {label}".strip()
            if citation.get("doi"):
                label = f"[{label}](https://doi.org/{citation['doi']})"
        text.append(
            "| "
            + " | ".join(
                str(x).replace("|", ";").replace("\n", " ")
                for x in (
                    label,
                    row["component_count"],
                    row["materials"],
                    row["statuses"],
                )
            )
            + " |"
        )
    return ("\n".join(text) + "\n").encode("utf-8")


def download_recorded_archives(directory, receipt_path):
    """Rehydrate the fixed receipt, rather than silently changing release scope."""
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    directory.mkdir(parents=True, exist_ok=True)
    urls = []
    for archive in receipt["archives"]:
        path = directory / archive["filename"]
        if not path.is_file():
            urllib.request.urlretrieve(archive["url"], path)
        if digest(path.read_bytes()) != archive["sha256"]:
            raise ValueError(f"Archive checksum mismatch: {path.name}")
        urls.append(
            {
                "filename": path.name,
                "packagetype": archive["kind"],
                "url": archive["url"],
                "digests": {"sha256": archive["sha256"]},
            }
        )
    metadata = {
        "info": {
            "version": receipt["release_version"],
            "requires_dist": receipt["python_requires_dist_not_bundled"],
        },
        "urls": urls,
    }
    metadata_path = directory / "pypi.json"
    if not metadata_path.is_file():
        metadata_path.write_bytes(
            (json.dumps(metadata, indent=2) + "\n").encode("utf-8")
        )
    else:
        existing = json.loads(metadata_path.read_text(encoding="utf-8"))
        expected = {(u["filename"], u["digests"]["sha256"]) for u in urls}
        actual = {(u["filename"], u["digests"]["sha256"]) for u in existing["urls"]}
        if actual != expected:
            raise ValueError(
                "Existing PyPI metadata differs from the fixed archive receipt"
            )


def inventory(directory, cargo_path):
    policy = json.loads(POLICY.read_text(encoding="utf-8"))
    scopes = policy["scopes"]
    config = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    exclude_block = re.search(r"(?ms)^exclude = (\[.*?^\])", config)
    exclusions = ast.literal_eval(exclude_block.group(1)) if exclude_block else []
    for path in policy["local_distribution_exclusions"]:
        if not any(fnmatch.fnmatchcase(path, pattern) for pattern in exclusions):
            raise ValueError(
                f"Inventory exclusion not present in build configuration: {path}"
            )
    metadata, archives, released, variants, wheel_differences = load_release(directory)
    if dependency_lock_fingerprint(
        released["Cargo.lock"]
    ) != dependency_lock_fingerprint((ROOT / "Cargo.lock").read_bytes()):
        raise ValueError(
            "Current dependency lock differs from the inventoried release; provide release-matched dependency metadata"
        )
    current = current_payloads(released)
    old_components = components(released, scopes, False)
    new_components = components(current, scopes, True)
    resources = []
    for state, items in (("release", old_components), ("local", new_components)):
        mapping = defaultdict(list)
        for item in items:
            if item["resource"]:
                mapping[item["resource"]].append(item)
        resources.append(mapping)
    files = []
    for (name, checksum), row in sorted(variants.items()):
        before = classify_file(
            name, released.get(name, b""), scopes, resources[0], False
        )
        after = classify_file(name, current.get(name, b""), scopes, resources[1], True)
        local = ROOT / name
        scope = match_scope(name, scopes, True)
        if after[0] in ("license_notice", "build_metadata"):
            terms = "Notice/metadata only"
        elif scope:
            terms = scope["declared_license"]
        elif after[0] == "project_mit":
            terms = "MIT"
        else:
            terms = (
                " | ".join(r["rights_statements"] for r in resources[1].get(name, []))
                or "See component/source scopes; no blanket source license established"
            )
        files.append(
            {
                "path": name,
                "release_sha256": checksum,
                "bytes": row["bytes"],
                "archive_ids": "|".join(row["archives"]),
                "release_status": before[0],
                "local_status": after[0],
                "declared_terms": terms,
                "scope": after[1] or before[1],
                "evidence": after[2],
                "required_action": after[3],
                "local_file_sha256": digest(local.read_bytes())
                if name in current and local.is_file()
                else "",
                "local_distribution_excluded": str(
                    name in policy["local_distribution_exclusions"]
                ).lower(),
            }
        )
    # New local data notices are not silently presented as part of 0.17.3.
    for name in sorted(set(current) - set(released)):
        if name.startswith("peritheos/data/"):
            status, scope, evidence, todo = classify_file(
                name, current[name], scopes, resources[1], True
            )
            files.append(
                {
                    "path": name,
                    "release_sha256": "",
                    "bytes": len(current[name]),
                    "archive_ids": "",
                    "release_status": "not_in_release",
                    "local_status": status,
                    "declared_terms": "Notice/metadata only"
                    if status == "license_notice"
                    else "See scoped statements",
                    "scope": scope,
                    "evidence": evidence,
                    "required_action": todo,
                    "local_file_sha256": digest(current[name]),
                    "local_distribution_excluded": "false",
                }
            )
    for scope in scopes:
        if not any(
            any(fnmatch.fnmatchcase(name, p) for p in scope["files"])
            for name in current
        ):
            raise ValueError(f"Scope matches no inventoried file: {scope['id']}")
    # Every archive member must be accounted for, including metadata and copies.
    for archive in archives:
        covered = sum(archive["id"] in row["archive_ids"].split("|") for row in files)
        if covered != archive["files"]:
            raise ValueError(f"Incomplete file inventory: {archive['filename']}")
    rights = []
    for state, payloads in (("release", released), ("local", current)):
        for name, content in sorted(payloads.items()):
            if not name.endswith((".json", ".eosmat")):
                continue
            try:
                value = json.loads(content)
            except (ValueError, UnicodeDecodeError):
                continue
            for item in assertions(value):
                rights.append(
                    {
                        "state": state,
                        "file": name,
                        "pointer": item["pointer"],
                        "statement": compact(item["statement"]),
                    }
                )
    cargo = cargo_inventory(cargo_path)
    summary = {
        "schema": "peritheos.license-inventory.v1",
        "reviewed_on": policy["reviewed_on"],
        "release_version": metadata["info"]["version"],
        "overlay_package_version": next(
            node.value.value
            for node in ast.parse(
                (ROOT / "peritheos/__init__.py").read_text(encoding="utf-8")
            ).body
            if isinstance(node, ast.Assign)
            and any(
                getattr(target, "id", None) == "__version__" for target in node.targets
            )
        ),
        "archives": archives,
        "published_wheels_checked": sum(a["kind"] == "bdist_wheel" for a in archives),
        "all_wheel_data_bytes_identical": not wheel_differences,
        "wheel_data_differences_from_first_wheel": wheel_differences,
        "all_wheel_data_identical_after_crlf_normalization": all(
            not x["differences_after_crlf_normalization"]
            for x in wheel_differences.values()
        ),
        "release_file_variants": len(variants),
        "release_paths": len({p for p, _ in variants}),
        "release_components": dict(Counter(r["kind"] for r in old_components)),
        "local_dataset_statuses": dict(
            sorted(
                Counter(
                    r["status"] for r in new_components if r["kind"] == "dataset"
                ).items()
            )
        ),
        "local_file_statuses": dict(
            sorted(Counter(r["local_status"] for r in files).items())
        ),
        "rust_dependency_components": len(cargo),
        "python_requires_dist_not_bundled": metadata["info"].get("requires_dist", []),
        "limits": policy["limits"],
    }
    combined = [
        {"state": state, **r}
        for state, rows in (("release", old_components), ("local", new_components))
        for r in rows
    ]
    outputs = {
        "files.csv": csv_bytes(files, list(files[0])),
        "catalog-components.csv": csv_bytes(combined, list(combined[0])),
        "rights-statements.csv": csv_bytes(rights, list(rights[0])),
        "cargo-components.csv": csv_bytes(cargo, list(cargo[0])),
        "sources.csv": csv_bytes(
            source_inventory(combined), list(source_inventory(combined)[0])
        ),
        "dataset-sources.md": dataset_source_report(combined),
        "summary.json": (
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
        ).encode("utf-8"),
    }
    return outputs, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archives",
        type=Path,
        required=True,
        help="Verified PyPI archives plus their original pypi.json metadata",
    )
    parser.add_argument(
        "--cargo-metadata",
        type=Path,
        required=True,
        help="cargo metadata --locked --format-version 1 output",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download the exact archives in the committed output/summary.json receipt, verifying their hashes",
    )
    args = parser.parse_args()
    if args.download:
        download_recorded_archives(args.archives, args.output / "summary.json")
    outputs, summary = inventory(args.archives, args.cargo_metadata)
    for name, content in outputs.items():
        path = args.output / name
        if args.check:
            if not path.is_file() or path.read_bytes() != content:
                raise SystemExit(f"Stale license inventory: {path}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
    print(
        f"Inventoried {len(summary['archives'])} release archives, {summary['release_paths']} paths and {sum(summary['release_components'].values())} catalog components; no scientific sources modified."
    )


if __name__ == "__main__":
    main()
