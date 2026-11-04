"""Generate/check Python-distribution Rust notices from locked registry sources.

Repository maintenance tool. Release wheels use the audited compiler pinned in
publish.yml; other source-build distributors must refresh matching std notices.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "peritheos/data/third_party/rust"
TARGETS = (
    "aarch64-apple-darwin",
    "aarch64-unknown-linux-gnu",
    "x86_64-apple-darwin",
    "x86_64-pc-windows-msvc",
    "x86_64-unknown-linux-gnu",
)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(value):
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode()


def dependency_roles(metadata):
    """Track target libraries separately from host build/procedural-macro tools."""
    packages = {p["id"]: p for p in metadata["packages"]}
    nodes = {n["id"]: n for n in metadata["resolve"]["nodes"]}
    pending = [
        (p["id"], False) for p in packages.values() if p["name"] == "peritheos-python"
    ]
    seen, result = set(), {}
    while pending:
        key, host = pending.pop()
        package = packages[key]
        host = host or any("proc-macro" in t["kind"] for t in package["targets"])
        if (key, host) in seen:
            continue
        seen.add((key, host))
        if package["source"] is not None:
            result.setdefault(key, set()).add("host-tool" if host else "target-library")
        for dep in nodes[key]["deps"]:
            for kind in dep["dep_kinds"]:
                if kind["kind"] in (None, "build"):
                    pending.append((dep["pkg"], host or kind["kind"] == "build"))
    return {key: sorted(value) for key, value in result.items()}


def source_notice_blocks(text):
    """Retain complete source comment blocks containing additional attributions."""
    blocks, current, start, multiline = [], [], 0, False
    for number, line in enumerate(text.splitlines(keepends=True), 1):
        stripped = line.lstrip()
        comment = multiline or stripped.startswith(("//", "/*", "# "))
        if comment or (current and not stripped.strip()):
            if not current:
                start = number
            current.append(line)
            if "/*" in line:
                multiline = True
            if "*/" in line:
                multiline = False
        else:
            if re.search(
                r"copyright|SPDX-License-Identifier|dual licensed|licensed under",
                "".join(current),
                re.I,
            ):
                blocks.append((start, "".join(current)))
            current, multiline = [], False
    if re.search(
        r"copyright|SPDX-License-Identifier|dual licensed|licensed under",
        "".join(current),
        re.I,
    ):
        blocks.append((start, "".join(current)))
    return blocks


def crate_notices(package, expected_checksum):
    base = Path(package["manifest_path"]).parent
    archive_path = (
        base.parent.parent.parent / "cache" / base.parent.name / (base.name + ".crate")
    )
    archive_bytes = archive_path.read_bytes()
    if sha256(archive_bytes) != expected_checksum:
        raise ValueError(
            f"Registry archive checksum differs from Cargo.lock: {package['name']}"
        )
    prefix = base.name + "/"
    files = {}
    with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:gz") as archive:
        for member in archive.getmembers():
            if member.isfile():
                if not member.name.startswith(prefix):
                    raise ValueError(
                        f"Unexpected registry archive member: {member.name}"
                    )
                files[member.name[len(prefix) :]] = archive.extractfile(member).read()
    output, attribution = {}, []
    for name, original in sorted(files.items()):
        path = Path(name)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"Unsafe registry path: {name}")
        is_notice = path.name.upper().startswith(
            ("LICENSE", "COPYING", "NOTICE", "COPYRIGHT")
        )
        if not is_notice and path.suffix.lower() not in (".rs", ".c", ".h", ".s"):
            continue
        content = (base / name).read_bytes()
        if content != original:
            raise ValueError(
                f"Registry source file checksum mismatch: {package['name']}/{name}"
            )
        if is_notice:
            content.decode("utf-8")
            output[
                "licenses/" + package["name"] + "-" + package["version"] + "/" + name
            ] = content
        else:
            for line, block in source_notice_blocks(content.decode("utf-8")):
                attribution.append(f"Source: {name}:{line}\n{block}\n")
    if not output:
        raise ValueError(f"No license files found for {package['name']}")
    if attribution:
        output[f"source-notices/{package['name']}-{package['version']}.txt"] = (
            "Additional source-header notices, retained verbatim below their source locations.\n\n"
            + "".join(attribution)
        ).encode()
    return output


def lock_checksums(content):
    result = {}
    for block in content.decode().split("[[package]]")[1:]:
        fields = dict(re.findall(r'^([a-z_-]+) = "([^"]*)"$', block, re.M))
        if "checksum" in fields:
            result[(fields["name"], fields["version"])] = fields["checksum"]
    return result


def compiler_identity():
    raw = subprocess.check_output(["rustc", "-Vv"], text=True)
    fields = dict(line.split(": ", 1) for line in raw.splitlines() if ": " in line)
    return {"version": fields["release"], "commit": fields["commit-hash"]}


def refresh_stdlib():
    sysroot = Path(
        subprocess.check_output(["rustc", "--print", "sysroot"], text=True).strip()
    )
    source = sysroot / "share/doc/rust/COPYRIGHT-library.html"
    if not source.is_file():
        raise ValueError(
            "Matching Rust COPYRIGHT-library.html is required; install the rust-docs component"
        )
    content = source.read_bytes()
    identity = compiler_identity()
    identity.update(
        {
            "notice": "rust-stdlib/COPYRIGHT-library.html",
            "sha256": sha256(content),
            "source": f"https://doc.rust-lang.org/{identity['version']}/COPYRIGHT-library.html",
            "scope": "Official Rust Standard Library copyright document, including its in-tree and third-party/build dependency notices; a conservative all-target notice set, not the compiler/toolchain being redistributed.",
        }
    )
    OUTPUT.joinpath("rust-stdlib").mkdir(parents=True, exist_ok=True)
    (OUTPUT / identity["notice"]).write_bytes(content)
    (OUTPUT / "stdlib.json").write_bytes(json_bytes(identity))


def metadata_for(target, directory):
    if directory:
        return json.loads((directory / ((target or "all") + ".json")).read_bytes())
    command = ["cargo", "metadata", "--locked", "--format-version", "1"]
    if target:
        command += ["--filter-platform", target]
    return json.loads(subprocess.check_output(command, cwd=ROOT))


def generate(directory=None):
    all_metadata = metadata_for(None, directory)
    packages = {p["id"]: p for p in all_metadata["packages"]}
    all_roles = dependency_roles(all_metadata)
    targets = {t: dependency_roles(metadata_for(t, directory)) for t in TARGETS}
    lock = (ROOT / "Cargo.lock").read_bytes()
    checksums = lock_checksums(lock)
    outputs, records = {}, []
    for key in sorted(
        all_roles, key=lambda k: (packages[k]["name"], packages[k]["version"])
    ):
        package = packages[key]
        checksum = checksums[(package["name"], package["version"])]
        notices = crate_notices(package, checksum)
        outputs.update(notices)
        records.append(
            {
                "crate": package["name"],
                "version": package["version"],
                "license_expression": package["license"],
                "archive_sha256": checksum,
                "source": f"https://crates.io/crates/{package['name']}/{package['version']}",
                "target_roles": {
                    t: roles[key] for t, roles in targets.items() if key in roles
                },
                "notice_files": [
                    {"path": p, "sha256": sha256(b)} for p, b in sorted(notices.items())
                ],
            }
        )
    std = json.loads((OUTPUT / "stdlib.json").read_bytes())
    std_bytes = (OUTPUT / std["notice"]).read_bytes()
    if sha256(std_bytes) != std["sha256"]:
        raise ValueError("Rust standard-library copyright document checksum mismatch")
    outputs[std["notice"]] = std_bytes
    outputs["stdlib.json"] = json_bytes(std)
    manifest = {
        "schema": "peritheos.rust-notices.v1",
        "cargo_lock_lf_sha256": sha256(lock.replace(b"\r\n", b"\n")),
        "targets": list(TARGETS),
        "stdlib": std,
        "packages": records,
        "scope": "Conservative all-target normal/build closure. Target roles distinguish Cargo-selected libraries from host build/procedural-macro tools; selection does not assert every byte survives linking. Additional conditional packages are included for completeness. Dev-only dependencies are excluded.",
        "source_builds": "Official wheels use the audited Rust version. Other-toolchain redistributors must refresh and review matching standard-library notices; source/MSRV builds remain supported without changing their compiler.",
        "external_libraries": "NumPy/SciPy/Python are separate installations; OS libraries are dynamically supplied by the platform. Extra bundled native libraries require an additional audit; installed wheel smoke checks reject unexpected native payloads.",
    }
    outputs["manifest.json"] = json_bytes(manifest)
    lines = [
        "# Rust third-party notices",
        "",
        "Peritheos code remains MIT. Each dependency retains the license expression below;",
        "OR alternatives remain alternatives. Supplied texts include all available options,",
        "COPYRIGHT/NOTICE files, nested notices and additional source-header attributions.",
        "",
        "The notice set conservatively covers all normal/build dependencies, including",
        "procedural macros and conditional packages not selected by the five wheel targets.",
        "Build-tool inclusion does not claim that tool binaries are bundled. `manifest.json`",
        "records exact crate archive checksums, target roles and every notice-file hash.",
        "",
        f"Official wheels use Rust {std['version']} ({std['commit']}).",
        "Its complete standard-library copyright document is in",
        "[COPYRIGHT-library.html](rust-stdlib/COPYRIGHT-library.html). This includes",
        "standard-library dependencies and build inputs without relicensing Peritheos.",
        "Audited Rust source: https://github.com/rust-lang/rust/tree/" + std["version"],
        "Dependency source archives are available at the exact crates.io versions below;",
        "Peritheos does not modify those dependency or standard-library sources.",
        "Other-toolchain redistributors must refresh matching notices before distribution.",
        "",
        "| Crate | Version | Declared license | Wheel-target selection |",
        "| --- | --- | --- | --- |",
    ]
    for r in records:
        selection = (
            "all five targets"
            if len(r["target_roles"]) == len(TARGETS)
            else "conditional extra; none of the five targets"
        )
        lines.append(
            f"| [{r['crate']}]({r['source']}) | {r['version']} | {r['license_expression']} | {selection} |"
        )
    lines += [
        "",
        "Source archive/license files are unchanged copies verified against original crate",
        "archive checksums in Cargo.lock. Additional source comment notices are grouped with",
        "their exact source locations. Do not omit these notices when redistributing the",
        "corresponding compiled code. This is not a scientific-data license grant.",
        "",
    ]
    outputs["NOTICE.md"] = "\n".join(lines).encode()
    return outputs, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--metadata-dir", type=Path)
    parser.add_argument("--refresh-stdlib", action="store_true")
    parser.add_argument("--verify-toolchain", action="store_true")
    args = parser.parse_args()
    if args.refresh_stdlib:
        if args.check:
            parser.error("--refresh-stdlib cannot be combined with --check")
        refresh_stdlib()
    outputs, manifest = generate(args.metadata_dir)
    if args.verify_toolchain and compiler_identity() != {
        k: manifest["stdlib"][k] for k in ("version", "commit")
    }:
        raise SystemExit(
            "Compiler differs from the audited standard-library notice bundle"
        )
    workflow = ROOT / ".github/workflows/publish.yml"
    if args.check and workflow.is_file():
        pin = re.search(r'RELEASE_RUST_TOOLCHAIN:\s*"([^"]+)"', workflow.read_text())
        if not pin or pin[1] != manifest["stdlib"]["version"]:
            raise SystemExit(
                "Publication toolchain pin differs from the audited Rust notices"
            )
    actual = {
        p.relative_to(OUTPUT).as_posix() for p in OUTPUT.rglob("*") if p.is_file()
    }
    if args.check:
        if actual != set(outputs) or any(
            (OUTPUT / p).read_bytes() != b for p, b in outputs.items()
        ):
            raise SystemExit(
                "Stale or modified Rust third-party notices; regenerate and review"
            )
    else:
        for name in actual - set(outputs):
            (OUTPUT / name).unlink()
        for name, content in outputs.items():
            path = OUTPUT / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
    print(
        f"Rust notices: {len(manifest['packages'])} locked packages, {len(TARGETS)} targets, {len(outputs)} files; standard-library document verified."
    )


if __name__ == "__main__":
    main()
