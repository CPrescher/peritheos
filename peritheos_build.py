"""PEP 517 build hooks and explicit distribution exports (standard library only).

Repository EOSMAT and evidence stay complete. Published builds use a temporary
runtime export; editable installs retain the source checkout for curation.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import copy
import csv
import hashlib
import io
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def version(root=ROOT):
    tree = ast.parse((root / "peritheos/__init__.py").read_text(encoding="utf-8"))
    return next(
        node.value.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(getattr(t, "id", None) == "__version__" for t in node.targets)
    )


def digest(value):
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode()
    ).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(
        (
            json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        ).encode("utf-8")
    )


def _summary(value):
    """Retain prose/scalar findings, excluding raw numeric sequences/checkpoints.

    Lists of strings retain source locations, constraints and qualifications.
    Model parameter arrays are outside evidence and never pass through here.
    """
    if isinstance(value, dict):
        return {
            k: _summary(v)
            for k, v in value.items()
            if not isinstance(v, list) or all(isinstance(x, str) for x in v)
        }
    return copy.deepcopy(value)


def runtime_document(original):
    """Explicit v4 runtime export, preserving the executable scientific model."""
    if original.get("format_version") != 4:
        raise ValueError("Distribution export requires native EOSMAT v4")
    doc = copy.deepcopy(original)
    datasets = doc.get("datasets", [])
    for dataset in datasets:
        if "rows" in dataset:
            # Keep a normal checksum-bound external declaration, not an empty
            # table that could be mistaken for a source reporting zero rows.
            rows = dataset.pop("rows")
            stream = io.StringIO(newline="")
            writer = csv.writer(stream)
            writer.writerow([c["name"] for c in dataset["columns"]])
            writer.writerows(rows)
            dataset["resource"] = {
                "path": f"datasets/inline/{doc['identifier']}-{dataset['identifier']}.csv",
                "sha256": hashlib.sha256(stream.getvalue().encode()).hexdigest(),
                "media_type": "text/csv",
            }
        dataset["distribution"] = {
            "storage": "optional",
            "package": "peritheos-datasets",
            "note": "Declaration only; observations are not in the runtime package. The original source declaration is retained by the optional data package where admitted.",
        }
    first = (
        doc["eos_records"][0]["identifier"] if doc["eos_records"] else doc["identifier"]
    )
    url = "https://peritheos.com/eos/" + first
    doc["data_availability"] = {
        "package": "peritheos-datasets",
        "version": version(),
        "dataset_identifiers": [d["identifier"] for d in datasets],
        "source_document_sha256": digest(original),
        "url": url,
        "note": "Experimental observations and detailed reconstruction inputs are separate. The optional package contains only explicitly selected source scopes; see the source record for other evidence.",
    }
    for record in doc["eos_records"]:
        evidence = record.get("evidence", {})
        validation = evidence.get("scientific_validation", {})
        validation.pop("numerical_reproduction", None)
        reproduction = evidence.get("reproduction", {})
        assessment = reproduction.get("assessment")
        evidence["reproduction"] = {
            "report_url": url,
            **({"assessment": _summary(assessment)} if assessment else {}),
        }
        record["evidence"] = _summary(evidence)
    if "source" in doc:
        doc["source"] = {
            k: _summary(v)
            for k, v in doc["source"].items()
            if k.startswith(("structure_", "pattern_"))
            or k in ("reference", "url", "doi")
        }
    return doc


def runtime_study_document(original):
    """Keep public study context and declarations, excluding private audit history."""
    keys = (
        "identifier",
        "record_kind",
        "formula",
        "reference",
        "determination_method",
        "model_description",
        "model_description_source",
        "temperature_k",
        "temperature_role",
        "phase",
        "phase_status",
        "volume_normalization",
        "reference_pressure_gpa",
        "reference_volume",
        "parameter_uncertainties",
        "simulation_pressure_range_gpa",
        "simulation_volume_range",
        "parameter_uncertainties_status",
        "experimental_comparison",
        "comparisons",
        "reproduction",
        "redistribution_review",
    )
    doc = _summary({k: original[k] for k in keys if k in original})
    doc["datasets"] = copy.deepcopy(original["datasets"])
    for dataset in doc["datasets"]:
        if "rows" in dataset:
            raise ValueError("Study export currently requires checksum-bound resources")
        dataset["distribution"] = {
            "storage": "optional",
            "package": "peritheos-datasets",
            "note": "Declaration only; numerical values are optional.",
        }
    doc["source_url"] = "https://doi.org/" + original["reference"]["doi"]
    return doc


def stage_core(destination, root=ROOT):
    """Create a standalone source tree; never mutate the scientific checkout."""
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError(f"Export destination already exists: {destination}")
    destination.mkdir(parents=True)
    for name in (
        "Cargo.toml",
        "Cargo.lock",
        "LICENSE",
        "README.md",
        "CHANGELOG.md",
        "CITATION.cff",
        "CODE_OF_CONDUCT.md",
        "SECURITY.md",
        "SUPPORT.md",
    ):
        shutil.copy2(root / name, destination / name)
    shutil.copytree(root / "LICENSES", destination / "LICENSES")
    for crate in (root / "crates").iterdir():
        if not crate.is_dir():
            continue
        target = destination / "crates" / crate.name
        target.mkdir(parents=True)
        for path in crate.iterdir():
            if path.name == "src":
                shutil.copytree(path, target / path.name)
            elif path.is_file() and path.suffix in (".toml", ".md"):
                shutil.copy2(path, target / path.name)
    shutil.copytree(
        root / "peritheos",
        destination / "peritheos",
        ignore=shutil.ignore_patterns("data", "__pycache__", "*.so", "*.pyd"),
    )
    data = destination / "peritheos/data"
    data.mkdir()
    for path in (root / "peritheos/data").iterdir():
        if path.is_file() and (
            path.suffix == ".json" or path.name in ("__init__.py", "LICENSING.md")
        ):
            if path.name not in ("primary-source-audit.json", "observation-views.json"):
                shutil.copy2(path, data / path.name)
    shutil.copy2(root / "packages/core-LICENSING.md", data / "LICENSING.md")
    registry = json.loads(
        (root / "peritheos/data/observation-views.json").read_text(encoding="utf-8")
    )
    registry["datasets"] = []
    write_json(data / "observation-views.json", registry)
    shutil.copytree(root / "peritheos/data/third_party", data / "third_party")
    materials = data / "materials"
    materials.mkdir()
    for path in (root / "peritheos/data/materials").iterdir():
        if path.suffix == ".eosmat":
            original = json.loads(path.read_text(encoding="utf-8"))
            write_json(materials / path.name, runtime_document(original))
        elif path.name in ("__init__.py", "manifest.json"):
            shutil.copy2(path, materials / path.name)
    for path in (root / "peritheos/data/studies").glob("*.json"):
        original = json.loads(path.read_bytes())
        datasets = original.get("datasets", [])
        if datasets and all("columns" in d for d in datasets):
            write_json(data / "studies" / path.name, runtime_study_document(original))
    # Structural CIFs and their scoped notices remain; no experimental tables.
    cifs = data / "datasets"
    cifs.mkdir()
    for name in (
        "brucite-zigan-1967-cod-1000054.cif",
        "vanadium-james-1960-9012770.cif",
        "superhydrous-phase-b-koch-muller-2005-mkm0105.cif",
        "mgsio3-post-perovskite-murakami-2004-cod-9009217.cif",
        "pattern-source-cod-9007338.cif",
        "pattern-source-cod-1566661.cif",
        "pattern-source-cod-9010029.cif",
        "pattern-source-cod-9006822.cif",
        "pattern-source-cod-9005757.cif",
        "pattern-source-cod-9001848.cif",
        "calculated-patterns.LICENSE.md",
        "cod-cif-sources.LICENSE.md",
        "dioptas-materials.LICENSE.md",
        "eos-structure-compilations.LICENSE.md",
        "cif-license-review.json",
    ):
        shutil.copy2(root / "peritheos/data/datasets" / name, cifs / name)
    project = (root / "pyproject.toml").read_text(encoding="utf-8")
    project = project.replace(
        'build-backend = "peritheos_build"\nbackend-path = ["."]',
        'build-backend = "maturin"',
    )
    project = re.sub(
        r"license-files = \[.*?\]\n",
        'license-files = ["LICENSE", "LICENSES/*.txt", "peritheos/data/LICENSING.md", "peritheos/data/datasets/*.LICENSE.md", "peritheos/data/third_party/rust/**/*"]\n',
        project,
        count=1,
        flags=re.S,
    )
    start = project.index("[tool.maturin]")
    end = project.index("[tool.ruff]", start)
    project = (
        project[:start]
        + """[tool.maturin]
manifest-path = "crates/peritheos-python/Cargo.toml"
module-name = "peritheos._rust"
python-source = "."
include = [{path = "peritheos/data/**/*", format = ["sdist", "wheel"]}, {path = "LICENSES/*", format = "sdist"}]
exclude = ["**/__pycache__/**", "**/*.pyc"]

"""
        + project[end:]
    )
    project = re.sub(r"\[tool.uv.sources\].*?(?=\n\[|\Z)", "", project, flags=re.S)
    (destination / "pyproject.toml").write_text(project, encoding="utf-8")
    return destination


def stage_datasets(destination, root=ROOT):
    """Export only the explicitly reviewed dataset declarations and resources."""
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError(f"Export destination already exists: {destination}")
    destination.mkdir(parents=True)
    package = destination / "peritheos_datasets"
    package.mkdir()
    (package / "__init__.py").write_text(
        f'"""Optional, source-scoped Peritheos observation resources."""\n__version__ = "{version(root)}"\n',
        encoding="utf-8",
    )
    policy = json.loads(
        (root / "packages/peritheos-datasets/datasets.json").read_text(encoding="utf-8")
    )
    review = None
    reviewed_entries = {}
    if "redistribution_review" in policy:
        receipt = policy["redistribution_review"]
        review_path = root / receipt["path"]
        payload = review_path.read_bytes()
        if hashlib.sha256(payload).hexdigest() != receipt["sha256"]:
            raise ValueError(
                "Numerical redistribution review changed; review the policy binding"
            )
        review = json.loads(payload)
        reviewed_entries = {
            e["review_key"]: e
            for e in review["entries"] + review.get("supporting_studies", [])
            if "review_key" in e
        }
        write_json(package / "redistribution-review.json", review)
    definitions = {}
    resource_hashes = {}
    for entry in policy["datasets"]:
        licensed = entry["source_status"] in {
            "open_license_declared",
            "open_source_declared",
            "author_permission_declared",
        }
        if entry.get("review_key"):
            assessment = reviewed_entries.get(entry.get("review_key"), {})
            if not (
                assessment.get("decision") == "admit"
                and assessment.get("basis") == entry.get("redistribution_basis")
                and assessment.get("declaration_sha256") == entry["declaration_sha256"]
                and assessment.get("material") == entry["material"]
                and assessment.get("dataset") == entry["dataset"]
                and assessment.get("source_status") == entry["source_status"]
            ):
                raise ValueError(
                    f"Optional dataset has an unapproved source scope: {entry['dataset']}"
                )
        elif not licensed:
            raise ValueError(
                f"Optional dataset has an unapproved source scope: {entry['dataset']}"
            )
        material = entry["material"]
        document = json.loads(
            (root / "peritheos/data/materials" / f"{material}.eosmat").read_text(
                encoding="utf-8"
            )
        )
        ds = next(
            d for d in document["datasets"] if d["identifier"] == entry["dataset"]
        )
        if digest(ds) != entry["declaration_sha256"]:
            raise ValueError(
                f"Optional dataset declaration changed; review {material}/{entry['dataset']}"
            )
        if entry.get("review_key"):
            if reviewed_entries[entry["review_key"]].get("resource") != ds.get(
                "resource"
            ):
                raise ValueError(
                    f"Numerical review resource binding changed: {entry['dataset']}"
                )
        definitions.setdefault(material, []).append(ds)
        if "resource" in ds:
            resource = ds["resource"]
            relative = resource["path"]
            source = root / "peritheos/data" / relative
            checksum = hashlib.sha256(source.read_bytes()).hexdigest()
            if checksum != resource["sha256"]:
                raise ValueError(f"Optional dataset resource changed: {relative}")
            target = package / "data" / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            resource_hashes[relative] = checksum
    for name, datasets in definitions.items():
        write_json(package / "declarations" / f"{name}.json", datasets)
    study_definitions = {}
    for entry in policy.get("study_datasets", []):
        # A study is its own scope; it must never borrow a material approval.
        assessment = reviewed_entries.get(entry.get("review_key"), {})
        if "material" in entry or not (
            assessment.get("decision") == "admit"
            and assessment.get("basis") == entry.get("redistribution_basis")
            and assessment.get("study") == entry["study"]
            and assessment.get("dataset") == entry["dataset"]
            and assessment.get("source_status") == entry["source_status"]
            and assessment.get("declaration_sha256") == entry["declaration_sha256"]
        ):
            raise ValueError(f"Unapproved supporting-study scope: {entry['dataset']}")
        documents = [
            json.loads(path.read_bytes())
            for path in (root / "peritheos/data/studies").glob("*.json")
        ]
        document = next(d for d in documents if d["identifier"] == entry["study"])
        ds = next(
            d for d in document["datasets"] if d["identifier"] == entry["dataset"]
        )
        if digest(ds) != entry["declaration_sha256"] or ds.get(
            "resource"
        ) != assessment.get("resource"):
            raise ValueError(
                f"Supporting-study declaration changed: {entry['dataset']}"
            )
        resource = ds["resource"]
        relative = resource["path"]
        source = root / "peritheos/data" / relative
        if hashlib.sha256(source.read_bytes()).hexdigest() != resource["sha256"]:
            raise ValueError(f"Supporting-study resource changed: {relative}")
        target = package / "data" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        resource_hashes[relative] = resource["sha256"]
        study_definitions.setdefault(entry["study"], []).append(ds)
    for name, datasets in study_definitions.items():
        write_json(package / "study-declarations" / f"{name}.json", datasets)
    registry = json.loads(
        (root / "peritheos/data/observation-views.json").read_text(encoding="utf-8")
    )
    admitted = {(e["material"], e["dataset"]) for e in policy["datasets"]}
    registry["datasets"] = [
        e
        for e in registry["datasets"]
        if (e["material_id"], e["dataset_id"]) in admitted
    ]
    write_json(package / "observation-views.json", registry)
    write_json(
        package / "manifest.json",
        {
            "version": version(root),
            "policy": policy,
            "resources": resource_hashes,
            "dataset_count": len(policy["datasets"]),
            "study_dataset_count": len(policy.get("study_datasets", [])),
            "total_dataset_count": len(policy["datasets"])
            + len(policy.get("study_datasets", [])),
        },
    )
    # Retain all scoped legal notices, never the neighboring source artifacts.
    notices = package / "notices"
    notices.mkdir()
    shutil.copy2(
        root / "peritheos/data/eos-structure-rights-review.json",
        package / "eos-structure-rights-review.json",
    )
    for source in (root / "peritheos/data/datasets").glob("*.LICENSE.md"):
        shutil.copy2(source, notices / source.name)
    shutil.copy2(
        root / "peritheos/data/datasets/ma_2024_sources/LICENSE.md",
        notices / "ma-2024.LICENSE.md",
    )
    if review is not None:
        shutil.copy2(
            root / "packages/peritheos-datasets/NUMERICAL-TRANSCRIPTIONS.LICENSE.md",
            notices / "numerical-transcriptions.LICENSE.md",
        )
    shutil.copytree(root / "LICENSES", destination / "LICENSES")
    shutil.copy2(root / "LICENSE", destination / "LICENSE")
    shutil.copy2(
        root / "packages/peritheos-datasets/README.md", destination / "README.md"
    )
    (destination / "pyproject.toml").write_text(
        f'''[build-system]
requires = ["setuptools>=77", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "peritheos-datasets"
version = "{version(root)}"
description = "Optional source-scoped observation datasets for Peritheos"
readme = "README.md"
requires-python = ">=3.9"
dependencies = ["peritheos=={version(root)}"]
license = "LicenseRef-Peritheos-Distribution"
license-files = ["LICENSE", "LICENSES/*.txt", "peritheos_datasets/notices/*"]

[tool.setuptools.packages.find]
include = ["peritheos_datasets*"]

[tool.setuptools.package-data]
peritheos_datasets = ["**/*"]
''',
        encoding="utf-8",
    )
    return destination


@contextlib.contextmanager
def _working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _build(hook, output, config_settings=None, metadata_directory=None, *, data=False):
    output = str(Path(output).resolve())
    if metadata_directory is not None:
        metadata_directory = str(Path(metadata_directory).resolve())
    with tempfile.TemporaryDirectory(prefix="peritheos-export-") as temporary:
        stage = Path(temporary) / "source"
        (stage_datasets if data else stage_core)(stage)
        if not data:
            # Reuse compilation outputs without copying the large target tree.
            os.environ.setdefault("CARGO_TARGET_DIR", str(ROOT / "target"))
        with _working_directory(stage):
            if data:
                from setuptools import build_meta as backend
            else:
                import maturin as backend
            fn = getattr(backend, hook)
            if hook == "build_sdist" or hook.startswith("prepare_metadata"):
                return fn(output, config_settings)
            return fn(output, config_settings, metadata_directory)


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    return _build("build_wheel", wheel_directory, config_settings, metadata_directory)


def build_sdist(sdist_directory, config_settings=None):
    return _build("build_sdist", sdist_directory, config_settings)


def prepare_metadata_for_build_wheel(metadata_directory, config_settings=None):
    return _build(
        "prepare_metadata_for_build_wheel", metadata_directory, config_settings
    )


def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    import maturin

    return maturin.build_editable(wheel_directory, config_settings, metadata_directory)


def get_requires_for_build_wheel(config_settings=None):
    return []


get_requires_for_build_sdist = get_requires_for_build_wheel
get_requires_for_build_editable = get_requires_for_build_wheel


def prepare_metadata_for_build_editable(metadata_directory, config_settings=None):
    import maturin

    return maturin.prepare_metadata_for_build_wheel(metadata_directory, config_settings)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=["core", "datasets"])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print((stage_core if args.kind == "core" else stage_datasets)(args.output))
