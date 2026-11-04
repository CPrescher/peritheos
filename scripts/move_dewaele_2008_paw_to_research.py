"""Move Table VI research calculations out of the main EOS catalog, losslessly.

Repository maintenance only. Normal catalog discovery uses data/materials and
never scans this collection. Published coefficients and scientific findings stay
unchanged; only report destinations move with the records.
"""

from __future__ import annotations

import copy
import hashlib
import json
import shutil
from collections import Counter

from peritheos._eosmat_format import format_eosmat
from peritheos.eosmat import validate_eosmat_document
from scripts.reproduce_dewaele_2008_metals import DATA, ROOT, Z
from scripts.validate_primary_eos_refits import render_markdown

COLLECTION = ROOT / "peritheos/data/research/dewaele-2008-paw"
RESEARCH_REPORT = ROOT / "docs/research/dewaele-2008-paw.md"


def is_paw(identifier):
    return "_dewaele_2008_paw_" in identifier


def relocate_record_links(record):
    reproduction = record["evidence"]["reproduction"]
    assessment = reproduction["assessment"]
    if "previous_catalog_report" not in reproduction:
        reproduction["previous_catalog_report"] = {
            k: assessment[k]
            for k in ("report", "report_path", "report_url")
            if k in assessment
        }
    anchor = "record-" + record["identifier"]
    assessment.update(
        report="reproduction.json",
        report_path="docs/research/dewaele-2008-paw.md#" + anchor,
        report_url="https://peritheos.readthedocs.io/en/latest/research/dewaele-2008-paw/#"
        + anchor,
    )
    if "source_equation_check" in reproduction:
        reproduction["source_equation_check"]["report"] = "reproduction.json"


def write_json(path, document):
    path.write_text(
        json.dumps(document, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )


def render_research_report():
    lines = [
        "# Dewaele 2008 PAW research EOSMAT collection",
        "",
        "The 54 Table VI PAW parameterizations are preserved as a separate research",
        "collection. They are excluded from the main EOS database, catalog searches,",
        "default selection and catalog ledgers. The 13 Table IV experimental EOS",
        "remain in the main catalog. Moving files does not change scientific findings.",
        "",
        "The collection contains nine native EOSMAT files, their Figure 4 checkpoint",
        "resources, a source-audit ledger and the complete numerical checks.",
        "Fifty sets have graphical source-equation parity; two bcc-Fe LDA sets lack",
        "an independent benchmark and two hcp-Fe LDA sets retain the unresolved",
        "table/semicore-curve identity discrepancy. Graph checks are not PAW refits.",
        "",
        "[Primary publication](https://doi.org/10.1103/PhysRevB.78.104102) ·",
        "[Source investigation](../literature-reproductions/dewaele-2008-metals.md) ·",
        "[Collection files](https://github.com/CPrescher/peritheos/tree/main/peritheos/data/research/dewaele-2008-paw)",
        "",
        "The full collection is shipped as package data and can be opened explicitly.",
        "These are conventional-cell EOSMAT volumes; the table below also shows the",
        "atomic normalization used in the publication. Parameter errors remain",
        "unreported. Static-lattice and corrected 298 K fits are distinct isotherms,",
        "not a continuous thermal model. The 298 K corrections include zero-point motion.",
        "",
        "```python",
        "from importlib.resources import as_file, files",
        "from peritheos import Dataset, Material, load_eosmat",
        "",
        'collection = files("peritheos.data").joinpath("research/dewaele-2008-paw")',
        'with as_file(collection.joinpath("silver.eosmat")) as path:',
        "    document = load_eosmat(path)",
        "    material = Material.from_eosmat(document)",
        'checkpoints = Dataset.from_mapping(document["datasets"][0], resource_root=collection)',
        "print(len(material.eos_records), len(checkpoints))",
        "```",
        "",
        "For an external copy, keep the nine files and their `datasets/` directory",
        "together and pass that collection directory as `resource_root`.",
        "Deferred and diagnostic Fe restrictions remain in force when opening files.",
        "",
    ]
    for material in Z:
        document = json.loads((COLLECTION / f"{material}.eosmat").read_text())
        z = document["formula_units_per_cell"]
        lines += [
            f"## {document['name']}",
            "",
            f"[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/{material}.eosmat)",
            "",
        ]
        for record in document["eos_records"]:
            assessment = record["evidence"]["reproduction"]["assessment"]
            q = record["parameters"]
            lines += [
                f'<a id="record-{record["identifier"]}"></a>',
                "",
                f"### {record['label']}",
                "",
                f"`{record['identifier']}`",
                "",
                "| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |",
                "| ---: | ---: | ---: | ---: | ---: |",
                f"| {q['V0']['value'] / z:.6g} | {q['V0']['value']:.6g} | {q['K0']['value']:.6g} | {q['K0_prime']['value']:.6g} | {record['temperature_ref']:g} |",
                "",
                f"**Outcome:** {assessment['status']}; **scope:** {assessment['scope']}.",
                "",
                assessment["basis"],
                "",
                *[f"- {x}" for x in assessment["qualifications"]],
                "",
            ]
    RESEARCH_REPORT.write_text("\n".join(lines).rstrip() + "\n")


def move():
    COLLECTION.mkdir(parents=True, exist_ok=True)
    (COLLECTION / "datasets").mkdir(exist_ok=True)
    moved = []
    for material in Z:
        path = ROOT / f"peritheos/data/materials/{material}.eosmat"
        document = json.loads(path.read_text())
        records = [r for r in document["eos_records"] if is_paw(r["identifier"])]
        if not records:
            continue
        research_path = COLLECTION / f"{material}.eosmat"
        if research_path.exists():
            raise ValueError(
                f"Research file already exists while main records remain: {research_path}"
            )
        research = copy.deepcopy(document)
        research["eos_records"] = copy.deepcopy(records)
        ids = {r["identifier"] for r in records}
        research["datasets"] = [
            copy.deepcopy(d)
            for d in document.get("datasets", [])
            if ids.intersection(d.get("used_by_eos_records", []))
        ]
        for dataset in research["datasets"]:
            dataset["used_by_eos_records"] = [
                x for x in dataset["used_by_eos_records"] if x in ids
            ]
            resource = dataset["resource"]["path"]
            destination = COLLECTION / resource
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / "peritheos/data" / resource, destination)
        for record in research["eos_records"]:
            relocate_record_links(record)
        validate_eosmat_document(research)
        research_path.write_text(format_eosmat(research))
        document["eos_records"] = [
            r for r in document["eos_records"] if r["identifier"] not in ids
        ]
        remaining = []
        for dataset in document.get("datasets", []):
            users = dataset.get("used_by_eos_records", [])
            if ids.intersection(users):
                dataset = copy.deepcopy(dataset)
                dataset["used_by_eos_records"] = [x for x in users if x not in ids]
                if not dataset["used_by_eos_records"]:
                    continue
            remaining.append(dataset)
        document["datasets"] = remaining
        validate_eosmat_document(document)
        path.write_text(format_eosmat(document))
        moved.extend(sorted(ids))
    archive_ids = {
        r["identifier"]
        for p in COLLECTION.glob("*.eosmat")
        for r in json.loads(p.read_text())["eos_records"]
    }
    assert len(archive_ids) == 54
    for source, target, key in [
        (
            ROOT / "docs/data/primary-eos-refits.json",
            COLLECTION / "reproduction.json",
            "record_identifier",
        ),
        (
            ROOT / "peritheos/data/primary-source-audit.json",
            COLLECTION / "source-audit.json",
            "record",
        ),
    ]:
        ledger = json.loads(source.read_text())
        selected = [r for r in ledger["records"] if r[key] in archive_ids]
        if selected:
            assert len(selected) == 54
            archive = {
                **ledger,
                "format": "peritheos.research-eos-audit",
                "catalog_scope": "research_only",
                "records": selected,
                "summary": dict(Counter(r["status"] for r in selected)),
            }
            write_json(target, archive)
            ledger["records"] = [
                r for r in ledger["records"] if r[key] not in archive_ids
            ]
            counts = dict(
                sorted(Counter(r["status"] for r in ledger["records"]).items())
            )
            total_key = "total" if key == "record_identifier" else "records"
            ledger["summary"] = {total_key: len(ledger["records"]), **counts}
            write_json(source, ledger)
            if key == "record_identifier":
                (ROOT / "docs/primary-eos-refits.md").write_text(
                    render_markdown(ledger)
                )
        assert target.exists()
    shutil.copy2(
        DATA / "dewaele-2008-source.json", COLLECTION / "publication-source.json"
    )
    manifest = {
        "format": "peritheos.research-eos-collection",
        "format_version": 1,
        "identifier": "dewaele_2008_paw",
        "doi": "10.1103/PhysRevB.78.104102",
        "catalog_scope": "research_only",
        "main_catalog_records": 13,
        "eos_records": 54,
        "material_files": [f"{m}.eosmat" for m in Z],
        "reproduction_report": "reproduction.json",
        "source_audit": "source-audit.json",
        "source_pdf_sha256": json.loads(
            (DATA / "dewaele-2008-source.json").read_text()
        )["source_pdf_sha256"],
        "files": {
            p.relative_to(COLLECTION).as_posix(): hashlib.sha256(
                p.read_bytes()
            ).hexdigest()
            for p in sorted(COLLECTION.rglob("*"))
            if p.is_file() and p.name not in {"manifest.json", "README.md"}
        },
    }
    write_json(COLLECTION / "manifest.json", manifest)
    (COLLECTION / "README.md").write_text(
        "# Dewaele 2008 PAW research collection\n\n"
        "54 published Table VI PAW parameterizations in nine native EOSMAT files.\n"
        "Research only: excluded from the main material database and catalog APIs.\n"
        "The 13 experimental Table IV EOS remain in the main catalog.\n\n"
        "All original coefficients, errors, source and validation findings are retained.\n"
        "Figure 4 graphical checkpoints are in datasets/; complete checks are in\n"
        "reproduction.json, primary-source findings in source-audit.json, and the\n"
        "publication transcription in publication-source.json. Fifty sets have\n"
        "source-equation parity; two bcc-Fe LDA sets are deferred and two hcp-Fe\n"
        "LDA sets retain their unresolved semicore-identity finding.\n\n"
        "Open an EOSMAT file explicitly with load_eosmat and Material.from_eosmat.\n"
        "Load its dataset metadata with Dataset.from_mapping(metadata,\n"
        "resource_root=collection_directory). The companion documentation is\n"
        "docs/research/dewaele-2008-paw.md in the repository.\n\n"
        "Primary source: https://doi.org/10.1103/PhysRevB.78.104102\n"
        "Source numerical facts retain attribution; the publisher PDF is not\n"
        "redistributed and no licence over third-party content is asserted.\n"
    )
    render_research_report()
    print(
        f"Moved {len(moved)} records; research collection has {len(archive_ids)} records."
    )


if __name__ == "__main__":
    move()
