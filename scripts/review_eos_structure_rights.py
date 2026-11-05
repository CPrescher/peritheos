"""Repository maintenance: bind reviewed EOS/structure rights to unchanged sources.

This inventories exact reviewed content. It does not infer source permissions,
alter scientific records, select datasets or clear new records automatically.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path

from peritheos_build import digest, runtime_document

ROOT = Path(__file__).resolve().parents[1]
POLICY = "docs/data/eos-structure-rights-policy.json"
REVIEW = "peritheos/data/eos-structure-rights-review.json"
INDEX = "docs/data/license-inventory/eos-structure-review.csv"
GROUPS = "docs/data/license-inventory/eos-source-groups.csv"
STRUCTURE_FIELDS = (
    "formula",
    "phase",
    "units",
    "formula_units_per_cell",
    "symmetry",
    "space_group",
    "space_group_number",
    "lattice",
    "atom_sites",
    "peaks",
)
STRUCTURE_PROVENANCE = (
    "reference",
    "source",
    "notes",
    "structure_reference",
    "structure_provenance",
    "phase_references",
)
RATIONALES = {
    "independently_calculated_pattern": (
        "Independent kinematic X-ray reflection array, not an inherited card or "
        "measured intensities. Input models, source checksums, proxy/occupancy "
        "limits and omitted corrections are explicit in calculated-pattern-inputs.json. "
        "Six CIFs carry COD CC0 declarations; one ideal model is independently "
        "authored. FeH2 and two SnO2 models use attributed numerical structural "
        "parameters under a scoped factual-content assessment, without a separate "
        "publication grant. CC0 covers contributor-held generated-output rights "
        "only; EOS/observation/source-publication rights remain separate."
    ),
    "numerical_parameter_assessment": (
        "Normalized numerical EOS coefficients, uncertainty/convention metadata and "
        "independently implemented equations, with source attribution. No source "
        "table layout, author software or source artwork is licensed by this "
        "assessment. Related dataset grants are not extended to these parameters. "
        "Source/collection/contract scope remains qualified in the notice."
    ),
    "thermodynamic_collection_review_required": (
        "Seven four-coefficient isothermal extracts share the official ds62 "
        "thermodynamic collection. The recovered artifact/README contain no "
        "explicit redistribution grant; academic availability is not general "
        "permission. Assess collection/acquisition rights or obtain an exact grant "
        "before claiming unrestricted clearance. No THERMOCALC software is shipped."
    ),
    "numerical_structure_assessment": (
        "Normalized lattice, symmetry, coordinate/occupancy or explicitly idealized "
        "structural facts. Preserve the cited source, proxy/normalization and "
        "missing-structure qualifications. COD references without a verified resource "
        "binding remain provenance leads, not new grants. No blanket publication, "
        "database or acquisition-contract clearance is asserted."
    ),
    "cod_source_and_numerical_changes": (
        "Named checksum-bound COD resource has a verified CC0 source "
        "declaration. The material's stated normalization and, for superhydrous "
        "phase B, separately attributed erratum correction remain explicit. "
        "This does not extend COD CC0 to other publications or to the EOS."
    ),
    "no_structural_values": (
        "No lattice, occupied atom-site array or peak list is supplied. Chemical "
        "identity/phase labels and source context remain; do not invent or infer a "
        "missing structural model. This is not clearance of all material metadata."
    ),
    "legacy_peak_collection_review_required": (
        "The exact imported hkl/d/intensity array matches the pinned MIT Dioptas "
        "source, whose notice is retained. Earlier beamline-card acquisition and "
        "collection ownership are unverified. A .jcpds extension is not proof of "
        "ICDD origin or an open license. Trace/recreate independently or obtain "
        "a grant for any protected collection before claiming clearance."
    ),
    "icsd_acquisition_review_required": (
        "Existing provenance cites ICSD 158959 and the Murakami publication. No "
        "authenticated acquisition path or freely licensed replacement binding "
        "was established. Published numerical facts and ICSD contract/database "
        "rights are separate. Obtain primary/open-source lineage or permission; "
        "do not claim that scaling the cell clears the copied topology."
    ),
}


def structural_source(document):
    return {k: document[k] for k in STRUCTURE_PROVENANCE if k in document}


def csv_bytes(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def build_review(root=ROOT):
    policy_bytes = (root / POLICY).read_bytes()
    policy = json.loads(policy_bytes)
    paths = sorted((root / "peritheos/data/materials").glob("*.eosmat"))
    actual = {p.relative_to(root).as_posix() for p in paths}
    if actual != set(policy["material_sha256"]):
        raise ValueError("Material set changed; a new rights review is required")
    entries, groups = [], defaultdict(list)
    migrations = []
    for path in paths:
        relative = path.relative_to(root).as_posix()
        content = path.read_bytes()
        if hashlib.sha256(content).hexdigest() != policy["material_sha256"][relative]:
            raise ValueError(f"Reviewed material changed: {relative}")
        document = json.loads(content)
        runtime = runtime_document(document)
        material = document["identifier"]
        for index, record in enumerate(document["eos_records"]):
            reference = record["reference"]
            source_key = reference.get("doi") or digest(reference)
            decision = (
                "thermodynamic_collection_review_required"
                if source_key == policy["thermodynamic_collection_doi"]
                else "numerical_parameter_assessment"
            )
            migration = record["evidence"]["scientific_validation"].get(
                "migration_source"
            )
            if migration:
                if migration["commit"] != policy["dioptas"]["commit"]:
                    raise ValueError("Unreviewed migration source")
                migrations.append(f"{material}/{record['identifier']}")
            entry = {
                "key": f"eos/{material}/{record['identifier']}",
                "kind": "eos",
                "material": material,
                "identifier": record["identifier"],
                "evidence": f"{relative}#/eos_records/{index}",
                "source_sha256": digest(record),
                "runtime_sha256": digest(runtime["eos_records"][index]),
                "assessment": decision,
                "source_license": "No new source grant established by this review",
                "reference": reference,
                "parameter_provenance_pointer": f"{relative}#/eos_records/{index}/evidence/parameter_provenance",
                "parameter_provenance_sha256": digest(
                    record["evidence"]["parameter_provenance"]
                ),
                "migration_source": migration,
            }
            entries.append(entry)
            groups[source_key].append(entry)
        fields = {k: document[k] for k in STRUCTURE_FIELDS if k in document}
        if material in policy["calculated_patterns"]:
            decision = "independently_calculated_pattern"
            if (
                digest(document["peaks"])
                != policy["calculated_patterns"][material]["peaks_sha256"]
            ):
                raise ValueError("Calculated pattern changed; review required")
        elif material in policy["legacy_patterns"]:
            decision = "legacy_peak_collection_review_required"
            if (
                digest(document["peaks"])
                != policy["legacy_patterns"][material]["peaks_sha256"]
            ):
                raise ValueError("Legacy peak source comparison changed")
        elif material in policy["icsd_lineage_materials"]:
            decision = "icsd_acquisition_review_required"
        elif material in policy["cod_structures"]:
            decision = "cod_source_and_numerical_changes"
        elif not any(document.get(k) for k in ("lattice", "atom_sites", "peaks")):
            decision = "no_structural_values"
        else:
            decision = "numerical_structure_assessment"
        entries.append(
            {
                "key": f"structure/{material}",
                "kind": "structure_and_diffraction",
                "material": material,
                "identifier": material,
                "evidence": relative + "#/",
                "source_sha256": digest(fields),
                "runtime_sha256": digest(
                    {k: runtime[k] for k in STRUCTURE_FIELDS if k in runtime}
                ),
                "assessment": decision,
                "source_license": policy["calculated_patterns"][material][
                    "source_terms"
                ]
                if material in policy["calculated_patterns"]
                else "Scoped COD CC0 declaration; corrections separate"
                if material in policy["cod_structures"]
                else "No new source grant established by this review",
                "site_count": len(document.get("atom_sites", [])),
                "peak_count": len(document.get("peaks", [])),
                "provenance_sha256": digest(structural_source(document)),
                "provenance_fields": [k for k in STRUCTURE_PROVENANCE if k in document],
                "cod_id": policy["cod_structures"].get(material),
                "legacy_source_comparison": policy["legacy_patterns"].get(material),
            }
        )
    counts = dict(sorted(Counter(e["assessment"] for e in entries).items()))
    review = {
        "schema": "peritheos.eos-structure-rights-review.v1",
        "reviewed_on": policy["reviewed_on"],
        "baseline_commit": policy["baseline_commit"],
        "policy_sha256": hashlib.sha256(policy_bytes).hexdigest(),
        "catalog_scientific_sha256": policy["catalog_scientific_sha256"],
        "scope": policy["scope"],
        "limits": policy["limits"],
        "counts": {
            "materials": len(paths),
            "eos_records": sum(e["kind"] == "eos" for e in entries),
            "structural_components": sum(
                e["kind"] == "structure_and_diffraction" for e in entries
            ),
            "eos_source_groups": len(groups),
            "migration_records": len(migrations),
            "assessments": counts,
        },
        "notice": "datasets/eos-structure-compilations.LICENSE.md",
        "dioptas": policy["dioptas"],
        "thermodynamic_source": policy["thermodynamic_source"],
        "primary_sources": policy["primary_sources"],
        "post_perovskite_cod_replacement": policy["post_perovskite_cod_replacement"],
        "calculated_pattern_replacement": policy["calculated_pattern_replacement"],
        "material_sha256": policy["material_sha256"],
        "rationales": RATIONALES,
        "entries": entries,
    }
    index = [
        {
            "key": e["key"],
            "kind": e["kind"],
            "material": e["material"],
            "identifier": e["identifier"],
            "assessment": e["assessment"],
            "source_terms": e["source_license"],
            "evidence": e["evidence"],
            "source_sha256": e["source_sha256"],
            "runtime_sha256": e["runtime_sha256"],
            "required_action": RATIONALES[e["assessment"]],
        }
        for e in entries
    ]
    sources = [
        {
            "source_key": key,
            "reference": json.dumps(values[0]["reference"], ensure_ascii=False),
            "entries": len(values),
            "materials": "|".join(sorted({v["material"] for v in values})),
            "assessment": values[0]["assessment"],
            "scope": "Exact attributed normalized parameter sets, not an entire source collection or a row-count exemption; repeated fits/calibration variants are not independent experiments.",
            "component_keys": "|".join(v["key"] for v in values),
        }
        for key, values in sorted(groups.items())
    ]
    return {
        REVIEW: (json.dumps(review, indent=2, ensure_ascii=False) + "\n").encode(),
        INDEX: csv_bytes(index),
        GROUPS: csv_bytes(sources),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    for name, content in build_review().items():
        path = ROOT / name
        if args.check:
            if not path.is_file() or path.read_bytes() != content:
                raise SystemExit(f"Stale EOS/structure review: {name}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
    print("EOS/structure rights review: exact source bindings and outputs match")


if __name__ == "__main__":
    main()
