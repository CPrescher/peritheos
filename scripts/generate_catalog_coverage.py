"""Generate catalog coverage from material documents and the public executable API.

Run as ``python -m scripts.generate_catalog_coverage``; ``--check`` is read-only.
Scientific audit and fit-reproduction statuses are reported independently.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from peritheos import list_eos_records, list_materials
from scripts.eosmat_io import read_material

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def collect_coverage() -> dict:
    """Join stored records to executable objects without hiding deferred records."""
    materials = {m.identifier: m for m in list_materials()}
    executable = {
        r.identifier: r for r in list_eos_records(include_reference_variants=True)
    }
    refits = json.loads(
        (DOCS / "data/primary-eos-refits.json").read_text(encoding="utf-8")
    )
    audits = json.loads(
        (ROOT / "peritheos/data/primary-source-audit.json").read_text(encoding="utf-8")
    )
    refit_ids = {r["record_identifier"] for r in refits["records"]}
    audit_ids = {r["record"] for r in audits["records"]}
    fingerprint = hashlib.sha256()
    rows = []
    records = []
    for path in sorted((ROOT / "peritheos/data/materials").glob("*.eosmat")):
        payload = path.read_bytes()
        fingerprint.update(path.name.encode() + b"\0" + payload + b"\0")
        document = read_material(path)
        identifier = document["identifier"]
        material = materials.get(identifier)
        kinds: Counter = Counter()
        validation: Counter = Counter()
        for raw in document["eos_records"]:
            record_id = raw["identifier"]
            live = executable.get(record_id)
            status = raw.get("scientific_validation", {}).get("status", "missing")
            validation[status] += 1
            if live is not None:
                kinds[live.equation_kind] += 1
            records.append(
                {
                    "identifier": record_id,
                    "material_identifier": identifier,
                    "executable": live is not None,
                    "equation_kind": live.equation_kind
                    if live
                    else raw.get("equation_kind"),
                    "scientific_validation": status,
                    "refit_ledger_present": record_id in refit_ids,
                    "source_audit_present": record_id in audit_ids,
                }
            )
        rows.append(
            {
                "identifier": identifier,
                "name": document["name"],
                "formula": document["formula"],
                "family_id": material.family_id
                if material
                else document.get("family_id"),
                "input_file": path.relative_to(ROOT).as_posix(),
                "stored_records": len(document["eos_records"]),
                "executable_records": sum(kinds.values()),
                "isothermal_records": kinds["isothermal"],
                "thermal_records": kinds["thermal"],
                "hugoniot_records": kinds["hugoniot"],
                "validation": dict(sorted(validation.items())),
                "datasets": len(document.get("datasets", [])),
            }
        )
    stored_ids = {r["identifier"] for r in records}
    if len(stored_ids) != len(records):
        raise ValueError("Duplicate record identifiers in material documents")
    if set(executable) - stored_ids:
        raise ValueError(
            "Executable API contains records absent from material documents"
        )
    calibration = Counter(
        (r.pressure_calibration or {}).get("status", "missing")
        for r in executable.values()
    )
    recalculation = Counter(
        (r.pressure_calibration or {}).get("recalculation", {}).get("status", "missing")
        for r in executable.values()
    )
    return {
        "format": "peritheos.catalog-coverage",
        "format_version": 1,
        "material_files_sha256": fingerprint.hexdigest(),
        "summary": {
            "stored_materials": len(rows),
            "stored_records": len(records),
            "executable_materials": len(materials),
            "executable_records": len(executable),
            "source_records": sum(
                not r.is_reference_variant for r in executable.values()
            ),
            "pressure_scale_reference_variants": sum(
                r.is_reference_variant for r in executable.values()
            ),
            "materials_with_thermal_eos": sum(r["thermal_records"] > 0 for r in rows),
            "executable_equation_kinds": dict(
                sorted(Counter(r.equation_kind for r in executable.values()).items())
            ),
            "stored_validation": dict(
                sorted(Counter(r["scientific_validation"] for r in records).items())
            ),
            "executable_calibration_status": dict(sorted(calibration.items())),
            "executable_recalculation_status": dict(sorted(recalculation.items())),
            "records_with_covariance": sum(
                r.parameter_covariance is not None for r in executable.values()
            ),
            "records_with_known_parameter_error_confidence": sum(
                r.parameter_error_confidence is not None for r in executable.values()
            ),
        },
        "ledgers": {
            "refit_outcomes": dict(
                sorted(Counter(r["status"] for r in refits["records"]).items())
            ),
            "missing_from_refit_ledger": sorted(stored_ids - refit_ids),
            "absent_from_catalog_refit_entries": sorted(refit_ids - stored_ids),
            "missing_from_source_audit": sorted(stored_ids - audit_ids),
            "absent_from_catalog_source_audit_entries": sorted(audit_ids - stored_ids),
        },
        "materials": rows,
        "records": records,
    }


def render_coverage(data: dict) -> str:
    """Render a compact overview plus a complete per-material coverage table."""
    s = data["summary"]
    lines = [
        "# Catalog coverage overview",
        "",
        "Generated from the current `.eosmat` documents and the public catalog API.",
        "This describes the checkout, including local material edits; it is not a release manifest.",
        "Regenerate with `python -m scripts.generate_catalog_coverage`; use `--check`",
        "to detect drift. The [machine-readable snapshot](data/catalog-coverage.json)",
        "includes every material and record, plus a fingerprint of the material files.",
        "",
        "Use the [active backlog](material-eos-candidates.md) for priorities and next actions.",
        "Scientific validation, fit reproduction, and availability answer different questions.",
        "A missing thermal EOS is a coverage gap, not a defect in a static record.",
        "",
        "## Current inventory",
        "",
        "| Measure | Count |",
        "|---|---:|",
    ]
    for label, key in [
        ("Stored material documents", "stored_materials"),
        ("Stored EOS records", "stored_records"),
        ("Executable materials", "executable_materials"),
        (
            "Executable EOS records, including references and qualified explicit-selection records",
            "executable_records",
        ),
        ("Source EOS choices in normal discovery", "source_records"),
        (
            "Published pressure-scale reference variants",
            "pressure_scale_reference_variants",
        ),
        ("Materials with an executable thermal EOS", "materials_with_thermal_eos"),
        ("Executable records with parameter covariance", "records_with_covariance"),
        (
            "Executable records with known parameter-error confidence",
            "records_with_known_parameter_error_confidence",
        ),
    ]:
        lines.append(f"| {label} | {s[key]} |")
    for key, title in [
        ("stored_validation", "Scientific validation of stored records"),
        ("executable_equation_kinds", "Executable equation kinds"),
        (
            "executable_calibration_status",
            "Pressure-calibration provenance of executable records",
        ),
        (
            "executable_recalculation_status",
            "Observation pressure-recalculation readiness",
        ),
    ]:
        lines.extend(["", f"## {title}", "", "| Status | Records |", "|---|---:|"])
        lines.extend(f"| `{status}` | {n} |" for status, n in s[key].items())
    lines.extend(
        [
            "",
            "Readiness values are the recorded assessments. A missing reference can mean a",
            "genuinely absent scale or a stale/unresolved link to an existing scale. `ready`",
            "does not establish original-fit parity or complete calibration-error propagation.",
            "P–T envelopes are marginal coverage, not measured rectangles or phase-stability fields.",
            "",
            "## Reproduction ledger and reconciliation",
            "",
            "Outcomes below describe the checked-in refit ledger, not an inferred assessment",
            "of records missing from that ledger. No scientific classifications are changed here.",
            "",
            "| Outcome | Records |",
            "|---|---:|",
        ]
    )
    lines.extend(
        f"| `{status}` | {n} |"
        for status, n in data["ledgers"]["refit_outcomes"].items()
    )
    for key, title in [
        ("missing_from_refit_ledger", "Catalog records absent from the refit ledger"),
        ("absent_from_catalog_refit_entries", "Refit entries absent from the catalog"),
        (
            "missing_from_source_audit",
            "Catalog records absent from the source-audit index",
        ),
        (
            "absent_from_catalog_source_audit_entries",
            "Source-audit entries absent from the catalog",
        ),
    ]:
        values = data["ledgers"][key]
        lines.extend(
            [
                "",
                f"**{title}:** "
                + (", ".join(f"`{v}`" for v in values) if values else "none."),
            ]
        )
    lines.extend(
        [
            "",
            "## Material coverage",
            "",
            "Counts of isothermal, thermal, and Hugoniot records include executable records only.",
            "Qualified means `not_reproduced`; deferred records remain stored source evidence.",
            "Dataset counts are declarations, not unique independent experiments or completeness claims.",
            "Multiple temperature slices do not imply a continuous thermal EOS.",
            "",
            "| Material | Formula | Stored | Executable | Isothermal | Thermal | Hugoniot | Qualified | Deferred | Datasets |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for r in data["materials"]:
        values = [
            r[k]
            for k in [
                "stored_records",
                "executable_records",
                "isothermal_records",
                "thermal_records",
                "hugoniot_records",
            ]
        ]
        values += [
            r["validation"].get("not_reproduced", 0),
            r["validation"].get("deferred", 0),
            r["datasets"],
        ]
        lines.append(
            f"| `{r['identifier']}` | {r['formula']} | "
            + " | ".join(map(str, values))
            + " |"
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = collect_coverage()
    outputs = {
        DOCS / "data/catalog-coverage.json": json.dumps(
            data, indent=2, ensure_ascii=False
        )
        + "\n",
        DOCS / "catalog-coverage.md": render_coverage(data),
    }
    for path, expected in outputs.items():
        if args.check:
            if not path.exists() or path.read_text(encoding="utf-8") != expected:
                raise SystemExit(
                    f"{path.relative_to(ROOT)} is stale; regenerate catalog coverage"
                )
        else:
            path.write_text(expected, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
