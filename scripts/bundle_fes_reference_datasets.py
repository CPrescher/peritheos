"""Package recovered FeS author inputs and literature cold data separately.

No fit coefficients or source observations are changed. The original data/EOS
files are retained byte-for-byte; CSVs expose them through the Dataset API.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_sata_2010_fes_vi import observations as sata_observations

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
CARD = ROOT / "peritheos/data/materials/fes_vi.eosmat"
AUTHOR = ROOT / "docs/data/fes-author-eosfit/sources"
SOURCE = DATA / "morard_2026_author_sources"
RECORD = "fes_vi_morard_2026_bm3_mgd"
PREFIX = "fes_morard_2026_author_"
REFERENCE = json.loads(CARD.read_text())["eos_records"][0]["reference"]
VALUE_COLUMNS = [
    ("pressure_gpa", "pressure", "GPa"),
    ("pressure_uncertainty_gpa", "pressure", "GPa"),
    ("temperature_k", "temperature", "K"),
    ("temperature_uncertainty_k", "temperature", "K"),
    ("molar_volume_cm3_mol", "molar_volume", "cm^3/mol"),
    ("molar_volume_uncertainty_cm3_mol", "molar_volume", "cm^3/mol"),
]
AUTHOR_NAMES = [
    "author_input_row",
    "observation_id",
    "source_group",
    "matched_sata2010_row",
    *[x[0] for x in VALUE_COLUMNS],
]
FOLLOWUP_NOTE = (
    "Fit reproduction is classified similar: the reproduced-cold staged coefficients are compatible with reported parameter uncertainties. "
    "The supplied Table S1 and recovered author inputs are bundled separately. "
    "The selected author input contains 146 thermal plus 21 literature cold rows. "
    "The 11 unsuitable quenched rows are excluded on the basis of user-reported author personal communication and are not database datasets. FeS6.eos confirms full "
    "BM3-MGD with q=1 and Tref=298 K, but differs from the printed parameters "
    "and was saved in March 2024. Actual EosFit replay gives thermal RMS "
    "1.632 GPa, maximum 5.058 GPa and 14 rows beyond +/-3 GPa. Conditional "
    "cold refits closely recover the saved cold coefficients. Candidate-input "
    "refits converge. The user confirms the supplied files are the complete author-used inputs; the selected rows are confirmed by author communication. "
    "Exact GUI regression controls and agreement with printed coefficients and residuals remain unresolved. Published coefficients stay unchanged and "
    "deferred. Selected-data, cold-then-thermal EosFit fits at Tr=300 K give "
    "gamma0=2.29252(3829) for fixed published cold coefficients and "
    "2.46418(3953) for fixed weighted-reproduced cold coefficients; both "
    "exceed the +/-3 GPa residual bound. Both reproduced-cold staged fits have cold coefficients within published error widths; equal-weight gamma0 is also within the published width, and weighted gamma0 error intervals overlap the published interval. Parameter agreement is compatible within reported uncertainties; exact source-output and residual-envelope reproduction are separate. See "
    "docs/literature-reproductions/fes-eosfit-direct.md."
)


def sha(payload):
    return hashlib.sha256(payload).hexdigest()


def encode_csv(names, rows):
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(names)
    writer.writerows(rows)
    return buffer.getvalue().encode()


def match_additional_cold_gsas(author):
    """Identify the 300 K author points by lattice volume and KCl pressure."""
    path = DATA / "morard_2026_sources/FeS-GSAS.csv"
    with path.open(newline="", encoding="utf-8-sig") as handle:
        raw = list(csv.reader(handle, delimiter=";"))
    blocks = []
    current = None
    for line, row in enumerate(raw, 1):
        if row and row[0].startswith("#"):
            current = {"scan": row[0], "source_line": line}
            blocks.append(current)
        elif current is not None and row and row[0] == "FeSVI":
            current["lattice_angstrom"] = [float(x) for x in row[3:6]]
            current["temperature_k"] = float(row[13])
        elif current is not None and row and row[0] == "KCl B2":
            current["marker_a_angstrom"] = float(row[3])
    matches = []
    for author_row in author:
        volume = float(author_row[8])
        candidates = [
            b
            for b in blocks
            if b.get("temperature_k") == 300
            and abs(math.prod(b["lattice_angstrom"]) * 0.15055 - volume) < 5.1e-9
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"Additional cold observation has {len(candidates)} GSAS matches"
            )
        b = candidates[0]
        x = b["marker_a_angstrom"] / 54.5 ** (1 / 3)
        pressure = (
            3 * 17.2 * (1 - x) / x**2 * math.exp(1.5 * (5.89 - 1) * (1 - x)) - 0.0028
        )
        difference = pressure - float(author_row[4])
        if abs(difference) > 5e-8:
            raise ValueError(
                "Additional cold pressure does not match paired KCl replay"
            )
        matches.append(
            {
                "observation_id": author_row[1],
                "author_input_row": author_row[0],
                **b,
                "volume_angstrom3_from_lattice": math.prod(b["lattice_angstrom"]),
                "replayed_pressure_gpa": pressure,
                "pressure_difference_gpa": difference,
            }
        )
    return {
        "source_resource": "datasets/morard_2026_sources/FeS-GSAS.csv",
        "source_sha256": sha(path.read_bytes()),
        "matches": matches,
        "marker_temperature_k": 298.75,
        "marker_temperature_law": "0.75*T_FeS + 0.25*295 K; numerically reconstructed, not explicitly printed in the paper",
        "marker_eos": "kcl_b2_dewaele_2012_vinet_3",
        "qualification": "All eleven points match the public GSAS export deposited with Morard's study, including 300 K lattice refinements and paired KCl pressures. This identifies the deposited data, not necessarily its original authorship. The paper states that Figure 6 cold data and the cold fit use Sata/Ohfuji literature observations; it does not individually identify these eleven additional input rows. Original-study attribution and measurement protocol remain unconfirmed.",
    }


def author_rows(filename):
    lines = (AUTHOR / filename).read_text().splitlines()
    if lines[0].split() != ["Format", "1", "P", "sigP", "T", "sigT", "V", "sigV"]:
        raise ValueError("Unexpected original column order")
    counts = Counter()
    sata = [x for x in sata_observations() if x["included"]]
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        tokens = line.split()
        if len(tokens) != 6:
            raise ValueError("Incomplete original author row")
        p, dp, t, dt, v, dv = map(float, tokens)
        group = (
            "morard_thermal"
            if t > 300
            else "additional_cold"
            if dp == 0.0056
            else "literature_cold"
        )
        counts[group] += 1
        matches = (
            [
                s["table_fes_row"]
                for s in sata
                if abs(s["volume_angstrom3"] * 0.15055 - v) < 1e-8
            ]
            if group == "literature_cold"
            else []
        )
        if len(matches) > 1:
            raise ValueError("Ambiguous literature volume match")
        rows.append(
            [
                len(rows) + 1,
                f"{group}_{counts[group]:03d}",
                group,
                matches[0] if matches else "",
                *tokens,
            ]
        )
    return rows


def metadata(identifier, filename, rows, label, source_filename, notes):
    columns = [
        {
            "name": name,
            "quantity": "source_order"
            if name in ("author_input_row", "matched_sata2010_row")
            else "source_identifier",
            "unit": "1",
            "role": "flag",
        }
        for name in AUTHOR_NAMES[:4]
    ]
    for i, (name, quantity, unit) in enumerate(VALUE_COLUMNS):
        column = {
            "name": name,
            "quantity": quantity,
            "unit": unit,
            "role": "uncertainty" if i % 2 else "value",
        }
        if i % 2:
            column["of"] = VALUE_COLUMNS[i - 1][0]
        columns.append(column)
    payload = encode_csv(AUTHOR_NAMES, rows)
    dataset = {
        "identifier": identifier,
        "kind": "observations",
        "description": label,
        "reference": REFERENCE,
        "source_location": f"Author-supplied {source_filename}; original data rows (header excluded), retrieved 2026-10-03.",
        "source_url": "https://doi.org/10.1103/h4pj-rvxx",
        "license": "unspecified",
        "used_by_eos_records": [],
        "related_eos_records": [RECORD],
        "columns": columns,
        "resource": {
            "path": "datasets/" + filename,
            "media_type": "text/csv",
            "sha256": sha(payload),
            "row_count": len(rows),
        },
        "uncertainty": {
            "type": "unspecified",
            "confidence_level": None,
            "notes": "Author-assigned EosFit uncertainty columns, preserved without an established one-sigma confidence. Literature cold rows use 1% P, 0.5% V and 5 K. Thermal errors match Table S1 apart from rounded volume conversion; dV uses linear lattice-error addition.",
        },
        "provenance": {
            "source_resource": "datasets/morard_2026_author_sources/" + source_filename,
            "source_sha256": sha((AUTHOR / source_filename).read_bytes()),
            "manifest": "datasets/morard_2026_author_sources/manifest.json",
            "observation_id_semantics": "Shared across the selected input and literature subset. Overlapping observations and the duplicate thermal rows are not independent additional measurements.",
            "final_author_fit_selection_confirmed": True,
            "selection_provenance": "User-reported personal communication from Guillaume Morard confirms exclusion of the 11 unsuitable quenched observations; recorded 2026-10-03.",
        },
        "notes": notes
        + " Associated with the Morard source audit; the 167-row selection is confirmed by user-reported author personal communication. The author files do not label each cold row's original study. matched_sata2010_row indicates a cell-volume match, not an independently verified original-study assignment.",
    }
    return dataset, payload


def build():
    short_name = "FeS6EOSfittot-SataOhfuji-300K.dat"
    short = author_rows(short_name)
    datasets, outputs = [], {}
    for suffix, filename, rows, label, original, notes in [
        (
            "without_additional_cold_pvt",
            "fes-morard-2026-author-without-additional-cold.csv",
            short,
            "Morard author EosFit input: 167 rows without additional cold data",
            short_name,
            "The -300K filename removes the 11 additional 300 K rows; this is not a cold-only dataset. All 146 thermal and 21 literature cold rows remain.",
        ),
        (
            "literature_cold_pv",
            "fes-morard-2026-author-literature-cold.csv",
            [x for x in short if x[2] == "literature_cold"],
            "Sata/Ohfuji combined cold reference data as used in Morard's author input",
            short_name,
            "21 literature cold observations at 300 K. Thirteen match the Sata (2010) re-reported VI cell volumes; eight further source assignments are unresolved. One matched P is 101.2 GPa here versus 101.1 GPa in the Sata table. Assigned author fit errors differ from the published table errors.",
        ),
    ]:
        dataset, payload = metadata(
            PREFIX + suffix, filename, rows, label, original, notes
        )
        datasets.append(dataset)
        outputs[DATA / filename] = payload
    sata = [x for x in sata_observations() if x["included"]]
    sata_names = [
        "source_row",
        "run",
        "pressure_gpa",
        "pressure_uncertainty_gpa",
        "temperature_k",
        "volume_angstrom3",
        "volume_uncertainty_angstrom3",
        "pressure_raw",
        "volume_raw",
    ]
    sata_rows = [
        [
            x["table_fes_row"],
            x["run"],
            x["pressure_gpa"],
            x["error_pressure_gpa"],
            300,
            x["volume_angstrom3"],
            x["error_volume_angstrom3"],
            x["pressure_raw"],
            x["volume_raw"],
        ]
        for x in sata
    ]
    payload = encode_csv(sata_names, sata_rows)
    filename = "fes-sata-2010-table1-vi.csv"
    outputs[DATA / filename] = payload
    columns = [
        {
            "name": n,
            "quantity": q,
            "unit": u,
            "role": r,
            **({"of": target} if target else {}),
        }
        for n, q, u, r, target in [
            ("source_row", "source_order", "1", "flag", None),
            ("run", "source_identifier", "1", "flag", None),
            ("pressure_gpa", "pressure", "GPa", "value", None),
            (
                "pressure_uncertainty_gpa",
                "pressure",
                "GPa",
                "uncertainty",
                "pressure_gpa",
            ),
            ("temperature_k", "temperature", "K", "value", None),
            ("volume_angstrom3", "volume", "angstrom^3", "value", None),
            (
                "volume_uncertainty_angstrom3",
                "volume",
                "angstrom^3",
                "uncertainty",
                "volume_angstrom3",
            ),
            ("pressure_raw", "source_token", "1", "flag", None),
            ("volume_raw", "source_token", "1", "flag", None),
        ]
    ]
    datasets.append(
        {
            "identifier": "fes_sata_2010_table1_vi_pv",
            "kind": "observations",
            "description": "Sata (2010) Table 1 re-reported FeS-VI cold observations with original table errors",
            "reference": {
                "authors": ["Sata"],
                "year": 2010,
                "doi": "10.1029/2009JB006975",
                "source": "Journal of Geophysical Research",
                "locator": "B09204",
            },
            "source_location": "Sata (2010), Table 1 FeS-VI rows, factual transcription from publisher HTML; original parenthesized errors and run IDs retained.",
            "source_url": "https://doi.org/10.1029/2009JB006975",
            "license": "unspecified",
            "used_by_eos_records": [],
            "related_eos_records": [RECORD],
            "columns": columns,
            "resource": {
                "path": "datasets/" + filename,
                "media_type": "text/csv",
                "sha256": sha(payload),
                "row_count": len(sata),
            },
            "uncertainty": {
                "type": "unspecified",
                "confidence_level": None,
                "notes": "Errors in last digits decoded from original tokens; confidence and cross-row covariance unknown. Temperature is the 300 K isotherm with no quoted uncertainty, not a measured zero-error temperature.",
            },
            "provenance": {
                "source_resource": "datasets/morard_2026_author_sources/fes-sata-2010-table1.csv",
                "manifest": "datasets/morard_2026_author_sources/manifest.json",
                "selection": "All 13 phase-labelled VI rows; six VII rows excluded. No rows excluded by pressure.",
            },
            "notes": "Separate from the 21-row Morard literature input and its assigned relative errors. Table 1 re-reports earlier Ohfuji (2007)/Sata (2008) data; complete original-study identity and per-row pressure-marker data remain unresolved. Z=4 cell volumes preserved. Coexisting VI measurements at high pressure retained.",
        }
    )
    source_entries = []
    for original in [
        AUTHOR / short_name,
        AUTHOR / "FeS6.eos",
        ROOT / "docs/data/fes-sata-2010-table1.csv",
    ]:
        payload = original.read_bytes()
        outputs[SOURCE / original.name] = payload
        source_entries.append(
            {
                "filename": original.name,
                "sha256": sha(payload),
                "bytes": len(payload),
                "role": "published numerical transcription"
                if original.name.endswith(".csv")
                else "user-supplied author file",
                "license": "unspecified",
                "bundled": True,
            }
        )
    outputs[SOURCE / "manifest.json"] = (
        json.dumps(
            {
                "format": "peritheos.primary-source-manifest",
                "audit_date": "2026-10-03",
                "sources": source_entries,
                "datasets": [x["identifier"] for x in datasets],
                "qualification": "Original source versions and their assigned uncertainty conventions are kept separate. Input variants and subsets overlap; do not concatenate them as independent observations. Saved EOS is March 2024, not confirmed final publication model.",
            },
            indent=2,
        )
        + "\n"
    ).encode()
    return datasets, outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    datasets, outputs = build()
    document = json.loads(CARD.read_text())
    old = [
        x
        for x in document["datasets"]
        if not x["identifier"].startswith(PREFIX)
        and x["identifier"] != "fes_sata_2010_table1_vi_pv"
    ]
    document["datasets"] = old + datasets
    record = document["eos_records"][0]
    validation = record["scientific_validation"]
    validation["audit_date"] = "2026-10-03"
    validation["reproduction_status"] = "similar"
    validation["residual_claim_reproduction_status"] = "not_reproduced"
    validation["note"] = FOLLOWUP_NOTE
    validation["unresolved"] = [
        "printed publication coefficients versus the supplied author-used saved model",
        "published +/-3 GPa residual bound versus 5.058 GPa with supplied author files",
        "GUI versus console regression controls, weighting and selected groups",
        "individual provenance and pressure calibration of all 21 literature cold rows",
    ]
    validation["primary_source_check"]["finding"] = FOLLOWUP_NOTE
    for location in [
        "Author-supplied FeS6.eos",
        "Both FeS6EOSfittot-SataOhfuji input variants",
        "Supplementary_Materials_Rev_140526.docx, Figures S1-S5",
    ]:
        if location not in validation["primary_source_check"]["locations"]:
            validation["primary_source_check"]["locations"].append(location)
    check = validation["primary_data_check"]
    check["audit_date"] = "2026-10-03"
    check["dataset_identifiers"] = [x["identifier"] for x in document["datasets"]]
    check["finding"] = (
        "Table S1, the selected 167-row author input, its 21 literature cold rows, and the 13-row Sata VI transcription are bundled separately with original uncertainty semantics. The 11 unsuitable quenched rows and the 178-row variant containing them are excluded from database datasets, based on user-reported author personal communication."
    )
    record["parameter_provenance"]["debye_temperature_law"] = (
        "Full integrated MGD confirmed for the supplied March 2024 FeS6.eos (Thermal=7, Param14=0, q=1). The printed source coefficients and 300 K reference are retained; the supplied author-used file has different V0, gamma0 and Tref. The user confirms these are the complete author-used inputs; numerical publication parity remains unresolved."
    )
    record["author_input_followup"] = {
        "manifest": "datasets/morard_2026_author_sources/manifest.json",
        "audit_report": "docs/data/fes-author-eosfit/manifest.json",
        "source_model_resource": "datasets/morard_2026_author_sources/FeS6.eos",
        "saved_author_parameters": {
            "V0_cm3_mol": 15.37,
            "K0_gpa": 115.49,
            "K0_prime": 4.99,
            "gamma0": 2.41671,
            "Tr_k": 298,
            "theta0_k": 417,
            "q": 1,
            "n": 2,
            "q_compromise": False,
        },
        "dataset_identifiers": [x["identifier"] for x in datasets],
        "publication_fit_mask_confirmed": True,
        "author_used_inputs_confirmed": True,
        "author_used_inputs_confirmation": "User confirms the supplied data and EOS files are all inputs used by Guillaume Morard; recorded 2026-10-03. These are treated as the author-used inputs, not missing-data candidates.",
        "publication_selection_provenance": "User-reported personal communication from Guillaume Morard confirms exclusion of the 11 unsuitable quenched observations; recorded 2026-10-03.",
        "staged_audit_report": "docs/data/fes-author-staged/manifest.json",
        "parameter_comparison_report": "docs/data/fes-author-staged/parameter-comparison.json",
        "parameter_reproduction_assessment": "Similar within reported uncertainties for reproduced-cold staged fits: all cold coefficients lie within published error widths; equal-weight gamma0 lies within the published width and weighted gamma0 intervals overlap. Printed-cold-fixed gamma0 refit differs beyond overlapping error intervals. The residual-envelope discrepancy is a separate source-output claim, not evidence that the data cannot be fitted.",
        "cold_row_attribution": "Figure 6 and Section III.D identify the paper's cold data/fit as Sata (2008) and Ohfuji (2007) literature data. The eleven additional 300 K points match the GSAS deposit but are not individually attributed in the paper; original-study authorship remains unconfirmed.",
        "reference_temperature_qualification": "The paper shows a 300 K isotherm and uses ambient-temperature literature data; it does not explicitly state the EosFit Tref setting. The printed-coefficient source record uses 300 K as its declared audit interpretation, whereas the supplied saved file explicitly sets Tref=298 K.",
    }
    outputs[CARD] = (json.dumps(document, indent=1, ensure_ascii=False) + "\n").encode()
    for path, payload in outputs.items():
        if args.check:
            if not path.exists() or path.read_bytes() != payload:
                raise SystemExit(f"Stale FeS dataset artifact: {path}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
    print(
        f"{'Checked' if args.check else 'Packaged'} {len(datasets)} FeS source datasets"
    )


if __name__ == "__main__":
    main()
