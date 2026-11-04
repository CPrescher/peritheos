"""Install only the audited original Sata (2010) static source record.

Append-only integration: existing records, defaults, structures and datasets
are preserved. Shared ledgers/reports are maintained separately.
"""

from __future__ import annotations

import copy
import csv
import json

from jsonschema import Draft202012Validator

from peritheos import convert_eosmat, eosmat_schema, validate_eosmat_document
from peritheos._eosmat_format import format_eosmat
from scripts.audit_sata_2010_fe3c import (
    DATASET_ID,
    NORMALIZED,
    PEAKS,
    QUALIFICATION,
    RECOMMENDATION_REASON,
    RECORD_ID,
    ROOT,
    SELECTED,
    SELECTED_ID,
    SOURCE,
    TOKENS,
    USAGE_RECOMMENDATION,
    VR,
    digest,
    ledger_outcome,
    observations,
    parse_token,
    reproduce,
    selected_rows,
)
from scripts.standardize_reconstruction_metadata import assessment_for_record

CARD = ROOT / "peritheos/data/materials/cementite_nonmagnetic.eosmat"
REFERENCE = {
    "authors": [
        "Sata, N.",
        "Hirose, K.",
        "Shen, G.",
        "Nakajima, Y.",
        "Ohishi, Y.",
        "Hirao, N.",
    ],
    "year": 2010,
    "title": "Compression of FeSi, Fe3C, Fe0.95O, and FeS under the core pressures and implication for light element in the Earth's core",
    "source": "Journal of Geophysical Research: Solid Earth",
    "volume": "115",
    "locator": "B09204",
    "doi": "10.1029/2009JB006975",
}


def source_record():
    finding = (
        "Final publisher Tables 1/6 and paragraph 15 verified against the original "
        "2009 manuscript for the common Fe3C tokens and fit triplet. Nine new "
        "states plus the three explicitly selected Li (2002) rows reproduce "
        "Pr/Kr/Kr-prime conditionally. Published coefficients preserved; no "
        "Litasov (2013) recalibrated coefficient or pressure is substituted."
        " Source verification and conditional coefficient parity coexist with "
        "poor agreement with some quoted coordinate-error estimates. Retained "
        "as a nondefault published EOS with explicit error/scatter qualifications; "
        "source confidence and covariance remain unspecified."
    )
    legacy = {
        "identifier": RECORD_ID,
        "label": "Sata et al. (2010), nonmagnetic Fe3C, original-scale 300 K BM3",
        "record_kind": "published",
        "determination_method": "experimental",
        "default": False,
        "temperature_ref": 300,
        "eos": {
            "type": "BM3",
            "model": "birch_murnaghan_3",
            "parameters": {"V0": VR, "K0": 290.0, "K0_prime": 3.76},
        },
        "parameter_errors": {"V0": None, "K0": 13.0, "K0_prime": 0.18},
        "parameter_error_confidence": None,
        "fixed_parameters": [],
        "experimental_pressure_range_gpa": [27.2, 186.6],
        "experimental_temperature_range_k": [300, 300],
        "pressure_range_status": "reported_exactly",
        "validity": {
            "pressure_gpa": [27.2, 186.6],
            "temperature_k": [300, 300],
            "notes": [
                "Source-assigned nonmagnetic branch. Combined coverage includes three Li (2002) points above 25 GPa; Sata's nine new states span 49.87–186.60 GPa.",
                "V0 is an extrapolated branch parameter, not a measured stable ambient nonmagnetic cell. No thermal extrapolation or magnetic-phase diagram.",
                "Minor Fe7C3 and Re peaks occur at the highest pressure; Fe3C peaks were refined separately. Existing Prescher Pnma structure is a distinct experiment.",
            ],
        },
        "range_provenance": {
            "pressure": "Final Table 1 nine rows plus paragraph 15's Li measured rows >25 GPa: 27.2, 30.3 and 30.5 GPa.",
            "temperature": "Final paragraph 6: XRD after thermal annealing at room temperature; Table 6/figures use 300 K. Annealing temperatures are not PVT observations.",
        },
        "parameter_provenance": {
            "determination_method": "Experimental room-temperature XRD compression of Fe3C with measured MgO cells plus separately calibrated Li (2002) measurements.",
            "equation": "Final Eq. 1 mBM3 at selected Vr; at the printed Pr=0, it equals Eq. 2 BM3. Final Table 6 prints the zero-pressure BM3 representation. These are one curve, not separate source models.",
            "V0": "Final Table 6 gives 9.341 A3/atom. Z=4 and four atoms/formula give 16 atoms/cell: 149.456 A3/cell. The paragraph's 149.46 is rounded. No V0 error or original ordinary-BM3 fixed-V0 constraint is reported; fixed_parameters is not inferred from the conditional diagnostic.",
            "K0": "Final paragraph 15: Kr=290 +/-13 GPa at Vr=9.341 A3/atom and fitted Pr=0.0 +/-1.6 GPa; Table 6 retains 290 GPa.",
            "K0_prime": "Final paragraph 15: Kr-prime=3.76 +/-0.18; final Table 6.",
            "normalization": "Table 1 Fe3C has Z=4. The whole conventional P cell contains Fe12C4. Source (a,b,c)=(standard c, standard a, standard b); volume is invariant. Existing structure/sites are unchanged.",
        },
        "published_reference_fit": {
            "model": "mBM3",
            "reference_volume_atom_a3": 9.341,
            "reference_volume_cell_a3": VR,
            "Pr_gpa": 0.0,
            "Pr_error_gpa": 1.6,
            "Kr_gpa": 290.0,
            "Kr_error_gpa": 13.0,
            "Kprime_r": 3.76,
            "Kprime_r_error": 0.18,
            "confidence": None,
            "covariance": None,
            "note": "Vr is selected; Pr is a fitted output, not a source-fixed zero. The catalog executes the publication's Table 6 BM3 representation at printed Pr=0.",
        },
        "parameter_uncertainty_convention": "Parenthetical source widths retained; confidence/covariance unspecified. V0 has no reported error. Reference-fit Pr error is retained separately, not silently assigned as a V0 error. Diagnostic V0 fixing is not a source constraint.",
        "fit_datasets": [SELECTED_ID],
        "reference": REFERENCE,
        "pressure_calibration": {
            "status": "partially_resolved",
            "audit_date": "2026-11-01",
            "methods": [
                {
                    "kind": "equation_of_state",
                    "reference_eos_record": "mgo_speziale_2001_bm3_2",
                    "reference": {
                        "authors": ["Speziale", "Zha", "Duffy", "Hemley", "Mao"],
                        "year": 2001,
                        "source": "Journal of Geophysical Research: Solid Earth",
                        "doi": "10.1029/2000JB900318",
                    },
                    "scope": "All nine new Sata Fe3C Table 1 rows only",
                    "source_location": "Final paragraphs 5–6, Table 1 footnote b (prints 2000), paragraph 26. Paired observed MgO a and V retained.",
                    "notes": "Source BM3 74.71 A3/cell,160.2 GPa,3.99 replay differs by up to 0.353 GPa. Exact pressure reduction not recovered; printed pressures are unchanged.",
                },
                {
                    "kind": "ruby_fluorescence",
                    "source_location": "Sata (2010), paragraph 15; Li (2002), Table 1 and footnote b, original pressure calibration described in the Li source audit.",
                    "reference": {
                        "authors": ["Li", "Mao", "Fei", "Gregoryanz", "Eremets", "Zha"],
                        "year": 2002,
                        "source": "Physics and Chemistry of Minerals",
                        "doi": "10.1007/s00269-001-0224-4",
                    },
                    "scope": "Three Li (2002) rows >25 GPa; ruby/neon calibration and reported pressure widths remain the original Li values.",
                    "notes": "No MgO or Litasov recalibration is applied. The existing Li audit retains its paired neon cells and missing raw ruby wavelengths.",
                },
            ],
            "recalculation": {
                "status": "not_possible",
                "notes": "Partial MgO equation replay is available on all nine paired cells. A complete exact author-pressure replay also needs the actual reducer/precision and Li ruby input/combination settings. No observation-level source pressure has been replaced.",
            },
        },
        "scientific_validation": {
            "status": "primary_source_validated",
            "report_path": "docs/literature-reproductions/sata-2010-fe3c.md",
            "audit_date": "2026-11-01",
            "note": finding,
            "data_agreement_status": "poor",
            "usage_recommendation": USAGE_RECOMMENDATION,
            "recommendation_date": "2026-11-01",
            "recommendation_audit": "docs/literature-reproductions/sata-2010-fe3c.md#conditional-fit-replay-and-residuals",
            "recommendation_reason": RECOMMENDATION_REASON,
            "warning": "Some quoted point errors are smaller than the overall experimental scatter. Retain the published EOS within 300 K and 27.2–186.6 GPa; the separate refit is diagnostic.",
            "primary_source_check": {
                "access_url": "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009JB006975",
                "doi": REFERENCE["doi"],
                "locations": [
                    "Tables 1 and 6",
                    "Equations 1–2",
                    "paragraphs 5–6, 12 and 14–15",
                    "Table 4 diffraction spacings",
                    "paragraphs 24–27, pressure-scale consistency",
                ],
                "finding": finding,
            },
            "verified_fields": [
                "original_Fe3C_rows_and_paired_MgO",
                "source_BM3_and_mBM3_equivalence_at_printed_Pr_zero",
                "atomic_to_cell_volume",
                "published_parameters_and_errors",
                "source_joint_selection",
                "conditional_parameter_parity",
                "source_pressure_scale_provenance",
                "temperature_and_pressure_limits",
            ],
            "primary_data_check": {
                "status": "bundled",
                "dataset_identifiers": [DATASET_ID, SELECTED_ID],
                "finding": "All nine original Fe3C/MgO rows and the unchanged three explicitly cited Li rows >25 GPa retained. Table 4 measured/calculated spacings are provenance only. Conditional source-width parameter agreement coexists with poor data agreement: six of nine Sata rows exceed first-order combined P/V coordinate widths. Published EOS retained with explicit uncertainty/data-scatter qualifications; exact weights/calibration workflow remain unresolved.",
            },
            "numerical_reproduction": {
                "reproduction_report": "docs/data/sata-2010-fe3c-audit.json",
                "author_fit_reproduced": False,
                "mBM3_joint_equal_pressure_coefficients": reproduce()[
                    "modified_bm3_replays"
                ]["equal_pressure"]["parameters_Pr_Kr_Kprime"],
                "notes": "Source-equation and conditional parameter checks, not independent physical/calibration validation or recovery of original weights.",
                "data_agreement": reproduce()["data_agreement"],
            },
            "unresolved": [
                "Original objective, weights, confidence, covariance and unrounded fit output unavailable.",
                "Printed MgO reduction leaves up to 0.353 GPa discrepancy; exact reducer and complete Li pressure replay unavailable.",
                "Final Table 4 heading says 193 GPa, versus 187 GPa in text/Figure 3 and 186.60 GPa for measurement 193; preserved as a source discrepancy.",
            ],
        },
        "notes": QUALIFICATION,
    }
    container = {
        "format_version": 3,
        "format": "peritheos.material",
        "identifier": "cementite_nonmagnetic",
        "name": "Cementite, nonmagnetic",
        "formula": "Fe3C",
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        },
        "datasets": datasets(),
        "eos_records": [legacy],
    }
    native = convert_eosmat(container, version=4, inline=True)["eos_records"][0]
    row = {
        "record_identifier": RECORD_ID,
        "material_identifier": "cementite_nonmagnetic",
        **ledger_outcome(legacy),
    }
    native["evidence"]["reproduction"] = {
        "assessment": assessment_for_record(native, row)
    }
    return native


def columns_for(rows):
    columns = []
    for name in rows[0]:
        quantity, unit = (
            ("temperature", "K")
            if name == "temperature_k"
            else ("pressure", "GPa")
            if "pressure" in name
            else ("conventional_unit_cell_volume", "angstrom^3")
            if "volume" in name
            else ("lattice_parameter", "angstrom")
            if "angstrom" in name
            else ("ratio", "1")
            if "vc_vfe" in name
            else ("source_order", "1")
        )
        col = {
            "name": name,
            "quantity": quantity,
            "unit": unit,
            "role": "uncertainty"
            if name.startswith("error_")
            else "flag"
            if name in {"source", "source_row", "run", "measurement_id"}
            else "value",
        }
        if name.startswith("error_"):
            col["of"] = name.removeprefix("error_")
        columns.append(col)
    return columns


def datasets():
    definitions = []
    for identifier, rows, path, description in [
        (
            DATASET_ID,
            observations()[1],
            NORMALIZED,
            "All nine Fe3C Table 1 observations, all lattice/volume/pressure error widths and paired MgO readings; source a,b,c setting retained.",
        ),
        (
            SELECTED_ID,
            selected_rows(),
            SELECTED,
            "Source-named joint selection: nine Sata MgO-referenced rows and three Li (2002) ruby/neon rows above 25 GPa; no pressure re-reduction.",
        ),
    ]:
        definitions.append(
            {
                "identifier": identifier,
                "kind": "pressure_volume_temperature",
                "description": description,
                "source_location": "Sata (2010) final Table 1 and paragraph 15; Li (2002) Table 1 rows T1-12/T1-13/T1-14 for source=1 in the selected input.",
                "reference": REFERENCE,
                "row_count": len(rows),
                "columns": columns_for(rows),
                "resource": {
                    "path": "datasets/" + path.name,
                    "sha256": digest(path),
                    "media_type": "text/csv",
                },
                "used_by_eos_records": [RECORD_ID],
                "uncertainty_convention": "Source last-digit widths; statistical confidence/covariance unspecified. No missing error is zero. Temperature is fixed at the source 300 K isotherm, without invented per-row T error.",
                "provenance": {
                    "source_tokens": "datasets/" + TOKENS.name,
                    "source_tokens_sha256": digest(TOKENS),
                    "manifest": "datasets/" + SOURCE.name,
                    "source_flags": {
                        "0": "Sata Table 1, MgO scale",
                        "1": "Li Table 1, original ruby/neon scale; same observations as cementite_li_2002_table1",
                    },
                    "missing_run_ids": "Li has no published run/measurement labels; run and measurement_id remain empty/null. source_row retains the original T1 row.",
                },
                "license": "Contributor transcription/normalization rights only CC0-1.0; third-party source numerical license unspecified. See cementite-sata-2010.LICENSE.md.",
            }
        )
    with PEAKS.open(newline="", encoding="utf-8") as stream:
        raw = list(csv.DictReader(stream))
    values = []
    for row in raw:
        d, error = parse_token(row["d_observed_raw"])
        values.append(
            [
                *[int(row[n]) for n in ("h", "k", "l")],
                d,
                error,
                float(row["d_calculated_raw"]),
                float(row["delta_d_raw"]),
            ]
        )
    definitions.append(
        {
            "identifier": "cementite_sata_2010_table4_diffraction",
            "kind": "diffraction_structure",
            "description": "All 12 observed and calculated Table 4 d-spacings; indexed Fe3C phase check, not EOS fit inputs.",
            "reference": REFERENCE,
            "source_location": "Final Table 4. Its heading says 193 GPa; text/Figure 3 say 187 GPa. No pressure correction inferred.",
            "columns": [
                {"name": n, "quantity": "miller_index", "unit": "1", "role": "value"}
                for n in ("h", "k", "l")
            ]
            + [
                {
                    "name": "d_observed_angstrom",
                    "quantity": "d_spacing",
                    "unit": "angstrom",
                    "role": "value",
                },
                {
                    "name": "error_d_angstrom",
                    "quantity": "d_spacing",
                    "unit": "angstrom",
                    "role": "uncertainty",
                    "of": "d_observed_angstrom",
                },
                {
                    "name": "d_calculated_source_angstrom",
                    "quantity": "d_spacing",
                    "unit": "angstrom",
                    "role": "value",
                },
                {
                    "name": "delta_d_source_angstrom",
                    "quantity": "d_spacing",
                    "unit": "angstrom",
                    "role": "value",
                },
            ],
            "rows": values,
            "row_count": 12,
            "used_by_eos_records": [],
            "provenance_only": True,
            "notes": "Measured and source-calculated d are separate; retain rounded source delta-d literally. Source axes/hkl are not relabelled as the existing Prescher Pnma structure. No observed diffraction intensity is supplied.",
            "source_resource": "datasets/" + PEAKS.name,
            "source_sha256": digest(PEAKS),
        }
    )
    return definitions


def integrate():
    original = json.loads(CARD.read_text())
    document = copy.deepcopy(original)
    if any(r["identifier"] == RECORD_ID for r in document["eos_records"]):
        raise ValueError(
            "Sata record already present; review any replacement separately"
        )
    document["eos_records"].append(source_record())
    document["datasets"].extend(datasets())
    validate_eosmat_document(document)
    Draft202012Validator(eosmat_schema(4)).validate(document)
    rendered = format_eosmat(document)
    restored = json.loads(rendered)
    assert restored["eos_records"][:-1] == original["eos_records"]
    assert restored["datasets"][: len(original["datasets"])] == original["datasets"]
    assert all(
        restored[key] == value
        for key, value in original.items()
        if key not in {"eos_records", "datasets"}
    )
    CARD.write_text(rendered, encoding="utf-8")
    print(
        "Appended one original-scale static Sata record; existing records, default and structure preserved"
    )


if __name__ == "__main__":
    integrate()
