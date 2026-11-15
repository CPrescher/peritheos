"""Repository maintenance: expose the recovered garnet tables, without fitting.

The original audit resources remain immutable. This creates typed, attributed
numerical tables; original author attachments, articles and figures stay out of
the optional distribution. ``--check`` verifies catalog declarations/resources.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/data/garnet-angel-2022"
PREFIX = "datasets/garnet-angel-2022"
PYROPE_IDS = ["pyrope_angel_2022_author_pvt", "pyrope_angel_2022_author_pkt"]
ANGEL = {
    "authors": ["Angel", "Gilio", "Mazzucchelli", "Alvaro"],
    "year": 2022,
    "source": "Contributions to Mineralogy and Petrology",
    "volume": "177",
    "locator": "54",
    "doi": "10.1007/s00410-022-01918-5",
}


AUTHOR_REFERENCES = {
    "zhang_helium": (
        "Zhang et al.",
        1998,
        "Physics and Chemistry of Minerals",
        "10.1007/s002690050118",
    ),
    "zhang_neon": (
        "Zhang et al.",
        1998,
        "Physics and Chemistry of Minerals",
        "10.1007/s002690050118",
    ),
    "zou_pv": ("Zou et al.", 2012, "Journal of Applied Physics", "10.1063/1.4736407"),
    "zou_ks": ("Zou et al.", 2012, "Journal of Applied Physics", "10.1063/1.4736407"),
    "du_pv": (
        "Du, Clark and Walker",
        2015,
        "American Mineralogist",
        "10.2138/am-2015-4752",
    ),
    "milani_pv": ("Milani et al.", 2015, "Lithos", "10.1016/j.lithos.2015.03.017"),
    "milani_vt": ("Milani et al.", 2015, "Lithos", "10.1016/j.lithos.2015.03.017"),
    "hartwig_vt": (
        "Hartwig and Galkin",
        2021,
        "Journal of Thermal Analysis and Calorimetry",
        "10.1007/s10973-020-09396-2",
    ),
    "chantel_pv": (
        "Chantel et al.",
        2016,
        "American Mineralogist",
        "10.2138/am-2016-5335",
    ),
    "chantel_ks": (
        "Chantel et al.",
        2016,
        "American Mineralogist",
        "10.2138/am-2016-5335",
    ),
    "sinogeikin_2000_ks": (
        "Sinogeikin and Bass",
        2000,
        "Physics of the Earth and Planetary Interiors",
        "10.1016/S0031-9201(00)00143-6",
    ),
    "sinogeikin_2002_ks": (
        "Sinogeikin and Bass",
        2002,
        "Earth and Planetary Science Letters",
        "10.1016/S0012-821X(02)00851-8",
    ),
}


def column(name, quantity, unit, role="value", of=None, **extra):
    return dict(
        name=name,
        quantity=quantity,
        unit=unit,
        role=role,
        **({"of": of} if of else {}),
        **extra,
    )


def csv_bytes(names, rows):
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(names)
    writer.writerows(rows)
    return stream.getvalue().encode()


def declaration(
    identifier,
    kind,
    columns,
    payload,
    reference,
    location,
    description,
    *,
    used_by=(),
    **extra,
):
    return dict(
        identifier=identifier,
        kind=kind,
        description=description,
        reference=reference,
        source_location=location,
        columns=columns,
        resource=dict(
            path=f"{PREFIX}/{identifier}.csv",
            sha256=hashlib.sha256(payload).hexdigest(),
            media_type="text/csv",
        ),
        used_by_eos_records=list(used_by),
        source_license=None,
        contributor_rights="Contributor transcription, typed metadata and arrangement only; no third-party source license inferred.",
        **extra,
    )


def author_tables():
    path = AUDIT / "pyrope/author-2026-11-14-inputs.csv"
    rows = list(csv.DictReader(io.StringIO(path.read_text())))
    result = []
    for kind, identifier, quantity, value, unit in [
        ("V", PYROPE_IDS[0], "molar_volume", "molar_volume_cm3_mol", "cm^3/mol"),
        (
            "KS",
            PYROPE_IDS[1],
            "adiabatic_bulk_modulus",
            "bulk_modulus_adiabatic_gpa",
            "GPa",
        ),
    ]:
        cols = [
            column("temperature_k", "temperature", "K"),
            column(
                "temperature_assigned_error_k",
                "temperature",
                "K",
                "uncertainty",
                "temperature_k",
            ),
            column("pressure_gpa", "pressure", "GPa"),
            column(
                "pressure_assigned_error_gpa",
                "pressure",
                "GPa",
                "uncertainty",
                "pressure_gpa",
            ),
            column(value, quantity, unit),
            column(
                value.replace("_cm3_mol", "_assigned_error_cm3_mol")
                if kind == "V"
                else "bulk_modulus_assigned_error_gpa",
                quantity,
                unit,
                "uncertainty",
                value,
            ),
        ]
        trace = [
            "source_file",
            "source_line",
            "source_dataset",
            "source_row",
            "fit_row",
            "use",
            "group",
            "scale_index",
        ]
        cols += [
            column(
                n,
                n,
                "identifier" if n in {"source_file", "source_dataset"} else "1",
                "flag",
            )
            for n in trace
        ]
        output = [
            [
                *[r[f"{k}_token"] for k in ["T", "sT", "P", "sP", "Y", "sY"]],
                r["source_file"],
                r["source_line"],
                r["dataset"],
                r["row"],
                r["fit_row"],
                r["use"],
                r["group"],
                r["scale_index"],
            ]
            for r in rows
            if r["kind"] == kind
        ]
        payload = csv_bytes([c["name"] for c in cols], output)
        ds = declaration(
            identifier,
            "pressure_volume_temperature"
            if kind == "V"
            else "pressure_temperature_adiabatic_bulk_modulus",
            cols,
            payload,
            ANGEL,
            f"Author-processed {'Pyrope_PVT_used.dat' if kind == 'V' else 'Pyrope_PKT_used.dat'}, supplied by Ross Angel via Clemens Prescher on 2026-11-14; exact original line numbers retained.",
            f"All {len(output)} author-processed {'volume' if kind == 'V' else 'adiabatic-modulus'} inputs, with assigned coordinate errors and final-log row/group identities.",
            used_by=["pyrope_angel_2022_bm3_q_compromise_debye"],
            provenance_type="author_processed_fit_inputs",
            uncertainty={
                "reported_as": "author-assigned EosFit coordinate error widths",
                "covariance": "not_supplied",
                "qualification": "Preserved fit weights; not established raw-measurement standard deviations. The error columns deliberately use the generic uncertainty role.",
            },
            processing={
                "volume_basis": "author common molar-volume scale, cm3/mol",
                "volume_normalization": "Already applied by the author; upstream reductions and exclusions remain partial.",
                "scale_convention": "Group 1: 115 volumes, fixed scale 1. Group 2: 64 moduli, fixed scale 1, including two ambient Sinogeikin-Bass (2000) points. Group 3: 15 pressure-dependent Sinogeikin-Bass (2000) moduli; one refined multiplier acts on CALCULATED modulus, not these observed values.",
                "reference_temperature_k": 298.0,
                "source_archive": "Original attachments remain private; only attributed numeric transcriptions are distributed.",
            },
            upstream_sources="Source-dataset labels identify the individual studies; see the source audit for full citations.",
        )
        ds["upstream_sources"] = [
            {
                "source_dataset": label,
                "reference": dict(
                    authors=[AUTHOR_REFERENCES[label][0]],
                    year=AUTHOR_REFERENCES[label][1],
                    source=AUTHOR_REFERENCES[label][2],
                    doi=AUTHOR_REFERENCES[label][3],
                ),
            }
            for label in dict.fromkeys(r["dataset"] for r in rows if r["kind"] == kind)
        ]
        result.append((ds, payload, [str(path.relative_to(ROOT))]))
    return result


# Each source table is kept in its native coordinates, including anchors,
# excluded rows and extra quantities. No missing measurement error is filled.
SOURCES = [
    (
        "milani-2017-table-2.csv",
        "grossular_milani_2017_table2_vt",
        "temperature_lattice_volume",
        "Milani et al.",
        2017,
        "10.2138/am-2017-5855",
        "Official deposit, Table 2, sheet2 rows 5:39",
        "Published temperature, lattice parameter and cell volume, including source tokens and reported standard deviations.",
    ),
    (
        "milani-2017-table-3.csv",
        "grossular_milani_2017_table3_pv",
        "pressure_lattice_volume",
        "Milani et al.",
        2017,
        "10.2138/am-2017-5855",
        "Official deposit, Table 3, sheet3 rows 5:12",
        "All eight source rows, including the ambient normalization anchor; exact experimental temperature and original quartz readings are unavailable.",
    ),
    (
        "du-2015-final-table-2.csv",
        "grossular_du_2015_table2_vt",
        "temperature_lattice_volume",
        "Du, Clark and Walker",
        2015,
        "10.2138/am-2015-4752",
        "Final journal Table 2, p. 217, grossular block; checked against the preprint",
        "Final published table. Temperatures already include the source NaCl correction; the 0.2 K error bar is not identified as an independent Gaussian standard deviation. Calculated expansivity tokens are not measured observations.",
    ),
    (
        "bosenick-geiger-1997-table-2.csv",
        "grossular_bosenick_geiger_1997_table2_vt",
        "temperature_molar_volume",
        "Bosenick and Geiger",
        1997,
        "10.1029/97JB01612",
        "Table 2, author-uploaded primary text indexed by the source archive",
        "All 20 molar-volume rows including repeats. The original study excludes selected low-temperature cooling data; the exact Angel mask/anchor covariance is unresolved. Source PDF was not pixel-verified.",
    ),
    (
        "rodehorst-2002-table-4a.csv",
        "grossular_rodehorst_2002_table4a",
        "temperature_lattice_parameter",
        "Rodehorst et al.",
        2002,
        "10.2138/am-2002-0417",
        "Table 4A, p. 545, primary text indexed by RRUFF",
        "Ten lattice-parameter rows, including the 298 K anchor. Uniform 0.001 angstrom source error retained without assuming its confidence convention; temperature errors unavailable. Source PDF was not pixel-verified.",
    ),
    (
        "isaak-1992-table-2.csv",
        "grossular_isaak_1992_table2_gr2b_strain",
        "temperature_strain",
        "Isaak et al.",
        1992,
        "10.1007/BF00198608",
        "Table 2, p. 108; specimen discussion pp. 107-109",
        "Observed strain is GR-2B (76.1 mol% grossular, 22.1 mol% andradite), not the GGD specimen named in Angel S1. Calculated strain/endmember expansivity remain separately labelled source-calculated tokens.",
    ),
    (
        "isaak-1992-table-3b.csv",
        "grossular_isaak_1992_table3b_ggd_moduli",
        "temperature_adiabatic_bulk_modulus",
        "Isaak et al.",
        1992,
        "10.1007/BF00198608",
        "Table 3b, p. 112",
        "GGD (96.53 mol% grossular). Only three of 22 modulus errors are printed; the other 19 remain missing. No temperature errors are inferred.",
    ),
    (
        "kono-2010-table-1.csv",
        "grossular_kono_2010_table1",
        "pressure_temperature_volume_acoustics",
        "Kono et al.",
        2010,
        "10.1007/s12583-010-0112-2",
        "Table 1, p. 787",
        "All 36 rows and all tabulated quantities/uncertainties retained. Angel uses the moduli but excludes these XRD volumes. Derived density, velocities and moduli share inputs; they are not independent constraints.",
    ),
    (
        "conrad-1999-table-6.csv",
        "grossular_conrad_1999_table6_stiffness",
        "pressure_elastic_stiffness",
        "Conrad et al.",
        1999,
        "10.2138/am-1999-0321",
        "Table 6, p. 379",
        "All five stiffness rows, including the two excluded by Angel. Cij standard deviations are reported; C11/C12 covariance is unavailable. No KS error or 2% diagnostic allowance is inserted into the source table.",
    ),
    (
        "bass-1989-ambient.csv",
        "grossular_bass_1989_ambient_moduli",
        "ambient_elastic_moduli",
        "Bass",
        1989,
        "10.1029/JB094iB06p07621",
        "Original-author abstract in the Illinois institutional repository; cross-checked against Isaak Table 4",
        "99.0 mol% grossular sample. The recovered source gives no exact numerical temperature. No full primary PDF was recovered.",
    ),
    (
        "gwanmesia-2014-table-1.csv",
        "grossular_gwanmesia_2014_table1",
        "pressure_temperature_volume_acoustics",
        "Gwanmesia, Wang, Heady and Liebermann",
        2014,
        "10.1016/j.pepi.2013.09.010",
        "Final publisher HTML, Table 1 and associated methods/caption",
        "All 52 rows and 624 literal numerical cells. Temperature is Celsius (30 C = 303.15 K). Angel uses KS and excludes XRD volumes. Source relative error estimates are metadata, not invented per-row standard deviations.",
    ),
]


def source_column(name, *, rodehorst=False):
    # Raw tokens, contributor evidence and specimen/selection labels remain text.
    if name.endswith("_token") or name in {
        "specimen",
        "source",
        "evidence",
        "observed_specimen",
        "role",
        "selection_angel_2022",
    }:
        return column(
            name,
            name,
            "text",
            "flag",
            observation_kind="source_calculated"
            if name.startswith(("calculated_", "estimated_endmember_"))
            else "source_token_or_label",
        )
    if name in {"source_row", "worksheet_row", "repeat_measurement"}:
        return column(name, name, "1", "flag")
    if "temperature" in name:
        q, unit, value = (
            "temperature",
            "degC" if name.endswith("_c") else "K",
            "temperature_c" if name.endswith("_c") else "temperature_k",
        )
    elif "pressure" in name:
        q, unit, value = "pressure", "GPa", "pressure_gpa"
    elif name.startswith("molar_volume"):
        q, unit, value = "molar_volume", "cm^3/mol", "molar_volume_cm3_mol"
    elif name.startswith("cell_volume"):
        q, unit, value = "conventional_unit_cell_volume", "angstrom^3", "cell_volume_a3"
    elif name.startswith("lattice_a"):
        q, unit, value = "lattice_parameter_a", "angstrom", "lattice_a_angstrom"
    elif name.startswith("density"):
        q, unit, value = "density", "g/cm^3", "density_g_cm3"
    elif name.startswith(("vp_", "vs_")):
        q, unit, value = (
            ("compressional_velocity" if name.startswith("vp") else "shear_velocity"),
            "km/s",
            name.split("_", 1)[0] + "_km_s",
        )
    elif name.startswith("source_error_volume_strain"):
        q, unit, value = "volume_strain", "1", "observed_volume_strain"
    elif name in {"tp_microsecond", "ts_microsecond", "length_mm"}:
        q, unit, value = name, "mm" if name == "length_mm" else "microsecond", name
    else:
        q = "adiabatic_bulk_modulus" if name.startswith("ks_") else "elastic_modulus"
        unit, value = (
            "GPa",
            name.replace("_reported_error", "")
            .replace("_sd_gpa", "_gpa")
            .replace("_error", ""),
        )
    error = any(s in name for s in ["_sd", "_error", "_uncertainty"])
    role = (
        "standard_deviation"
        if ("_sd" in name and not rodehorst)
        else "uncertainty"
        if error
        else "value"
    )
    return column(name, q, unit, role, value if error else None)


def source_tables():
    result = []
    for (
        filename,
        identifier,
        kind,
        authors,
        year,
        doi,
        location,
        description,
    ) in SOURCES:
        source = AUDIT / "grossular" / filename
        payload = source.read_bytes()
        rows = list(csv.DictReader(io.StringIO(payload.decode())))
        cols = [source_column(n, rodehorst="rodehorst" in identifier) for n in rows[0]]
        if "table2_gr2b" in identifier:
            cols = (
                cols[:2]
                + [
                    column(
                        "observed_linear_strain",
                        "linear_strain",
                        "1",
                        derivation="source observed linear strain token times 1e-3",
                    ),
                    column(
                        "observed_volume_strain",
                        "volume_strain",
                        "1",
                        derivation="source observed volume strain token times 1e-2",
                    ),
                ]
                + cols[2:]
            )
            payload = csv_bytes(
                [c["name"] for c in cols],
                [
                    [
                        *list(r.values())[:2],
                        float(r["observed_linear_strain_x1e3_token"]) * 1e-3,
                        float(r["observed_volume_strain_x1e2_token"]) * 1e-2,
                        *list(r.values())[2:],
                    ]
                    for r in rows
                ],
            )
        ds = declaration(
            identifier,
            kind,
            cols,
            payload,
            dict(
                authors=[authors],
                year=year,
                source={
                    "10.1029/97JB01612": "Journal of Geophysical Research",
                    "10.1029/JB094iB06p07621": "Journal of Geophysical Research",
                    "10.1007/BF00198608": "Physics and Chemistry of Minerals",
                    "10.1007/s12583-010-0112-2": "Journal of Earth Science",
                    "10.1016/j.pepi.2013.09.010": "Physics of the Earth and Planetary Interiors",
                }.get(doi, "American Mineralogist"),
                doi=doi,
            ),
            location,
            description,
            provenance_type="source_table_transcription",
            related_eos_records=["grossular_angel_2022_bm3_q_compromise_debye"],
            usage_recommendation="Supporting upstream observations; not the recovered author-processed joint fit inputs. No missing uncertainty or normalization is inferred.",
            uncertainty={
                "covariance": "not_recovered",
                "missing_values": "unknown; blank cells load as NaN",
                "role": "Reported errors retain source definitions; assigned/estimated errors in diagnostic fits are excluded.",
            },
        )
        if "gwanmesia" in identifier:
            evidence = json.loads(
                (AUDIT / "grossular/gwanmesia-2014-source.json").read_bytes()
            )
            ds["uncertainty"].update(
                reported_relative_error_estimates_percent=evidence[
                    "reported_error_estimates_percent"
                ],
                qualification=evidence["error_qualification"],
            )
            ds["quantity_origins"] = {
                "cell_volume_a3": "XRD observation, excluded from Angel fit",
                "density_g_cm3": "source-derived from volume",
                "vp_km_s": "source-derived from length/transit time",
                "vs_km_s": "source-derived from length/transit time",
                "ks_gpa": "source-derived acoustic modulus; correlated inputs",
            }
        result.append((ds, payload, [str(source.relative_to(ROOT))]))
    # Numerical figure coordinates only: omit PDF-space geometry and artwork.
    sources = [
        AUDIT / "grossular" / f"figure4-thieblot-{suffix}.csv"
        for suffix in ["selected", "excluded-visible"]
    ]
    cols = [
        column("temperature_k", "temperature", "K"),
        column("molar_volume_cm3_mol", "molar_volume", "cm^3/mol"),
        column("selection", "selection", "identifier", "flag"),
        column("source_marker", "source_marker", "1", "flag"),
    ]
    rows = [
        [r["temperature_k"], r["molar_volume_cm3_mol"], selection, r["drawing_index"]]
        for source, selection in zip(sources, ["selected", "excluded_visible"])
        for r in csv.DictReader(io.StringIO(source.read_text()))
    ]
    payload = csv_bytes([c["name"] for c in cols], rows)
    identifier = "grossular_thieblot_1998_angel_2022_figure4"
    ds = declaration(
        identifier,
        "digitized_temperature_normalized_molar_volume",
        cols,
        payload,
        ANGEL,
        "Figure 4, PDF p. 8; Thiéblot series, seven selected and six excluded visible markers",
        "Digitized, already normalized plotted volumes attributed to Thieblot, not recovered primary measurement values or author fit inputs.",
        provenance_type="digitized_normalized_plot",
        related_eos_records=["grossular_angel_2022_bm3_q_compromise_debye"],
        digitization={
            "temperature_resolution_k": 2.0,
            "volume_resolution_cm3_mol": 0.003,
            "qualification": "Graphical resolution bounds, not measurement standard deviations or independent error bars.",
        },
        uncertainty={
            "measurement_errors": "not_recovered",
            "covariance": "not_recovered",
        },
    )
    result.append((ds, payload, [str(s.relative_to(ROOT)) for s in sources]))
    return result


def tables():
    return {"pyrope": author_tables(), "grossular": source_tables()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not args.check:
        parser.error(
            "Use --check; initial integration is a separate reviewed operation"
        )
    for material, entries in tables().items():
        document = json.loads(
            (ROOT / f"peritheos/data/materials/{material}.eosmat").read_bytes()
        )
        declarations = {d["identifier"]: d for d in document["datasets"]}
        for ds, payload, _ in entries:
            assert declarations[ds["identifier"]] == ds, ds["identifier"]
            assert (
                ROOT / "peritheos/data" / ds["resource"]["path"]
            ).read_bytes() == payload
    print(
        "Garnet declarations and numerical resources match the recovered audit sources."
    )


if __name__ == "__main__":
    main()
