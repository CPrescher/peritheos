#!/usr/bin/env python3
"""Import primary iron papers from explicitly selected local Zotero attachments.

Run with the PDF-capable Python runtime. Only factual tables are redistributed.
Only the enumerated import batch is updated; unrelated records are preserved.
Subsequent validation uses bundled CSVs.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
MATERIALS = ROOT / "peritheos/data/materials"
NA = 6.02214076e23
DATE = "2026-09-10"
SOURCES = {
    "brown2000": (
        "G37FEUSQ",
        "XZ8K7CI5",
        "10.1063/1.1319320",
        ["Brown", "Fritz", "Hixson"],
        2000,
        "Hugoniot data for iron",
        "Journal of Applied Physics 88, 5496–5498",
    ),
    "dubrovinsky2000": (
        "64IHSEQS",
        "TLKYMW2B",
        "10.1103/PhysRevLett.84.1720",
        ["Dubrovinsky", "Saxena", "Tutti", "Rekhi", "LeBehan"],
        2000,
        "In Situ X-Ray Study of Thermal Expansion and Phase Transition of Iron at Multimegabar Pressure",
        "Physical Review Letters 84, 1720–1723",
    ),
    "yamazaki2012": (
        "YEQFKK5P",
        "P7XFU7V7",
        "10.1029/2012GL053540",
        [
            "Yamazaki",
            "Ito",
            "Yoshino",
            "Yoneda",
            "Guo",
            "Zhang",
            "Sun",
            "Shimojuku",
            "Tsujino",
            "Kunimoto",
            "Higo",
            "Funakoshi",
        ],
        2012,
        "P-V-T equation of state for ε-iron up to 80 GPa and 1900 K using the Kawai-type high pressure apparatus equipped with sintered diamond anvils",
        "Geophysical Research Letters 39, L20308",
    ),
    "sakai2014": (
        "SHIU4GRC",
        "LZ2CTDCN",
        "10.1016/j.pepi.2013.12.010",
        ["Sakai", "Takahashi", "Nishitani", "Mashino", "Ohtani", "Hirao"],
        2014,
        "Equation of state of pure iron and Fe0.9Ni0.1 alloy up to 3 Mbar",
        "Physics of the Earth and Planetary Interiors 228, 114–126",
    ),
    "dewaele2006": (
        "9UAMT6VU",
        "AQ8NJAAW",
        "10.1103/PhysRevLett.97.215504",
        ["Dewaele", "Loubeyre", "Occelli", "Mezouar", "Dorogokupets", "Torrent"],
        2006,
        "Quasihydrostatic Equation of State of Iron above 2 Mbar",
        "Physical Review Letters 97, 215504",
    ),
}


def ref(source):
    _, _, doi, authors, year, title, journal = SOURCES[source]
    return dict(authors=authors, year=year, title=title, source=journal, doi=doi)


def dump(path, value):
    path.write_text(
        json.dumps(value, indent=1, ensure_ascii=False, allow_nan=False) + "\n"
    )


def number(token):
    m = re.fullmatch(r"([\d.]+)(?:\s*[~(]([\d.]+)[!)]?)?", token.strip())
    if not m:
        raise ValueError(token)
    value = float(m[1])
    error = ""
    if m[2]:
        error = (
            float(m[2])
            if "." in m[2]
            else int(m[2]) * 10 ** (-len(m[1].partition(".")[2]))
        )
    return value, error


def put_number(row, name, token):
    row[name], row[name + "_error"] = number(token)


def write_rows(stem, rows):
    p = DATA / (stem + ".csv")
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with p.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return p


def extract(zotero):
    from pypdf import PdfReader

    provenance = {}

    def pdf(key):
        p = next((zotero / key).glob("*.pdf"))
        provenance[key] = {
            "filename": p.name,
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        }
        return PdfReader(p)

    # Whitespace normalization joins uncertainty tokens split across PDF text lines.
    y = []
    run = None
    for page_no, page in enumerate(pdf("UDTQ8W5D").pages, 1):
        text = re.sub(r"\s+", " ", page.extract_text())
        pattern = r"M\d{4}|(\d{3,4})\s+" + r"\s+".join([r"([\d.]+\s*\(\d+\))"] * 5)
        for m in re.finditer(pattern, text):
            if m[0].startswith("M"):
                run = m[0]
                continue
            row = dict(
                run=run,
                source_page=page_no,
                source_row=len(y) + 1,
                temperature_k=float(m[1]),
            )
            for name, token in zip(
                [
                    "au_volume_a3",
                    "pressure_gpa",
                    "a_angstrom",
                    "c_angstrom",
                    "volume_a3_conventional_cell",
                ],
                m.groups()[1:],
            ):
                put_number(row, name, token)
            y.append(row)
    write_rows("iron-yamazaki-2012-table-s1", y)
    # PDF coordinate-based columns preserve the missing P1/P4/P6/P7 cells.
    import pdfplumber

    p = next((zotero / "SHIU4GRC").glob("*.pdf"))
    pdf("SHIU4GRC")
    s = []
    with pdfplumber.open(p) as doc:
        for page_no in (3, 4):
            page = doc.pages[page_no - 1]
            words = page.extract_words(x_tolerance=1, y_tolerance=2)
            starts = [w for w in words if re.fullmatch(r"(?:FNC|F10N)\w+", w["text"])]
            for start in starts:
                same = [w for w in words if abs(w["top"] - start["top"]) < 2]
                row = {
                    "run": start["text"],
                    "material": "iron"
                    if start["text"].startswith("FNC")
                    else "fe09ni01_hcp",
                    "source_page": page_no,
                    "source_row": len(s) + 1,
                }
                # Column boundaries proportional to A4 source page width.
                edges = [85, 127, 168, 209, 247, 292, 328, 364, 400, 436, 476, 512]
                if page_no == 4:
                    edges = [x + 9.8 for x in edges]
                names = [
                    "pressure_p1_gpa",
                    "pressure_p4_gpa",
                    "pressure_p6_gpa",
                    "pressure_p7_gpa",
                    "temperature_k",
                    "a_angstrom",
                    "c_angstrom",
                    "volume_a3_conventional_cell",
                    "nacl_a_angstrom",
                    "nacl_volume_a3",
                    "mgo_a_angstrom",
                    "mgo_volume_a3",
                ]
                for i, name in enumerate(names):
                    tokens = [
                        w["text"]
                        for w in same
                        if edges[i]
                        <= w["x0"]
                        < (edges[i + 1] if i + 1 < len(edges) else page.width)
                    ]
                    token = "".join(tokens).rstrip("a")
                    if token:
                        put_number(row, name, token)
                    else:
                        row[name] = row[name + "_error"] = ""
                s.append(row)
    write_rows("sakai-2014-table2-pvt", s)
    b = []
    page = pdf("G37FEUSQ").pages[2].extract_text()
    pat = r"\s+".join([r"([\d.]+~[\d.]+!)"] * 4)
    for m in re.finditer(pat, page):
        row = dict(source_table="I" if len(b) < 24 else "II", source_row=len(b) + 1)
        for name, token in zip(
            [
                "particle_velocity_km_s",
                "shock_velocity_km_s",
                "pressure_gpa",
                "density_g_cm3",
            ],
            m.groups(),
        ):
            put_number(row, name, token)
        b.append(row)
    write_rows("iron-brown-2000-tables-i-ii", b)
    for v in SOURCES.values():
        pdf(v[0])
    for key in ["H9BL483X", "EQ4VCIFW", "EKR3V3ZM"]:
        p = next(p for p in (zotero / key).iterdir() if not p.name.startswith("."))
        provenance[key] = {
            "filename": p.name,
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        }
    dump(
        DATA / "iron-source-papers-zotero.json",
        {
            "retrieved_on": DATE,
            "attachments": provenance,
            "items": {k: {"item_key": v[1], "doi": v[2]} for k, v in SOURCES.items()},
        },
    )
    print("Extracted rows", len(y), len(s), len(b))


def read_rows(stem):
    with (DATA / (stem + ".csv")).open() as f:
        return list(csv.DictReader(f))


def dataset(stem, source, location, ids, notes):
    p = DATA / (stem + ".csv")
    rows = read_rows(stem)
    columns = []
    for name in rows[0]:
        base = name.removesuffix("_error")
        quantity = "source_field"
        unit = "dimensionless"
        role = "flag"
        if "pressure_" in base:
            quantity = "pressure"
            unit = "GPa"
            role = "value"
        elif "volume_" in base:
            quantity = "volume"
            unit = "angstrom^3/conventional_unit_cell"
            role = "value"
        elif base == "temperature_k":
            quantity = "temperature"
            unit = "K"
            role = "value"
        elif "angstrom" in base:
            quantity = "lattice_parameter"
            unit = "angstrom"
            role = "value"
        elif "velocity" in base:
            quantity = "velocity"
            unit = "km/s"
            role = "value"
        elif base == "density_g_cm3":
            quantity = "density"
            unit = "g/cm^3"
            role = "value"
        col = dict(name=name, quantity=quantity, unit=unit, role=role)
        if name.endswith("_error"):
            col.update(role="uncertainty", of=base)
        columns.append(col)
    return dict(
        identifier=stem.replace("-", "_"),
        kind="experimental_measurements",
        description=f"{len(rows)} primary table rows. " + notes,
        reference=ref(source),
        source_location=location,
        source_url="https://doi.org/" + SOURCES[source][2],
        columns=columns,
        resource=dict(
            path="datasets/" + p.name,
            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            media_type="text/csv",
        ),
        used_by_eos_records=ids,
        notes=notes,
        license="CC0-1.0",
        license_url="https://creativecommons.org/publicdomain/zero/1.0/",
        license_scope="Contributor-created factual transcription and arrangement only; excludes publisher PDF and all third-party rights.",
        provenance=dict(
            type="primary_pdf_table_transcription",
            source_manifest="datasets/iron-source-papers-zotero.json",
            transcribed_on=DATE,
            quality_control={"rows": len(rows)},
        ),
    )


def record(
    identifier,
    source,
    typ,
    pars,
    errors,
    location,
    datasets,
    prange,
    notes,
    fixed=(),
    scale="Source pressure scale",
):
    reference = ref(source)
    r = dict(
        identifier=identifier,
        label=f"{reference['authors'][0]} et al. ({reference['year']}), {identifier.split(str(reference['year']) + '_')[-1]}",
        reference=reference,
        record_kind="published",
        equation_kind="isothermal",
        volume_basis=dict(
            kind="formula_units",
            formula_units=2,
            molar_mass_g_mol=55.845 if not identifier.startswith("fe09") else 56.12984,
        ),
        eos=dict(
            type=typ,
            model={
                "BM3": "birch_murnaghan_3",
                "Vinet": "vinet",
                "LinearUsUpHugoniot": "linear_us_up_hugoniot",
            }[typ],
            parameters=pars,
        ),
        parameter_errors={k: errors.get(k) for k in pars},
        parameter_error_confidence=None,
        fixed_parameters=list(fixed),
        temperature_ref=300,
        experimental_pressure_range_gpa=prange,
        experimental_temperature_range_k=[300, 300],
        pressure_range_status="reported_exactly",
        validity=dict(
            pressure_gpa=prange,
            temperature_k=[300, 300],
            notes=[
                "Observation envelope; extensions to inner-core conditions are extrapolations."
            ],
        ),
        notes=notes,
        pressure_calibration=dict(
            status="partially_resolved",
            audit_date=DATE,
            methods=[
                dict(
                    kind="equation_of_state",
                    reference=reference,
                    source_location=location,
                    scope=scale,
                )
            ],
            recalculation=dict(
                status="not_possible",
                notes="Primary published pressures retained. " + scale,
            ),
        ),
        parameter_provenance={
            k: location + ": " + ("fixed/adopted" if k in fixed else "published fit")
            for k in pars
        },
        scientific_validation=dict(
            status="primary_source_validated",
            audit_date=DATE,
            note="Primary publication equations, coefficients and scope checked; regression parity is evaluated separately.",
            verified_fields=[
                "equation",
                "parameters",
                "units",
                "phase",
                "reference_state",
                "validity",
            ],
            primary_source_check=dict(
                doi=reference["doi"],
                access_url="https://doi.org/" + reference["doi"],
                locations=[location],
                finding=notes,
            ),
            primary_data_check=dict(
                status="bundled" if datasets else "not_available",
                audit_date=DATE,
                dataset_identifiers=datasets,
                source_locations=[location],
                finding="Primary observations bundled with errors and calibrant values where published."
                if datasets
                else "The 188 input P-V-T points are not tabulated in the attached publication.",
            ),
        ),
    )
    return r


def thermal(r, params, errors, fixed, temps=(300, 1900), kind="Dewaele2006"):
    r["equation_kind"] = "thermal"
    r["thermal"] = dict(
        type=kind,
        model="dewaele_2006" if kind == "Dewaele2006" else "thermal_reference_state",
        parameters=params,
        parameter_errors={k: errors.get(k) for k in params},
        parameter_error_confidence=None,
        fixed_parameters=fixed,
    )
    r["experimental_temperature_range_k"] = list(temps)
    r["validity"]["temperature_k"] = list(temps)
    r["parameter_provenance"] = {
        "reference_isotherm": r["parameter_provenance"],
        "thermal_correction": {
            k: "Primary thermal equation and parameter table; "
            + ("fixed/adopted" if k in fixed else "fitted")
            for k in params
        },
    }
    return r


def debye(g, gi, b, t):
    return dict(
        Tr=300,
        theta0=t,
        gamma0=g,
        gamma_inf=gi,
        beta=b,
        anharmonic_a=3.7e-5,
        anharmonic_m=1.87,
        electronic_e=1.95e-4,
        electronic_g=1.339,
        n=1,
    )


FIXED = ["Tr", "anharmonic_a", "anharmonic_m", "electronic_e", "electronic_g", "n"]
# Table 3/4: each tuple is V0,error,K_BM,error,Kp_BM,error,K_Vinet,error,Kp_Vinet,error.
SAKAI = {
    "iron": [
        (22.18, 0.20, 179.6, 2.2, 4.91, 0.12, 174.7, 2.2, 5.37, 0.12),
        (22.35, None, 164.8, 2.2, 5.18, 0.13, 160.1, 2.1, 5.65, 0.12),
        (22.468, None, 155.3, 2.2, 5.37, 0.14, 150.8, 2.0, 5.84, 0.12),
        (22.16, 0.19, 184.2, 2.2, 4.78, 0.12, 179.2, 2.2, 5.24, 0.12),
        (22.14, 0.19, 185.0, 2.3, 4.94, 0.12, 180.0, 2.3, 5.40, 0.13),
        (22.07, 0.20, 195.1, 2.4, 4.93, 0.12, 189.9, 2.4, 5.39, 0.13),
        (22.10, 0.17, 196.0, 2.5, 4.93, 0.13, 190.8, 2.4, 5.39, 0.13),
        (22.21, 0.23, 171.7, 2.4, 5.49, 0.15, 167.2, 2.2, 5.93, 0.13),
    ],
    "fe09ni01_hcp": [
        (22.18, 0.13, 195.3, 1.6, 4.37, 0.05, 189.6, 1.7, 4.81, 0.06),
        (22.36, None, 179.8, 1.6, 4.59, 0.05, 173.9, 1.6, 5.07, 0.06),
        (22.468, None, 171.0, 1.7, 4.73, 0.06, 165.2, 1.5, 5.23, 0.06),
        (22.16, 0.13, 199.4, 1.6, 4.28, 0.05, 193.9, 1.7, 4.70, 0.06),
        (22.14, 0.13, 200.5, 1.7, 4.42, 0.05, 194.7, 1.7, 4.86, 0.06),
        (22.04, 0.13, 213.0, 1.8, 4.41, 0.05, 207.0, 1.8, 4.85, 0.06),
        (22.10, 0.14, 211.5, 1.8, 4.42, 0.05, 205.4, 1.9, 4.86, 0.06),
        (22.17, 0.15, 191.2, 1.8, 4.82, 0.06, 185.1, 1.7, 5.31, 0.06),
    ],
}
SCALES = {
    "p1": "NaCl-B2 calibrated to Matsui (2009) Pt",
    "p2": "NaCl-B2 calibrated to Fei (2007) Pt",
    "p3": "NaCl-B2 calibrated to Dorogokupets and Oganov (2007) Pt",
    "p4": "NaCl-B2 calibrated to Tange (2009) MgO, Dorfman (2012)",
    "p5": "NaCl-B2 calibrated to Yokoo (2009) Pt",
    "p6": "NaCl-B2 calibrated to Holmes (1989) Pt",
}


def build():
    iron = json.loads((MATERIALS / "iron.eosmat").read_text())
    alloy_path = MATERIALS / "fe09ni01_hcp.eosmat"
    if alloy_path.exists():
        alloy = json.loads(alloy_path.read_text())
    else:
        alloy = {
            k: copy.deepcopy(v)
            for k, v in iron.items()
            if k not in ["eos_records", "datasets"]
        }
        alloy.update(
            identifier="fe09ni01_hcp",
            name="Fe0.9Ni0.1 alloy (hcp)",
            formula="Fe0.9Ni0.1",
            phase="hcp",
            eos_records=[],
            datasets=[],
            notes="Sakai et al. (2014) Fe0.9Ni0.1 hcp alloy. Mixed Fe/Ni occupancy on the hcp 2c site. Reference lattice is the observed F10N06_012 cell at 27.1 GPa, 300 K, not an ambient cell.",
            lattice=dict(a=2.4295, b=None, c=3.917, alpha=90.0, beta=90.0, gamma=120.0),
        )
        alloy["atom_sites"][0]["occupancy"] = 0.9
        ni = copy.deepcopy(alloy["atom_sites"][0])
        ni.update(element="Ni", occupancy=0.1)
        alloy["atom_sites"].append(ni)
        alloy["structure_provenance"] = {
            "reference": ref("sakai2014"),
            "source_location": "Table 2, F10N06_012; Section 3.1 hcp phase assignment",
            "notes": "Conventional hcp Wyckoff 2c ideal substitutional alloy representation; occupancies follow nominal atomic composition.",
        }
    new = {"iron": [], "fe09ni01_hcp": []}
    ds = "sakai_2014_table2_pvt"
    for mat, table in SAKAI.items():
        for spec, vals in zip(
            ["p1", "p1_v0_mao", "p1_v0_dewaele", "p2", "p3", "p4", "p5", "p6"], table
        ):
            v, ev, *rest = vals
            for typ, offset in [("BM3", 0), ("Vinet", 4)]:
                k, ek, kp, ekp = rest[offset : offset + 4]
                note = (
                    "Table "
                    + ("3" if mat == "iron" else "4")
                    + " source fit. V0 from a preliminary g-G plot is held during the K0/K0_prime fit; alternate Mao/Dewaele V0 choices are explicitly fixed. Printed ± confidence and covariance are not specified. Low-pressure NaCl-B1 rows use Brown (1999); other scale variants retain distinct NaCl-B2 reductions."
                )
                r = record(
                    f"{mat}_sakai_2014_{spec}_{typ.lower()}",
                    "sakai2014",
                    typ,
                    dict(V0=v, K0=k, K0_prime=kp),
                    dict(V0=ev, K0=ek, K0_prime=ekp),
                    "Tables 1–4; Sections 2 and 3.3–3.4",
                    [ds],
                    [
                        24.2 if mat == "iron" else 27.1,
                        304.7 if mat == "iron" else 295.8,
                    ],
                    note,
                    ["V0"],
                    SCALES[spec[:2]],
                )
                scale = spec[:2]
                if scale in ("p1", "p4", "p6"):
                    selected = [
                        x
                        for x in read_rows("sakai-2014-table2-pvt")
                        if x["material"] == mat and float(x["temperature_k"]) == 300
                    ]
                    ps = [
                        float(x["pressure_" + scale + "_gpa"] or x["pressure_p1_gpa"])
                        for x in selected
                    ]
                    r["experimental_pressure_range_gpa"] = [min(ps), max(ps)]
                    r["validity"]["pressure_gpa"] = [min(ps), max(ps)]
                else:
                    r["pressure_range_status"] = "reported_qualitatively"
                    r["validity"]["notes"].append(
                        "No exact pressure extrema for this scale are tabulated; bounds use the envelope of published P1/P4/P6 measurements."
                    )
                new[mat].append(r)
    yds = "iron_yamazaki_2012_table_s1"
    for typ, v, ev, k, ek, kp, t, et in [
        ("BM3", 22.15, 0.05, 202, 7, 4.5, 1173, 62),
        ("Vinet", 22.17, 0.06, 196, 8, 4.8, 1168, 61),
    ]:
        r = record(
            f"iron_yamazaki_2012_{typ.lower()}_thermal",
            "yamazaki2012",
            typ,
            dict(V0=v, K0=k, K0_prime=kp),
            dict(V0=ev, K0=ek, K0_prime=0.2),
            "Equations (1)–(4), Table 1, Table S1",
            [yds],
            [19.8, 82.9],
            "Full Debye, anharmonic and electronic pressure, reference-subtracted at 300 K. The power-law gamma is represented exactly by gamma_inf=0, beta=q. Table 1 errors are standard deviations; thermal correction constants adopted from Alfè et al. (2001).",
            scale="Au, Tsuchiya (2003); original scale, not later Fei (2007) re-reduction.",
        )
        thermal(
            r,
            debye(3.2, 0, 0.8, t),
            dict(gamma0=0.2, beta=0.3, theta0=et),
            FIXED + ["gamma_inf"],
        )
        r["parameter_error_confidence"] = r["thermal"]["parameter_error_confidence"] = (
            0.682689492137
        )
        new["iron"].append(r)
    # Exact reciprocal quadratic K(T) represented by the existing cubic law with cubic coefficient zero.
    b1, b2, b3 = 0.005973, 1.38e-6, 4.6e-10
    k = 1 / (b1 + b2 * 300 + b3 * 300**2)
    ek = (
        k * k * (0.000046**2 + (300 * 0.017e-6) ** 2 + (300**2 * 0.080e-10) ** 2) ** 0.5
    )
    r = record(
        "iron_dubrovinsky_2000_bm3_thermal",
        "dubrovinsky2000",
        "BM3",
        dict(V0=6.73 * 2e24 / NA, K0=k, K0_prime=5.81),
        dict(V0=0.01 * 2e24 / NA, K0=None, K0_prime=0.06),
        "Equation (1), Table II",
        [],
        [18, 305],
        "K0 derived exactly from 1/(b1+b2*300+b3*300²). No covariance is available to propagate a rigorous K0 uncertainty. Table II gives one expansivity value; constant alpha=6.93e-5/K is the reconstructable law. No dhcp EOS is provided. Original 109+79 P-V-T rows are not tabulated.",
        scale="Pt, Holmes et al. (1989), with Nellis et al. (1988) also cited.",
    )
    thermal(
        r,
        dict(Tr=300, alpha0=6.93e-5, dK_dT=0, beta1=b2, beta2=b3, beta3=0),
        dict(alpha0=0.37e-5, beta1=0.017e-6, beta2=0.080e-10),
        ["Tr", "dK_dT", "beta3"],
        (300, 1700),
        "AlphaKT",
    )
    r["thermal"]["bulk_modulus_law"] = "reciprocal_cubic"
    r["source_compressibility_coefficients"] = {
        "b1": b1,
        "b1_error": 0.000046,
        "b2": b2,
        "b2_error": 0.017e-6,
        "b3": b3,
        "b3_error": 0.080e-10,
        "conditional_uncorrelated_K0_error_gpa": ek,
    }
    new["iron"].append(r)
    # Five thermal variants use the P4 cold curve and P7 MgO hot pressures (Table 8).
    for name, g, gi, b, t, errs, fixed in [
        ("type1", 2.169, 2.169, 1, 1162, dict(gamma0=0.207, theta0=204), ["beta"]),
        (
            "type2",
            2.609,
            0,
            1.309,
            600,
            dict(gamma0=0.350, beta=0.679, theta0=474),
            ["gamma_inf"],
        ),
        (
            "type4_1",
            2.882,
            1.087,
            3.289,
            577,
            dict(gamma0=0.394, gamma_inf=0.487, theta0=471),
            ["beta"],
        ),
        (
            "type4_2",
            2.883,
            0.968,
            3.289,
            417,
            dict(gamma0=0.374, gamma_inf=0.242),
            ["beta", "theta0"],
        ),
        (
            "type4_3",
            1.875,
            1.305,
            1.161,
            417,
            dict(beta=0.237),
            ["gamma0", "gamma_inf", "theta0"],
        ),
    ]:
        r = copy.deepcopy(
            next(r for r in new["fe09ni01_hcp"] if r["identifier"].endswith("p4_bm3"))
        )
        r.update(
            identifier=f"fe09ni01_hcp_sakai_2014_{name}_thermal",
            label=f"Sakai et al. (2014), Fe0.9Ni0.1 {name} thermal EOS (P4 + P7)",
        )
        thermal(r, debye(g, gi, b, t), errs, FIXED + fixed, (300, 2300))
        r["notes"] = (
            "Tables 6–8; Section 3.5. P4 BM3 reference isotherm, P7 MgO thermal pressures. Anharmonic/electronic constants adopted from pure Fe. "
            + (
                "Type 1 constant gamma: gamma_inf tied to gamma0 and theta=theta0*(V/V0)^(-gamma0), following the Type 4 limiting definition in Section 3.5; printed Table 6 q is not independently tabulated. "
                if name == "type1"
                else ""
            )
            + "Full cold-curve parameters and their uncertainty retained, fixed in thermal regression."
        )
        r["scientific_validation"]["primary_source_check"].update(
            locations=["Tables 4, 6–8; Section 3.5"], finding=r["notes"]
        )
        hot = [
            x
            for x in read_rows("sakai-2014-table2-pvt")
            if x["material"] == "fe09ni01_hcp" and float(x["temperature_k"]) > 300
        ]
        r["experimental_pressure_range_gpa"] = [
            min(float(x["pressure_p7_gpa"]) for x in hot),
            max(float(x["pressure_p7_gpa"]) for x in hot),
        ]
        r["validity"]["pressure_gpa"] = [27.1, 290.7]
        r["validity"]["notes"].append(
            "Thermal observations cover 92.7–148.9 GPa; broader cold curve does not imply high-temperature data at those pressures."
        )
        r["pressure_calibration"]["methods"][0]["scope"] = (
            "P4 at 300 K, P7 Tange (2009) MgO at high temperature."
        )
        new["fe09ni01_hcp"].append(r)
    # Dewaele preferred Vinet thermal already exists; distinct BM3 alternative is absent.
    r = record(
        "iron_dewaele_2006_bm3",
        "dewaele2006",
        "BM3",
        dict(V0=22.468, K0=165.0, K0_prime=4.97),
        dict(V0=0.024, K0=None, K0_prime=0.04),
        "Table I, BM row; atomic volume 11.234(12) doubled to conventional cell",
        ["iron_dewaele_2006_epaps_compression"],
        [17, 197],
        "Distinct Table I Birch–Murnaghan alternative to the existing Vinet thermal record. K0=165 GPa fixed. Source table errors are at 95% confidence. Atomic-volume error is doubled, not copied from later secondary tables.",
        ["K0"],
        "Ruby/W, Dorogokupets–Oganov (2006).",
    )
    r["parameter_error_confidence"] = 0.95
    new["iron"].append(r)
    # Brown provides a nominal-iron standard measured on low-carbon steel, not pure Fe.
    r = record(
        "iron_brown_2000_linear_hugoniot",
        "brown2000",
        "LinearUsUpHugoniot",
        dict(V0=2 * 55.845 / (NA * 1e-24 * 7.850), rho0=7.850, c0=3.935, s=1.578, P0=0),
        {},
        "Table III linear fit to all 37 Tables I–II points",
        ["iron_brown_2000_tables_i_ii"],
        [39.6, 191.8],
        "Published full-range linear Us-up fit, operationally restricted to the solid-iron subset used by Fei (2016), below 200 GPa. Samples are ~99% low-carbon steel (up to 0.7% Mn), not pure iron; density 7.850 g/cm³ must not be replaced by pure-Fe density. No shock temperatures reported. The full 442.1 GPa data, quadratic fit and fit range are retained in the source audit; high-pressure mixed-phase path is not represented as an hcp EOS.",
        ["V0", "rho0", "P0"],
        "Shock velocities and Rankine–Hugoniot relations; no static calibrant.",
    )
    r.update(
        label="Brown et al. (2000), nominal-iron steel linear Hugoniot, solid subset",
        sample_composition="Approximately 99% low-carbon steel; nominal Fe volume normalization. See source for impurity limits.",
        equation_kind="hugoniot",
        loading_path="principal",
        branch_kind="transformed",
        initial_state=dict(
            phase="bcc",
            material_identifier="fe",
            temperature_k=300.0,
            temperature_basis="Room-temperature normalization, not a reported shock temperature.",
            pressure_gpa=0.0,
            density_g_cm3=7.850,
        ),
        branch_domain=dict(
            particle_velocity_km_s=[0.948, 2.885],
            kind="experimental_coverage",
            boundary_status="inferred",
            notes=[
                "Operational subset below 200 GPa used by Fei 2016; not a phase boundary."
            ],
        ),
        source_fit_domain=dict(
            particle_velocity_km_s=[0.948, 4.894], pressure_gpa=[39.6, 442.1]
        ),
        source_quadratic_fit=dict(
            c0=3.691, s=1.788, q=-0.038, shock_velocity_rmse_km_s=0.039
        ),
    )
    r["validity"]["notes"].append(
        "Temperature bounds describe only the normalized initial state, not a constant-temperature shock path. The sample is steel, not pure Fe."
    )
    r["pressure_calibration"]["status"] = "not_applicable"
    r["pressure_calibration"]["recalculation"]["status"] = "not_applicable"
    new["iron"].append(r)
    for mat, doc in [("iron", iron), ("fe09ni01_hcp", alloy)]:
        existing = {r["identifier"] for r in doc["eos_records"]}
        for r in new[mat]:
            if r["identifier"] not in existing:
                doc["eos_records"].append(r)
            else:
                doc["eos_records"] = [
                    r if old["identifier"] == r["identifier"] else old
                    for old in doc["eos_records"]
                ]
        ids = [r["identifier"] for r in new[mat] if "sakai_2014" in r["identifier"]]
        datasets = [
            dataset(
                "sakai-2014-table2-pvt",
                "sakai2014",
                "Table 2, journal pp. 116–117",
                ids,
                "Both compositions, all published pressure columns and marker lattice values retained. Missing cells remain blank. Filter material and source pressure scale explicitly. No P2/P3/P5 pressures are tabulated.",
            )
        ]
        if mat == "iron":
            datasets += [
                dataset(
                    "iron-yamazaki-2012-table-s1",
                    "yamazaki2012",
                    "Table S1",
                    [
                        r["identifier"]
                        for r in new[mat]
                        if "yamazaki" in r["identifier"]
                    ],
                    "All runs and temperatures retained on original Tsuchiya (2003) Au scale.",
                ),
                dataset(
                    "iron-brown-2000-tables-i-ii",
                    "brown2000",
                    "Tables I and II, journal p. 5497",
                    ["iron_brown_2000_linear_hugoniot"],
                    "Full 37-row shock dataset; 1-sigma errors. Temperature is absent, not inferred. Initial density 7.850 g/cm³.",
                ),
            ]
        for d in datasets:
            if d["identifier"] not in {x["identifier"] for x in doc["datasets"]}:
                doc["datasets"].append(d)
            else:
                doc["datasets"] = [
                    d if old["identifier"] == d["identifier"] else old
                    for old in doc["datasets"]
                ]
        if mat == "iron":
            d = next(
                d
                for d in doc["datasets"]
                if d["identifier"] == "iron_dewaele_2006_epaps_compression"
            )
            if "iron_dewaele_2006_bm3" not in d["used_by_eos_records"]:
                d["used_by_eos_records"].append("iron_dewaele_2006_bm3")
        dump(MATERIALS / (mat + ".eosmat"), doc)
    nacl_path = MATERIALS / "nacl_b2.eosmat"
    nacl = json.loads(nacl_path.read_text())
    r = record(
        "nacl_b2_sakai_2014_yokoo_pt_bm3",
        "sakai2014",
        "BM3",
        dict(V0=38.34, K0=45.18, K0_prime=4.22),
        dict(V0=4.69, K0=0.48, K0_prime=0.02),
        "Section 2, journal p. 115: new P5 NaCl-B2 calibration",
        [],
        [30, 330],
        "New Yokoo-Pt calibration of Sakai (2011a) NaCl-B2 observations. V0=38.34(4.69) A³ from g-G analysis, then fixed. The unusually large printed V0 uncertainty is retained verbatim. Exact new calibration pressure extrema and original paired NaCl/Pt rows are not given in this article; the validity envelope is approximate.",
        ["V0"],
        "NaCl-B2 recalibrated against Yokoo et al. (2009) Pt.",
    )
    r["volume_basis"] = dict(nacl["eos_records"][0]["volume_basis"])
    r["pressure_range_status"] = "reported_qualitatively"
    del r["experimental_pressure_range_gpa"]
    r["scientific_validation"]["primary_data_check"]["finding"] = (
        "Original Sakai (2011a) paired NaCl/Pt calibration observations are not attached to this 2014 paper. Table 2 NaCl values are downstream measurements, not independent P5 calibration constraints."
    )
    nacl["eos_records"] = [
        x for x in nacl["eos_records"] if x["identifier"] != r["identifier"]
    ] + [r]
    dump(nacl_path, nacl)
    print("Source records", {**{m: len(v) for m, v in new.items()}, "nacl_b2": 1})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zotero-storage", type=Path)
    args = parser.parse_args()
    if args.zotero_storage:
        extract(args.zotero_storage)
    build()
