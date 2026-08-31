"""Apply Peritheos's primary-source audit to migrated ``.eosmat`` records.

This script is intentionally separate from the mechanical Dioptas importer.
The importer establishes file provenance; this script records the independent
scientific review.  A record is executable only after the cited primary paper
or its official supplement has established the equation, every stored EOS
parameter, units/reference state, phase, and the represented data range.

Run from the repository root::

    python scripts/apply_primary_source_audit.py

The generated audit report is committed so releases do not depend on network
access.  BurnMan and Pytheos source code were not consulted during this audit.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

AUDIT_DATE = "2026-08-31"
ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos" / "data" / "materials"
REPORT = ROOT / "peritheos" / "data" / "primary-source-audit.json"

VERIFIED_FIELDS = [
    "equation",
    "parameters",
    "units",
    "reference_state",
    "phase",
    "published_uncertainties",
    "validity",
]


def source(url: str, locations: list[str], note: str = "") -> dict[str, Any]:
    """Construct one compact primary-evidence description."""
    result: dict[str, Any] = {"access_url": url, "locations": locations}
    if note:
        result["finding"] = note
    return result


# Every entry below was compared with the original article or an author/
# publisher-hosted copy.  The location is deliberately precise enough for a
# later reviewer to repeat the comparison without using another software
# catalog as an authority.
VALIDATED_SOURCES: dict[str, dict[str, Any]] = {
    "10.1007/bf00203299": source(
        "https://doi.org/10.1007/BF00203299",
        ["300 K fit reported in the article; EOS table and P-V-T discussion"],
    ),
    "10.1007/s004100050471": source(
        "https://doi.org/10.1007/s004100050471",
        ["Abstract and EOS results for the single-crystal 298 K data"],
    ),
    "10.1016/j.cageo.2016.06.002": source(
        "https://doi.org/10.1016/j.cageo.2016.06.002",
        ["Equations 1-12", "Table 1", "Figure 2 pressure-calculation table"],
        "Table 1 parameters and the Figure 2 calculations were recomputed "
        "independently; the record names the generic cold curve and "
        "multi-oscillator thermal model rather than the paper.",
    ),
    "10.1016/j.epsl.2011.02.025": source(
        "https://doi.org/10.1016/j.epsl.2011.02.025",
        ["Equation-of-state section", "Table 1", "Figures 2 and 3"],
        "Fixed values in italics and fitted uncertainties were preserved.",
    ),
    "10.1016/j.jallcom.2005.04.008": source(
        "https://people.iith.ac.in/kanchana/publications/2005/17.pdf",
        [
            "Equation (1), pages 57-58",
            "section 4, pages 58-59",
            "Tables 1 and 2",
            "Figure 3",
        ],
        "The stored records are the room-temperature experimental BM3 fits for "
        "the cubic fluorite phases, not the separate SIC-LSD calculations.",
    ),
    "10.1016/j.gca.2003.12.007": source(
        "https://lweb.cfa.harvard.edu/~lzeng/papers/frank_fei_hu_2004.pdf",
        [
            "Equation (2), pages 2782-2783",
            "section 3.1, pages 2782-2783",
            "Table 1",
            "Figure 2",
        ],
        "The stored record is the simultaneous three-parameter 300 K BM3 fit; "
        "its molar V0 is converted to the two-formula-unit cubic cell.",
    ),
    "10.1016/j.lithos.2015.03.017": source(
        "https://doi.org/10.1016/j.lithos.2015.03.017",
        ["Primary article abstract and reported BM3 fits for pyrope and almandine"],
    ),
    "10.1016/j.pepi.2021.106786": source(
        "https://doi.org/10.1016/j.pepi.2021.106786",
        ["Equation-of-state discussion", "Table 1"],
    ),
    "10.1016/0022-3697(81)90074-3": source(
        "https://www.researchgate.net/publication/256267527_Bulk_Moduli_and_High-Pressure_Crystal_Structures_of_Rutile-Type_Compounds",
        [
            "EOS fit discussion, journal page 144",
            "Tables 1 and 2, journal pages 145-146",
            "Table 3, journal page 147",
        ],
        "The author-uploaded primary article reports K0 with K0'=7 assumed. "
        "The SnO2 and GeO2 V0 values are calculated from the Table 3 ambient "
        "lattice constants; no EOS-fit covariance or V0 error is published.",
    ),
    "10.1016/s0031-9201(00)00154-0": source(
        "https://duffy.princeton.edu/sites/g/files/toruqf616/files/shim_pepi_00.pdf",
        ["section 3", "Table 1", "Table 2, pages 334-335"],
        "The adopted unconstrained BM3 fit fixes V0=45.58 A^3 in Table 2, "
        "while the abstract explicitly reports V0=45.58+/-0.05 A^3. The "
        "reported V0 uncertainty is retained even though V0 was fixed in "
        "this fit.",
    ),
    "10.1029/2000jb900318": source(
        "https://doi.org/10.1029/2000JB900318",
        ["Primary abstract", "constrained third-order Birch-Murnaghan fit"],
    ),
    "10.1029/2005gl024955": source(
        "https://doi.org/10.1029/2005GL024955",
        ["Equation-of-state section", "Table 1"],
    ),
    "10.1029/2008jb005900": source(
        "https://doi.org/10.1029/2008JB005900",
        ["Equation-of-state formulation", "Table 4", "combined shock/static fit"],
        "The stored 300 K record is the static BM3 component of the published model.",
    ),
    "10.1029/2009gl038107": source(
        "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009GL038107",
        ["paragraphs 5-9", "Table 1", "Figures 1 and 2"],
        "The paper supplies the BM3 reference isotherm and Mie-Gruneisen-Debye "
        "parameters. Theta0=814 K is explicitly fixed; K0=169.2(9) GPa is "
        "fixed in the fit but its source uncertainty is retained.",
    ),
    "10.1029/96gl03769": source(
        "https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/96GL03769",
        ["Equation (1)", "Table 1", "Thermoelastic Equation of State, pages 6-7"],
        "The 300 K BM3 reference isotherm uses V0=403.32(8) A^3, "
        "K0=124.5(4.0) GPa, and K0'=5 assumed. The full paper also reports "
        "the separate high-temperature coefficients, which this static record does not claim.",
    ),
    "10.1029/2011jb008988": source(
        "https://doi.org/10.1029/2011JB008988",
        ["EOS formulation", "Table 2, BM3 and Vinet rows"],
    ),
    "10.1029/2018je005582": source(
        "https://doi.org/10.1029/2018JE005582",
        ["Equation-of-state section", "Table 1, B3 and B1 rows"],
    ),
    "10.1029/2019jb017853": source(
        "https://doi.org/10.1029/2019JB017853",
        ["Thermal EOS formulation", "Table 1, CaCl2 and seifertite rows"],
        "The stored BM2 records are the 300 K reference isotherms of the MGD fits.",
    ),
    "10.1029/jb095ib13p21737": source(
        "https://doi.org/10.1029/JB095iB13p21737",
        ["Primary abstract", "combined hcp-Fe BM3 fit to 300 GPa"],
        "The molar volume was converted to the conventional two-atom hcp cell.",
    ),
    "10.1029/jb093ib12p15279": source(
        "https://www.researchgate.net/publication/248791203_Static_compression_and_equation_of_state_of_CaO_to_135_Mbar",
        ["Equation (5)", "Tables 1-4", "pages 15283-15285"],
        "The author-uploaded primary article reports separate BM3 fits for "
        "B1 and B2 CaO; Mbar and molar-volume values are converted explicitly.",
    ),
    "10.1029/jb094ib12p17889": source(
        "https://okayama.elsevierpure.com/en/publications/stability-and-equation-of-state-of-casiosub3sub-perovskite-to-134/",
        ["Author/institution record", "primary abstract", "47-point BM3 fit"],
        "The authors report V0=45.37(8) A^3 and K0=281(4) GPa with K0'=4 assumed.",
    ),
    "10.1038/ngeo2370": source(
        "https://doi.org/10.1038/ngeo2370",
        ["Main text", "Figure 1 caption and EOS fit"],
    ),
    "10.1038/s41598-019-51037-8": source(
        "https://doi.org/10.1038/s41598-019-51037-8",
        ["Equation-of-state section", "Table 3"],
    ),
    "10.1038/s41598-019-51931-1": source(
        "https://doi.org/10.1038/s41598-019-51931-1",
        ["Equation-of-state discussion", "Table 8"],
        "Per-atom volumes were converted using each phase's conventional-cell Z.",
    ),
    "10.1038/s41598-022-10523-2": source(
        "https://doi.org/10.1038/s41598-022-10523-2",
        ["Equation-of-state results", "EOS parameter table"],
    ),
    "10.1038/s41598-024-78006-0": source(
        "https://doi.org/10.1038/s41598-024-78006-0",
        ["Equation-of-state results", "Table 3"],
    ),
    "10.1038/s43246-025-00963-4": source(
        "https://doi.org/10.1038/s43246-025-00963-4",
        ["Model description", "Table 1, combined experimental fit"],
        "The stored record is the 300 K BM3 subset of a combined thermal fit.",
    ),
    "10.1038/srep22652": source(
        "https://doi.org/10.1038/srep22652",
        ["Methods and EOS results", "EOS parameter table"],
    ),
    "10.1039/d4nr00093e": source(
        "https://doi.org/10.1039/D4NR00093E",
        ["Equation 4", "Table 3", "Figure 5"],
    ),
    "10.1063/1.3644969": source(
        "https://doi.org/10.1063/1.3644969",
        ["Static EOS formulation", "wurtzite and rocksalt EOS table"],
    ),
    "10.1063/1.4863300": source(
        "https://doi.org/10.1063/1.4863300",
        ["Equation 2", "Table I, 300 K Vinet fit", "Figure 3"],
    ),
    "10.1063/1.4894421": source(
        "https://doi.org/10.1063/1.4894421",
        ["Primary abstract", "BM2 P-V-T fits for ice VI and ice VII"],
        "Molar volumes were converted with Z=10 for ice VI and Z=2 for ice VII.",
    ),
    "10.1063/1.344177": source(
        "https://www.researchgate.net/profile/John-Moriarty-7/publication/224513875_The_equation_of_state_of_platinum_to_660_GPa_66_Mbar/links/0deec51d6192cc65e2000000/The-equation-of-state-of-platinum-to-660-GPa-66-Mbar.pdf",
        [
            "Equations (7)-(12)",
            "Table IV",
            "Figures 4 and 5",
            "summary, page 2966",
        ],
        "The stored pressure scale is the theoretical 300 K universal (Vinet) "
        "isotherm qualified against shock data, plus the paper's approximate "
        "constant thermal-pressure extension below 2000 K.",
    ),
    "10.1063/5.0179469": source(
        "https://doi.org/10.1063/5.0179469",
        ["Equation 1", "Table II, own 300 K data", "Figure 4"],
    ),
    "10.1073/pnas.0609013104": source(
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC1890468/",
        ["Equations 1-3", "Table 1", "Ne BM3 values in main text"],
        "Equation 3 uses the published variable-exponent Debye-temperature law.",
    ),
    "10.1080/08957950212807": source(
        "https://doi.org/10.1080/08957950212807",
        ["Primary abstract", "BM3 fit and zero-pressure molar volume"],
        "The molar volume was converted to the conventional unit cell.",
    ),
    "10.1088/0953-8984/13/11/303": source(
        "https://electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20Physics%20Condensed%20Matter/2001%20Volume%2013/0953-8984_13_11_303.pdf",
        [
            "section 2, ambient lattice constants and pressure methods",
            "section 4.1, pages 2449-2450",
            "Figure 2, Birch-Murnaghan fit through 46 GPa",
        ],
        "The experimental EOS is fitted to relative volume. The measured ambient "
        "hexagonal subcell was converted to the equivalent four-formula-unit "
        "orthorhombic conventional-cell volume used by the record.",
    ),
    "10.1098/rsta.2022.0331": source(
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC10069115/",
        ["Abstract", "Figure 1 caption", "Table 1"],
        "The paper is internally inconsistent: the abstract says third-order "
        "Birch-Murnaghan, while the Figure 1 caption says second-order despite "
        "reporting K0'=3.3(1). Peritheos retains BM3 and records the conflict.",
    ),
    "10.1103/physrevb.103.014101": source(
        "https://doi.org/10.1103/PhysRevB.103.014101",
        ["Equation-of-state section", "EOS fit table"],
    ),
    "10.1103/physrevb.39.12598": source(
        "https://doi.org/10.1103/PhysRevB.39.12598",
        ["Page 12599", "Murnaghan fit paragraph", "Figure 1"],
    ),
    "10.1103/physrevb.30.6045": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.30.6045/fulltext",
        [
            "Equation (1)-(5), pages 6047-6048",
            "Table I, B2 compression data",
            "page 6048, second-order finite-strain fit",
        ],
        "The reported B2/B1 zero-pressure volume ratio was converted using "
        "the ambient B1 formula-unit volume; the second-order fit fixes K0'=4.",
    ),
    "10.1103/physrevb.49.12540": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.49.12540/fulltext",
        [
            "Equation (1), page 12543",
            "Table I, 300 K row",
            "page 12544, stated 0-10 GPa limit",
        ],
        "The EOS and crystallographic data are for D2O ice VIII, not ordinary H2O; "
        "the material metadata was corrected accordingly. The authors state that "
        "no pressure/volume error can be quoted because it depends on the adopted EOS.",
    ),
    "10.1103/physrevb.54.5": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.54.5/fulltext",
        ["page 6, hcp-volume paragraph", "Figure 3"],
        "The primary paper fits hcp cobalt with the Vinet EOS; the migrated BM3 "
        "model label was corrected.",
    ),
    "10.1103/physrevb.59.8526": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.59.8526/fulltext",
        ["Equation (1), page 8528", "Table II", "conclusion, page 8528"],
        "The fixed V0=10.865 cm3/mol was converted to a two-atom bcc cell.",
    ),
    "10.1103/physrevb.70.012101": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.70.012101/fulltext",
        ["Table I", "page 2 EOS discussion", "Figure 2"],
        "The paper fixes relative V0 to unity; the absolute hcp cell follows its stated ambient lattice parameters.",
    ),
    "10.1103/physrevb.70.094112": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.70.094112/fulltext",
        ["Vinet equation below Figure 2", "Table II, revised-ruby-scale rows"],
        "The records use the revised pressure scale (PR') rather than the classical ruby-scale fits.",
    ),
    "10.1103/physrevb.73.224119": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.73.224119/fulltext",
        ["Table I", "Table II", "section III.B, pages 224119-5 to -6"],
        "The adopted BM3 fit is the 0-71.5 GPa quasihydrostatic helium subset; "
        "the wider errors span the alternative fit through 134 GPa.",
    ),
    "10.1103/physrevb.75.214104": source(
        "https://doi.org/10.1103/PhysRevB.75.214104",
        ["Equation 2", "Table II", "Figures 2 and 3"],
    ),
    "10.1103/physrevb.77.094106": source(
        "https://doi.org/10.1103/PhysRevB.77.094106",
        ["Vinet formulation", "Table I", "300 K neon data"],
    ),
    "10.1103/physrevb.78.104102": source(
        "https://doi.org/10.1103/PhysRevB.78.104102",
        ["Vinet formulation", "Table II, Ni and Ag rows"],
    ),
    "10.1103/physrevb.85.214105": source(
        "https://doi.org/10.1103/PhysRevB.85.214105",
        ["Vinet formulation", "Table II, KCl/KBr B1 and B2 rows"],
    ),
    "10.1103/physrevb.88.064107": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.88.064107/fulltext",
        ["section IV, page 064107-2", "Table I", "Figure 3", "section V"],
        "All three experimental Vinet parameters were fitted; the quoted errors "
        "are 95% confidence intervals. The migrated constrained-fit metadata was corrected.",
    ),
    "10.1103/physrevb.90.134105": source(
        "https://link.aps.org/accepted/10.1103/PhysRevB.90.134105",
        ["300 K EOS section", "Vinet fit table for hcp and bcc Mg"],
    ),
    "10.1103/physrevb.91.134108": source(
        "https://doi.org/10.1103/PhysRevB.91.134108",
        ["Section IV", "Table III, experimental alpha and omega Vinet rows"],
        "Per-atom V0 values were converted to conventional-cell volumes.",
    ),
    "10.1103/physrevlett.113.265504": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.113.265504/fulltext",
        ["Vinet discussion", "Figure 3", "Table I, experimental FeH2 and FeH3 rows"],
        "Molar volumes were converted to the represented conventional cells.",
    ),
    "10.1103/physrevlett.74.1371": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.74.1371/fulltext",
        ["Figure 2", "page 1372, third-order Birch-Murnaghan fit"],
        "The represented record is the authors' nonhydrostatic static-compression "
        "fit, not their preferred combined fourth-order hydrostatic isotherm.",
    ),
    "10.1103/zp3m-kjpc": source(
        "https://doi.org/10.1103/zp3m-kjpc",
        ["Vinet formulation", "Table I, beta row", "Table II, alpha and omega rows"],
        "The beta-phase V0 is explicitly an extrapolated fixed fit parameter.",
    ),
    "10.1126/science.235.4789.668": source(
        "https://www.science.org/doi/epdf/10.1126/science.235.4789.668",
        ["Equations (1)-(3)", "Tables 1 and 2", "Figures 1 and 2"],
        "The primary PDF identifies the sample as (Mg0.88Fe0.12)SiO3 and "
        "reports the unconstrained Eulerian finite-strain BM3 fit.",
    ),
    "10.1140/epjb/e2003-00034-6": source(
        "https://arxiv.org/abs/cond-mat/0210013",
        ["ambient-structure paragraph", "Table of BM EOS parameters", "EOS discussion"],
        "The authors' preprint reports a=5.9240(4) A and the B1-PbS BM3 fit "
        "K0=51.0(1.2) GPa, K0'=4.3(9). V0 and its uncertainty are propagated "
        "from the published cubic lattice parameter.",
    ),
    "10.1107/s0021889897000861": source(
        "https://doi.org/10.1107/S0021889897000861",
        ["Primary abstract", "third-order Birch-Murnaghan fit"],
    ),
    "10.1107/s0108768197005739": source(
        "https://doi.org/10.1107/S0108768197005739",
        ["Primary article", "hcp Ar EOS fit and reported uncertainty"],
    ),
    "10.2138/am-2000-11-1229": source(
        "https://doi.org/10.2138/am-2000-11-1229",
        ["Primary abstract", "K0'=4 constrained BM3 alternative"],
    ),
    "10.2138/am-2002-0419": source(
        "https://web.gps.caltech.edu/~jackson/pdf/Angel_Jackson_AmMin02_87_558.pdf",
        ["Abstract", "Table 1", "Summary, page 561"],
        "The recommended BM3 parameters combine two compression experiments, "
        "Brillouin data, and high-pressure ultrasonic data.",
    ),
    "10.2138/am-2002-0701": source(
        "https://www.researchgate.net/publication/231182054_Thermal_equations_of_state_for_B1_and_B2_KCl",
        ["Equation BE1", "Table 2", "Table 3", "B2 KCl Thermal EOS, pages 807-810"],
        "The author-uploaded primary article's preferred B2 fit is the bold "
        "Table 3 row. It explicitly declines individual parameter errors "
        "because of covariance, but reports alpha0*K0=0.0275+/-0.0009 kbar/K.",
    ),
    "10.2138/am-2002-2-316": source(
        "https://doi.org/10.2138/am-2002-2-316",
        ["Equation-of-state section", "EOS parameter table"],
    ),
    "10.2138/am-2003-2-307": source(
        "https://doi.org/10.2138/am-2003-2-307",
        ["Primary abstract", "combined single-crystal and powder BM3 fit"],
    ),
    "10.2138/am-1997-7-805": source(
        "https://rruff.info/doclib/am/vol82/AM82_682.pdf",
        [
            "Abstract",
            "Table 1, measured P-V data",
            "Unit-cell parameters and equation of state, pages 684-685",
        ],
        "The unconstrained BM3 fit is represented. The alternative fit with "
        "K0'=4 fixed is not stored in this record.",
    ),
    "10.2138/am-2020-7279": source(
        "https://doi.org/10.2138/am-2020-7279",
        ["Table 1, static BM3 row", "Table 2, MGD reference-isotherm row"],
    ),
    "10.2138/am.2006.2347": source(
        "https://doi.org/10.2138/am.2006.2347",
        ["Primary abstract", "orthorhombic FeS VI constrained BM3 fit"],
    ),
    "10.2138/am.2008.2942": source(
        "https://doi.org/10.2138/am.2008.2942",
        ["EOS method", "goethite and epsilon-FeOOH parameter tables"],
    ),
    "10.2138/am.2011.3775": source(
        "https://doi.org/10.2138/am.2011.3775",
        ["Equation-of-state section", "BM2 fit table"],
    ),
    "10.2138/am.2014.4526": source(
        "https://doi.org/10.2138/am.2014.4526",
        ["Primary abstract", "BM3 fit"],
    ),
    "10.3103/s1063457614010092": source(
        "https://arxiv.org/abs/1407.1677",
        ["Vinet formulation", "EOS comparison table", "Figure 2"],
        "The author manuscript states that V0=93.2061 A^3 was fixed.",
    ),
    "10.3390/min9110684": source(
        "https://doi.org/10.3390/min9110684",
        ["Vinet formulation", "Table 2", "Figure 3"],
    ),
    "10.5194/ejm-33-519-2021": source(
        "https://doi.org/10.5194/ejm-33-519-2021",
        ["Equation-of-state section", "Table 3"],
    ),
}


# Some older articles and technical reports have no DOI but do have stable,
# unambiguous publisher or institutional copies.  Key these by record ID so a
# title-only match can never promote a different migrated record.
VALIDATED_RECORD_SOURCES: dict[str, dict[str, Any]] = {
    "coesite_levien_1981_bm3_1": source(
        "https://msaweb.org/AmMin/AM66/AM66_324.pdf",
        [
            "Abstract",
            "Table 7, P-V data",
            "page 328, Elasticity section, Birch-Murnaghan fit",
        ],
        "The conventional-cell V0 is calculated from the paper's ambient "
        "lattice parameters; K0 and K0' are the reported unweighted fit.",
    ),
    "lead_fcc_fortes_2019_bm4_1": source(
        "https://epubs.stfc.ac.uk/manifestation/40740885/RAL-TR-2019-002.pdf",
        [
            "Equations 1-5",
            "sections 2.1-2.3",
            "Table 1, report page 9",
        ],
        "The stored record is the 300 K isotherm of the report's temperature-dependent BM4 model.",
    ),
    "zircon_hazen_1979_bm3_1": source(
        "https://msaweb.org/AmMin/AM64/AM64_196.pdf",
        [
            "Table 1",
            "page 198, Zircon compressibility section",
        ],
        "The reported K0=227(2) GPa fit assumes K0'=6.5; the migrated BM2 "
        "representation was corrected to BM3 with K0' fixed.",
    ),
}


# A few valid DOI groups have one intentionally nonvalidated record.  This is
# keyed by record ID so DOI-level evidence can never accidentally promote it.
FORCED_DEFERRED: dict[str, str] = {
    "aragonite_martinez_1996_bm3_1": (
        "The primary paper was inspected, but the migrated AlphaKT component "
        "uses the reported mean 298-1000 K expansivity as a constant. The "
        "published global high-temperature BM3 model instead uses equations "
        "(2), (4), and (5), including alpha(T)=alpha0+alpha1*T. Peritheos does "
        "not silently treat that approximation as the published model; the "
        "record remains non-executable until the exact thermal law is represented."
    ),
    "gold_anderson_1989_bm3_1": (
        "The primary paper was identified, but the complete recommended "
        "parameter derivation and its reference-state conventions were not "
        "available in an accessible primary copy during this audit."
    ),
}


DEFERRED_BY_DOI: dict[str, str] = {
    "10.1007/s002690000108": (
        "The Ono et al. primary abstract establishes the cubic phases and fitted "
        "bulk moduli, but an accessible primary copy was not found that establishes "
        "the stored reference volumes, full fit constraints, and uncertainties."
    ),
    "10.1016/0022-3697(91)90181-x": (
        "The primary abstract reports the B2/B1 zero-pressure volume ratio, K0, "
        "and fixed K0'=4, but the stored absolute B2 V0 additionally depends on an "
        "ambient B1 volume that is not established by the accessible primary text."
    ),
    "10.1016/0022-3697(93)90106-2": (
        "The primary abstract establishes the NiS BM3 K0 and K0' errors, but no "
        "accessible primary copy was found that establishes the stored V0, its "
        "reference state, and the complete measured pressure range."
    ),
    "10.1016/0031-9201(92)90063-2": (
        "The authors' institutional bibliography and primary abstract establish "
        "K0=161 GPa with K0'=4 assumed, but the full primary tables needed to "
        "verify V0=1513.1 A^3, its constraint, and the exact data range were inaccessible."
    ),
    "10.1016/s0012-821x(00)00033-9": (
        "The primary abstract establishes a second-order Birch-Murnaghan K0=134(5) "
        "GPa for Fe-bearing phase D, but the full primary table needed to verify "
        "the stored V0 and represented range was not accessible."
    ),
    "10.1016/s0022-3697(98)00296-0": (
        "The primary abstract establishes the NiO shock-derived K0 and K0' and a "
        "147 GPa upper bound, but an accessible primary copy was not found for the "
        "stored V0, uncertainty convention, or exact static-isotherm reduction."
    ),
    "10.1016/s0038-1098(99)00322-1": (
        "The identified primary fit is Vinet over combined bcc/fcc lithium data, "
        "whereas the migrated record is BM3. The complete primary article was not "
        "accessible, so Peritheos does not change the equation or promote the record."
    ),
    "10.1029/2000gl012606": (
        "The primary abstract reports the cementite BM3 K0 and K0', but an "
        "accessible full primary copy was not found to establish the stored V0, "
        "its error, and the exact represented pressure interval."
    ),
    "10.1029/94jb00127": (
        "The primary abstract establishes the CsCl Eulerian finite-strain K0 and "
        "K0' errors, but the full primary density/reference-volume convention and "
        "all stored metadata were not accessible."
    ),
    "10.1029/jb078i035p08470": (
        "The accessible primary abstract does not establish all stored B1-SrO "
        "reference-volume, equation-order, fit-constraint, and error metadata."
    ),
    "10.1029/jb079i008p01165": (
        "The primary abstract establishes the magnetite K0, fixed K0', and pressure "
        "range, but an accessible full primary copy was not found for the stored "
        "reference volume and complete uncertainty provenance."
    ),
    "10.1029/jb086ib12p11773": (
        "The primary abstract verifies the B1-B2 transition and a 25 GPa B2 data "
        "span, but does not expose the B2 reference volume, K0, fit equation, or "
        "their uncertainty/constraint conventions; the full article was inaccessible."
    ),
    "10.1029/jb094ib03p03037": (
        "The primary abstract verifies the Birch-Murnaghan framework and composition "
        "trends, but an accessible full primary copy was not found that establishes "
        "the stored Mg0.4Fe0.6O V0, fit coefficients, constraints, and their errors."
    ),
    "10.1029/rf002p0029": (
        "The cited DOI is a handbook chapter on thermal expansion rather than an "
        "unambiguous original source for the stored FeO static EOS parameter set. "
        "The record is deferred until its actual primary EOS source is identified."
    ),
    "10.1063/1.1726610": (
        "The accessible primary abstract does not establish all stored CoO V0, "
        "Birch-Murnaghan coefficients, errors, and validity metadata."
    ),
    "10.1063/1.351203": (
        "The cited Hixson-Fritz work is a shock-compression standard. An accessible "
        "primary copy was not found that justifies the migrated standalone 300 K "
        "BM3 reduction, its parameters, errors, and stated range."
    ),
    "10.1088/0953-8984/5/33/010": (
        "The primary article is an ab initio calculation. The accessible primary "
        "metadata do not establish the migrated experimental-structure V0 together "
        "with the theoretical BM3 coefficients and a defensible validity interval."
    ),
    "10.1103/fxgq-96sg": (
        "APS lists publication on 24 April 2026, but the accepted manuscript is "
        "under CHORUS embargo until 24 April 2027. The primary equation and "
        "tables could not be inspected; catalog values are therefore not executable."
    ),
}


# Equation (1) and Table 1 of Sokolova et al. require both the atom count n
# and atomic/effective atomic number Z.  The mechanical Dioptas migration did
# not carry all of these model inputs, so the audit restores them directly
# from the primary table rather than treating a software catalog as authority.
SOKOLOVA_COMPOSITION: dict[str, tuple[float, float]] = {
    "aluminum_sokolova_2016_holzapfel_2": (1.0, 13.0),
    "copper_sokolova_2016_holzapfel_2": (1.0, 29.0),
    "diamond_sokolova_2016_holzapfel_3": (1.0, 6.0),
    "gold_sokolova_2016_holzapfel_4": (1.0, 79.0),
    "mgo_sokolova_2016_holzapfel_4": (2.0, 10.34),
    "molybdenum_sokolova_2016_holzapfel_2": (1.0, 42.0),
    "niobium_sokolova_2016_holzapfel_2": (1.0, 41.0),
    "platinum_sokolova_2016_holzapfel_3": (1.0, 78.0),
    "silver_sokolova_2016_holzapfel_2": (1.0, 47.0),
    "tantalum_sokolova_2016_holzapfel_3": (1.0, 73.0),
    "tungsten_sokolova_2016_holzapfel_4": (1.0, 74.0),
}

SUN_SILICA_RECORDS = {
    "silica_cacl2_sun_2019_bm2_1",
    "seifertite_sun_2019_bm2_1",
}

BEZACIER_ICE_RECORDS = {
    "ice_vi_bezacier_2014_bm2_1",
    "ice_vii_bezacier_2014_bm2_1",
}


def set_primary_parameter(component: dict[str, Any], name: str, value: float) -> None:
    """Set one primary-sourced component parameter and its metadata."""
    component.setdefault("parameters", {})[name] = value
    component.setdefault("parameter_errors", {})[name] = None
    fixed = component.setdefault("fixed_parameters", [])
    if name not in fixed:
        fixed.append(name)


def append_correction(record: dict[str, Any], correction: dict[str, Any]) -> None:
    """Append one correction exactly once, making the audit idempotent."""
    corrections = record.setdefault("audit_corrections", [])
    path = correction["path"]
    corrections[:] = [item for item in corrections if item.get("path") != path]
    corrections.append(correction)


def restore_primary_model_inputs(record: dict[str, Any]) -> None:
    """Restore inputs omitted by the mechanical interchange migration."""
    identifier = record["identifier"]
    if record.get("thermal", {}).get("type") == "AlphaKT":
        record["thermal"]["model"] = "thermal_reference_state"
    if identifier in SOKOLOVA_COMPOSITION:
        n, atomic_number = SOKOLOVA_COMPOSITION[identifier]
        for component, name, value in (
            ("eos", "n", n),
            ("eos", "Z", atomic_number),
            ("thermal", "n", n),
        ):
            equation_component = record[component]
            set_primary_parameter(equation_component, name, value)
            append_correction(
                record,
                {
                    "path": f"{component}.parameters.{name}",
                    "source_value": None,
                    "value": value,
                    "reason": (
                        "Required model input omitted by the Dioptas migration; "
                        "restored from Sokolova et al. (2016) Table 1."
                    ),
                    "primary_reference": {
                        "doi": "10.1016/j.cageo.2016.06.002",
                        "location": "Equation (1) and Table 1",
                    },
                },
            )

    if identifier in SUN_SILICA_RECORDS:
        set_primary_parameter(record["thermal"], "n", 3.0)
        append_correction(
            record,
            {
                "path": "thermal.parameters.n",
                "source_value": None,
                "value": 3.0,
                "reason": (
                    "The molar Debye energy requires the three atoms in the SiO2 "
                    "formula unit; this model input was omitted in migration."
                ),
                "primary_reference": {
                    "doi": "10.1029/2019JB017853",
                    "location": "thermal EOS formulation and Table 1",
                },
            },
        )

    if identifier in BEZACIER_ICE_RECORDS:
        set_primary_parameter(record["thermal"], "Tr", 300.0)
        append_correction(
            record,
            {
                "path": "thermal.parameters.Tr",
                "source_value": None,
                "value": 300.0,
                "reason": (
                    "The thermal reference temperature in equations (2)-(3) was "
                    "stored only at record level by the migration."
                ),
                "primary_reference": {
                    "doi": "10.1063/1.4894421",
                    "location": "equations (1)-(3) and Table II",
                },
            },
        )

    if identifier in {"zircon_hazen_1979_bm2_1", "zircon_hazen_1979_bm3_1"}:
        record["identifier"] = "zircon_hazen_1979_bm3_1"
        record["label"] = "Hazen and Finger (1979), K0' fixed at 6.5"
        record["eos"] = {
            "type": "BM3",
            "parameters": {"V0": 260.79, "K0": 227.0, "K0_prime": 6.5},
            "model": "birch_murnaghan_3",
        }
        record["parameter_errors"] = {
            "V0": 0.04,
            "K0": 2.0,
            "K0_prime": None,
        }
        record["fixed_parameters"] = ["K0_prime"]
        record["experimental_pressure_range_gpa"] = [0.0, 4.81]
        record["notes"] = (
            "Hazen and Finger report V0=260.79(4) A^3 and K0=227(2) GPa "
            "for Murnaghan and Birch-Murnaghan fits with K0'=6.5 assumed. "
            "The paper's fitted range is 0-48.1 kbar."
        )
        for correction in (
            {
                "path": "identifier",
                "source_value": "zircon_hazen_1979_bm2_1",
                "value": "zircon_hazen_1979_bm3_1",
                "reason": "The primary fit requires an explicit fixed K0'=6.5.",
            },
            {
                "path": "eos",
                "source_value": "BM2(V0=260.803, K0=227.0)",
                "value": "BM3(V0=260.79, K0=227.0, K0_prime=6.5)",
                "reason": (
                    "The primary paper reports V0=260.79(4) A^3 and "
                    "K0=227(2) GPa with K0'=6.5 assumed."
                ),
            },
            {
                "path": "fixed_parameters.K0_prime",
                "source_value": None,
                "value": 6.5,
                "reason": "K0'=6.5 is assumed, not fitted, in the primary paper.",
            },
            {
                "path": "parameter_errors.V0",
                "source_value": None,
                "value": 0.04,
                "reason": "The primary paper reports V0=260.79(4) A^3.",
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [0.0, 4.8],
                "value": [0.0, 4.81],
                "reason": "The highest measured pressure in Table 1 is 48.1 kbar.",
            },
        ):
            correction["primary_reference"] = {
                "location": "Hazen and Finger (1979), page 198",
                "url": "https://msaweb.org/AmMin/AM64/AM64_196.pdf",
            }
            append_correction(record, correction)

    if identifier == "alumina_dewaele_2013_vinet_1":
        record["parameter_errors"] = {
            "V0": 0.52,
            "K0": 7.6,
            "K0_prime": 0.15,
        }
        record["fixed_parameters"] = []
        record["experimental_pressure_range_gpa"] = [0.56, 165.0]
        record["temperature_ref"] = 300.0
        record["label"] = "Dewaele and Torrent (2013), Vinet"
        record["notes"] = (
            "Quasihydrostatic 300 K Vinet fit to corundum data from 0.56 to "
            "165 GPa. V0=255.45+/-0.52 A^3, K0=254.1+/-7.6 GPa, and "
            "K0'=4.00+/-0.15 were all fitted; errors are 95% confidence intervals."
        )
        for correction in (
            {
                "path": "parameter_errors",
                "source_value": {"V0": None, "K0": None, "K0_prime": None},
                "value": dict(record["parameter_errors"]),
                "reason": "The primary paper reports 95% confidence intervals for all three fit coefficients.",
            },
            {
                "path": "fixed_parameters",
                "source_value": ["K0_prime"],
                "value": [],
                "reason": "K0'=4.00+/-0.15 was fitted, not constrained to 4.",
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [0.0, 165.0],
                "value": [0.56, 165.0],
                "reason": "Table I spans measured pressures from 0.56 to 165 GPa.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1103/physrevb.88.064107",
                "location": "section IV, Table I, Figure 3, and section V",
            }
            append_correction(record, correction)

    if identifier == "cobalt_hcp_fujihisa_1996_bm3_1":
        record["identifier"] = "cobalt_hcp_fujihisa_1996_vinet_1"
        record["label"] = "Fujihisa and Takemura (1996), hcp Vinet"
        record["eos"]["type"] = "Vinet"
        record["eos"]["model"] = "vinet"
        record["fixed_parameters"] = ["V0"]
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature hcp-cobalt Vinet fit to 79 GPa: K0=199+/-6 GPa "
            "and K0'=3.6+/-0.2. Volumes were normalized to the ambient hcp "
            "volume, retained here as the fixed conventional-cell V0."
        )
        for correction in (
            {
                "path": "identifier",
                "source_value": "cobalt_hcp_fujihisa_1996_bm3_1",
                "value": "cobalt_hcp_fujihisa_1996_vinet_1",
                "reason": "The primary paper explicitly identifies the fit as Vinet.",
            },
            {
                "path": "eos",
                "source_value": "BM3",
                "value": "Vinet",
                "reason": "The primary paper explicitly identifies the fit as Vinet.",
            },
            {
                "path": "fixed_parameters",
                "source_value": ["V0"],
                "value": ["V0"],
                "reason": "The fitted data are normalized to the ambient hcp volume.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1103/physrevb.54.5",
                "location": "page 6, hcp-volume paragraph and Figure 3",
            }
            append_correction(record, correction)

    if identifier == "niobium_takemura_2006_bm3_1":
        record["eos"]["parameters"]["V0"] = 35.9599
        record["parameter_errors"]["V0"] = 0.0098
        record["experimental_pressure_range_gpa"] = [0.0, 71.5]
        record["notes"] = (
            "Adopted 300 K BM3 fit to the quasihydrostatic helium data through "
            "71.5 GPa: K0=168+/-4 GPa and K0'=3.4+/-0.3. V0 is fixed from "
            "a=3.3007(3) A at atmospheric pressure; the parameter errors span "
            "the alternative fit using data through 134 GPa."
        )
        for correction in (
            {
                "path": "parameter_errors.V0",
                "source_value": None,
                "value": 0.0098,
                "reason": "Propagated from the primary paper's a=3.3007(3) A atmospheric lattice parameter.",
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [0.0, 134.0],
                "value": [0.0, 71.5],
                "reason": "The adopted fit uses the best quasihydrostatic helium subset through 71.5 GPa.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1103/physrevb.73.224119",
                "location": "Table I, Table II, and section III.B",
            }
            append_correction(record, correction)

    if identifier == "tantalum_cynn_1999_bm3_1":
        record["eos"]["parameters"]["V0"] = 36.0835
        record["fixed_parameters"] = ["V0"]
        record["experimental_pressure_range_gpa"] = [0.0, 174.23]
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Ambient-temperature BM3 fit to bcc tantalum through 174.23 GPa: "
            "K0=194.7+/-4.8 GPa and K0'=3.4+/-0.1. The paper fixes "
            "V0=10.865 cm^3/mol, converted here to the conventional two-atom cell."
        )
        for correction in (
            {
                "path": "eos.parameters.V0",
                "source_value": 36.1268,
                "value": 36.0835,
                "reason": "Converted directly from the paper's fixed V0=10.865 cm^3/mol.",
            },
            {
                "path": "fixed_parameters",
                "source_value": [],
                "value": ["V0"],
                "reason": "The primary fit uses the stated ambient molar volume as V0.",
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [0.0, 174.0],
                "value": [0.0, 174.23],
                "reason": "Table II reports the highest compression point at 174.23 GPa.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1103/physrevb.59.8526",
                "location": "Equation (1), Table II, and page 8528 EOS paragraph",
            }
            append_correction(record, correction)

    if identifier == "iceviii_besson_1994_bm3_1":
        record["label"] = "Besson et al. (1994), D2O ice VIII, 300 K BM3"
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Besson et al. give the 300 K D2O ice-VIII molar V0=12.45 "
            "cm^3/mol, converted here to the Z=8 tetragonal cell. Their BM3 fit "
            "gives K0=20.4 GPa and K0'=4.7 for 0-10 GPa. They state that no "
            "pressure/volume error can be quoted because the result depends on "
            "the adopted 300 K EOS, and that it must not be used above 10 GPa."
        )
        append_correction(
            record,
            {
                "path": "material.formula",
                "source_value": "H2O",
                "value": "D2O",
                "reason": "The primary neutron-diffraction EOS is explicitly for deuterated ice VIII.",
                "primary_reference": {
                    "doi": "10.1103/physrevb.49.12540",
                    "location": "title, abstract, section II.B, and Table I",
                },
            },
        )

    if identifier == "magnesite_ross_1997_bm3_1":
        record["eos"]["parameters"]["V0"] = 279.41
        record["parameter_errors"]["V0"] = 0.08
        record["temperature_ref"] = 298.0
        record["notes"] = (
            "Room-temperature unconstrained BM3 fit to single-crystal data "
            "from 0 to 7 GPa: V0=279.41+/-0.08 A^3, K0=117+/-3 GPa, and "
            "K0'=2.3+/-0.7. The paper also reports an alternative BM2 fit "
            "with K0=111+/-1 GPa; this record represents the BM3 result."
        )
        append_correction(
            record,
            {
                "path": "eos.parameters.V0",
                "source_value": 279.28,
                "value": 279.41,
                "reason": "Use the V0 fitted jointly with the stored unconstrained BM3 parameters.",
                "primary_reference": {
                    "doi": "10.2138/am-1997-7-805",
                    "location": "abstract and Unit-cell parameters and equation of state, pages 684-685",
                },
            },
        )
        append_correction(
            record,
            {
                "path": "parameter_errors.V0",
                "source_value": 0.03,
                "value": 0.08,
                "reason": "The primary unconstrained BM3 fit reports V0=279.41(8) A^3.",
                "primary_reference": {
                    "doi": "10.2138/am-1997-7-805",
                    "location": "abstract and Unit-cell parameters and equation of state, pages 684-685",
                },
            },
        )

    if identifier == "molybenum_carbide_mo2c_haines_2001_bm3_1":
        record["eos"]["parameters"]["V0"] = 148.9071
        record["parameter_errors"] = {
            "V0": 0.049,
            "K0": 5.0,
            "K0_prime": 0.3,
        }
        record["experimental_pressure_range_gpa"] = [0.0, 46.0]
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature BM3 fit to relative-volume X-ray diffraction data "
            "through 46 GPa: K0=307+/-5 GPa and K0'=6.2+/-0.3. V0 is calculated "
            "from this sample's measured ambient hexagonal subcell a=3.0128(4) A "
            "and c=4.7357(9) A, converted to the equivalent four-formula-unit "
            "orthorhombic conventional-cell volume; its error is propagated from "
            "the reported lattice-parameter errors."
        )
        for correction in (
            {
                "path": "eos.parameters.V0",
                "source_value": 147.4587,
                "value": 148.9071,
                "reason": (
                    "Use the measured ambient volume of the sample normalized in "
                    "the primary relative-volume EOS, not unrelated literature "
                    "orthorhombic lattice constants."
                ),
            },
            {
                "path": "parameter_errors.V0",
                "source_value": None,
                "value": 0.049,
                "reason": (
                    "Propagated from a=3.0128(4) A and c=4.7357(9) A for "
                    "V=4*(sqrt(3)/2)*a^2*c."
                ),
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [0.0, 50.0],
                "value": [0.0, 46.0],
                "reason": "Section 2 states that the second experimental run reached 46 GPa.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1088/0953-8984/13/11/303",
                "location": "section 2, section 4.1, and Figure 2",
            }
            append_correction(record, correction)

    if identifier in {
        "platinum_holmes_1989_bm3_1",
        "platinum_holmes_1989_vinet_1",
    }:
        record["identifier"] = "platinum_holmes_1989_vinet_1"
        record["label"] = "Holmes et al. (1989), theoretical 300 K Vinet scale"
        record["eos"] = {
            "type": "Vinet",
            "parameters": {
                "V0": 60.4000884,
                "K0": 266.0,
                "K0_prime": 5.81,
            },
            "model": "vinet",
        }
        record["parameter_errors"] = {
            "V0": None,
            "K0": None,
            "K0_prime": None,
        }
        record["fixed_parameters"] = ["V0"]
        record["thermal"] = {
            "type": "LinearThermalPressure",
            "model": "linear_thermal_pressure",
            "parameters": {"Tr": 300.0, "alpha_KT": 0.0069426},
            "parameter_errors": {"Tr": None, "alpha_KT": None},
            "fixed_parameters": [],
        }
        record["temperature_ref"] = 300.0
        record["experimental_pressure_range_gpa"] = [0.0, 550.0]
        record["experimental_temperature_range_k"] = [300.0, 2000.0]
        record["notes"] = (
            "The 300 K pressure scale is Holmes et al.'s theoretical universal "
            "(Vinet) EOS, qualified against shock data to at least 10% in pressure: "
            "V0=101.9 bohr^3/atom (converted to the four-atom fcc cell), K0=266 "
            "GPa, and K0'=5.81. Their approximate finite-temperature extension "
            "adds alpha*K_T*(T-300 K), with alpha=0.261e-4 K^-1 and K_T=266 GPa, "
            "and is stated to be adequate below 2000 K. The 300 K isotherm itself "
            "is represented through 550 GPa; 32-660 GPa is the shock-Hugoniot "
            "range, not the static-isotherm validity range."
        )
        for correction in (
            {
                "path": "identifier",
                "source_value": "platinum_holmes_1989_bm3_1",
                "value": "platinum_holmes_1989_vinet_1",
                "reason": "Equation (11) is explicitly the universal (Vinet) EOS, not BM3.",
            },
            {
                "path": "eos",
                "source_value": "BM3(V0=60.3793086, K0=266, K0_prime=5.81)",
                "value": "Vinet(V0=60.4000884, K0=266, K0_prime=5.81)",
                "reason": (
                    "Use Equation (11) and convert its V0=101.9 bohr^3/atom "
                    "to the four-atom fcc conventional cell."
                ),
            },
            {
                "path": "thermal",
                "source_value": "AlphaKT(alpha0=1.35e-5, dK_dT=0)",
                "value": "LinearThermalPressure(Tr=300, alpha_KT=0.0069426)",
                "reason": (
                    "Equation (12) is the constant thermal-pressure increment "
                    "obtained from Table IV alpha=0.261e-4 K^-1 and K_T=266 GPa."
                ),
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [32.0, 660.0],
                "value": [0.0, 550.0],
                "reason": (
                    "The previous range describes the experimental Hugoniot; the "
                    "represented theoretical 300 K isotherm is stated for 0-550 GPa."
                ),
            },
            {
                "path": "experimental_temperature_range_k",
                "source_value": None,
                "value": [300.0, 2000.0],
                "reason": "The paper states that Equation (12) is adequate below 2000 K.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1063/1.344177",
                "location": "Equations (7)-(12), Table IV, and summary",
            }
            append_correction(record, correction)

    if identifier == "ca_perovskite_shim_2000_bm3_1":
        record["label"] = "Shim et al. (2000), cubic CaSiO3 BM3"
        record["parameter_errors"]["V0"] = 0.05
        record["experimental_pressure_range_gpa"] = [0.0, 108.0]
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature unconstrained BM3 fit for cubic CaSiO3 "
            "perovskite: V0=45.58+/-0.05 A^3 (fixed in the Table 2 fit), "
            "K0=236+/-4 GPa, and "
            "K0'=3.9+/-0.2. Uncertainties are one standard deviation; the "
            "published V0 uncertainty is retained despite its fixed fit "
            "status. Data extend to 108 GPa."
        )
        for correction in (
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [0.0, 134.0],
                "value": [0.0, 108.0],
                "reason": "The primary experiment and Table 1 extend to 108 GPa, not 134 GPa.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1016/s0031-9201(00)00154-0",
                "location": "section 3 and Tables 1-2, pages 334-335",
            }
            append_correction(record, correction)

    if identifier == "ca_perovskite_mao_1989_bm3_2":
        record["label"] = "Mao et al. (1989), cubic CaSiO3 BM3"
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature third-order Birch-Murnaghan fit to 47 P-V "
            "measurements through 134 GPa: V0=45.37+/-0.08 A^3, "
            "K0=281+/-4 GPa, and K0'=4 assumed."
        )

    if identifier == "cao_richet_1988_bm3_1":
        record["label"] = "Richet et al. (1988), CaO B1 BM3"
        record["eos"]["parameters"]["V0"] = 111.3225
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature B1-CaO BM3 fit: molar V0=16.76 cm^3/mol "
            "(converted to the four-formula-unit cell), K0=1.112 Mbar, and "
            "K0'=4.2+/-0.2. The abstract gives K0=1.11+/-0.01 Mbar; the "
            "stored 1 GPa error reflects that precision. V0 is fixed and no "
            "independent V0 error is reported."
        )
        append_correction(
            record,
            {
                "path": "eos.parameters.V0",
                "source_value": 111.3256,
                "value": 111.3225,
                "reason": "Converted from the primary Table 4 molar V0=16.76 cm^3/mol using Z=4.",
                "primary_reference": {
                    "doi": "10.1029/jb093ib12p15279",
                    "location": "Equation (5), Table 4, and abstract",
                },
            },
        )

    if identifier == "cao_b2_richet_1988_bm3_1":
        record["label"] = "Richet et al. (1988), CaO B2 BM3"
        record["experimental_pressure_range_gpa"] = [52.7, 135.0]
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature B2-CaO BM3 fit: molar V0=14.8+/-0.2 "
            "cm^3/mol (converted to the one-formula-unit cell), K0=1.3+/-0.2 "
            "Mbar, and K0'=3.5+/-0.5. B2 observations span 52.7-135 GPa."
        )
        append_correction(
            record,
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [52.7, 134.9],
                "value": [52.7, 135.0],
                "reason": "The paper states an upper pressure of 1.35 Mbar.",
                "primary_reference": {
                    "doi": "10.1029/jb093ib12p15279",
                    "location": "Tables 2-4 and abstract",
                },
            },
        )

    if identifier in {"geo2_rutile_hazen_1981_bm3_1", "sno2_hazen_1981_bm3_1"}:
        is_ge = identifier.startswith("geo2")
        old_v0 = 55.3268 if is_ge else 71.5521
        new_v0 = 55.3318 if is_ge else 71.4957
        new_range = [0.0, 3.70] if is_ge else [0.0, 4.96]
        record["label"] = (
            "Hazen and Finger (1981), rutile GeO2 BM3"
            if is_ge
            else "Hazen and Finger (1981), rutile SnO2 BM3"
        )
        record["eos"]["parameters"]["V0"] = new_v0
        record["experimental_pressure_range_gpa"] = new_range
        record["temperature_ref"] = 293.15
        record["notes"] = (
            f"20 degC BM3 representation with K0' fixed at 7: K0="
            f"{258 if is_ge else 218}+/-{5 if is_ge else 2} GPa. V0={new_v0} "
            "A^3 is calculated from the Table 3 ambient a0 and c0; the paper "
            "does not publish an EOS-fit V0 error or parameter covariance."
        )
        for correction in (
            {
                "path": "eos.parameters.V0",
                "source_value": old_v0,
                "value": new_v0,
                "reason": "Calculated directly from the Table 3 ambient rutile lattice constants.",
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [0.0, 5.0],
                "value": new_range,
                "reason": "Use the highest pressure tabulated for this material.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1016/0022-3697(81)90074-3",
                "location": "Tables 1-3 and EOS discussion, pages 144-147",
            }
            append_correction(record, correction)

    if identifier == "pbs_b1_knorr_2003_bm3_1":
        record["label"] = "Knorr et al. (2003), PbS B1 BM3"
        record["eos"]["parameters"]["V0"] = 207.8955
        record["parameter_errors"] = {
            "V0": 0.0421,
            "K0": 1.2,
            "K0_prime": 0.9,
        }
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature B1-PbS Birch-Murnaghan fit: K0=51.0+/-1.2 "
            "GPa and K0'=4.3+/-0.9. V0=207.8955+/-0.0421 A^3 is calculated "
            "from the reported cubic a=5.9240(4) A. The B1/B16 coexistence "
            "interval ends near 6 GPa."
        )
        for correction in (
            {
                "path": "eos.parameters.V0",
                "source_value": 209.2,
                "value": 207.8955,
                "reason": "Use a^3 from the primary paper's B1-PbS a=5.9240(4) A.",
            },
            {
                "path": "parameter_errors.V0",
                "source_value": 0.1,
                "value": 0.0421,
                "reason": "First-order propagation of the published lattice-parameter error through V=a^3.",
            },
            {
                "path": "parameter_errors.K0",
                "source_value": 1.0,
                "value": 1.2,
                "reason": "The primary EOS table reports K0=51.0(1.2) GPa.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1140/epjb/e2003-00034-6",
                "location": "ambient-structure paragraph, EOS table, and EOS discussion",
            }
            append_correction(record, correction)

    if identifier == "wadsleyite_katsura_2009_bm3_1":
        record["eos"]["parameters"]["V0"] = 538.49
        record["parameter_errors"] = {
            "V0": 0.02,
            "K0": 0.9,
            "K0_prime": 0.1,
        }
        record["thermal"]["parameter_errors"] = {
            "Tr": None,
            "theta0": None,
            "gamma0": 0.02,
            "q": 0.1,
            "n": None,
        }
        record["thermal"]["fixed_parameters"] = ["Tr", "theta0", "n"]
        record["experimental_temperature_range_k"] = [300.0, 2100.0]
        record["notes"] = (
            "Katsura et al. fit 11-20 GPa and 300-2100 K data with a BM3 "
            "reference isotherm plus Mie-Gruneisen-Debye thermal pressure. "
            "Their measured V0=538.49(2) A^3 and K0=169.2(9) GPa are fixed "
            "during the fit; K0'=4.1(1), gamma0=1.64(2), and q=1.5(1) are "
            "fitted, while theta0=814 K is fixed. Fixed source parameters "
            "retain their published uncertainties."
        )
        for correction in (
            {
                "path": "eos.parameters.V0",
                "source_value": 538.185,
                "value": 538.49,
                "reason": "Use the primary paper's measured sample V0=538.49(2) A^3.",
            },
            {
                "path": "parameter_errors",
                "source_value": {"V0": None, "K0": None, "K0_prime": 0.1},
                "value": dict(record["parameter_errors"]),
                "reason": "Restore every explicitly published reference-isotherm error.",
            },
            {
                "path": "thermal.parameter_errors",
                "source_value": None,
                "value": dict(record["thermal"]["parameter_errors"]),
                "reason": "Restore gamma0=1.64(2) and q=1.5(1), and mark fixed thermal inputs.",
            },
            {
                "path": "experimental_temperature_range_k",
                "source_value": None,
                "value": [300.0, 2100.0],
                "reason": "The primary experiment spans 300-2100 K.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1029/2009gl038107",
                "location": "paragraphs 5-9 and Table 1",
            }
            append_correction(record, correction)

    if identifier == "naalsi2o6_zhao_1997_bm3_1":
        record["label"] = "Zhao et al. (1997), jadeite 300 K BM3"
        record["eos"]["parameters"].update({"V0": 403.32, "K0": 124.5})
        record["parameter_errors"] = {"V0": 0.08, "K0": 4.0, "K0_prime": None}
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "The 300 K reference isotherm of Zhao et al.'s high-temperature "
            "BM3 model uses V0=403.32(8) A^3, K0=124.5(4.0) GPa, and K0'=5 "
            "assumed. This record intentionally represents only that reference "
            "isotherm, not the paper's temperature-dependent K and expansion."
        )
        for correction in (
            {
                "path": "eos.parameters.V0",
                "source_value": 401.19,
                "value": 403.32,
                "reason": "Table 1 reports the ambient 300 K cell volume as 403.32(8) A^3.",
            },
            {
                "path": "eos.parameters.K0",
                "source_value": 125.0,
                "value": 124.5,
                "reason": "The EOS fit reports K0=124.5(4.0) GPa; 125 GPa is only its rounded abstract value.",
            },
            {
                "path": "parameter_errors.V0",
                "source_value": None,
                "value": 0.08,
                "reason": "Restore the Table 1 V0=403.32(8) A^3 measurement error.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1029/96gl03769",
                "location": "Equation (1), Table 1, and pages 6-7",
            }
            append_correction(record, correction)

    if identifier == "perovskite_orthorhombic_knittle_1987_bm3_1":
        record["label"] = "Knittle and Jeanloz (1987), (Mg0.88Fe0.12)SiO3 BM3"
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature unconstrained Eulerian finite-strain BM3 fit "
            "for orthorhombic (Mg0.88Fe0.12)SiO3 perovskite from 25 to 127 "
            "GPa: V0=162.77+/-0.39 A^3 (from the reported ambient lattice "
            "parameters), K0=266+/-6 GPa, and K0'=3.9+/-0.4."
        )

    if identifier == "kcl_walker_2002_bm3_2":
        record["temperature_ref"] = 296.15
        record["experimental_pressure_range_gpa"] = [3.18, 8.14]
        record["experimental_temperature_range_k"] = [296.15, 873.15]
        record["thermal"] = {
            "type": "LinearThermalPressure",
            "parameters": {"Tr": 296.15, "alpha_KT": 0.00275},
            "model": "linear_thermal_pressure",
            "parameter_errors": {"Tr": None, "alpha_KT": 0.00009},
            "fixed_parameters": ["Tr"],
        }
        record["audit_corrections"] = [
            correction
            for correction in record.get("audit_corrections", [])
            if correction.get("path") != "thermal.parameters.alpha0"
        ]
        record["notes"] = (
            "Walker et al.'s preferred B2-KCl BE1 fit uses fictive V0=53.53 "
            "A^3, K0=23.7 GPa, and K0'=4.4. Equation BE1 adds the directly "
            "reported alpha0*K0=0.0275+/-0.0009 kbar/K thermal pressure "
            "coefficient relative to the 23 degC isotherm. The paper states "
            "that individual V0, K0, K0', and alpha0 errors are not meaningful "
            "because of parameter correlation, so only the published error of "
            "the identifiable alpha0*K0 product is retained. The represented "
            "data span 3.18-8.14 GPa and 23-600 degC."
        )
        for correction in (
            {
                "path": "thermal",
                "source_value": {
                    "type": "AlphaKT",
                    "parameters": {"alpha0": 0.00018, "dK_dT": 0.0},
                    "model": "thermal_reference_state",
                },
                "value": dict(record["thermal"]),
                "reason": (
                    "Equation BE1 is BM3 plus the additive alpha0*K0*delta-T "
                    "term, not a temperature-dependent V0 and K0 reference "
                    "state. Store the directly reported product and its error "
                    "after converting kbar/K to GPa/K."
                ),
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": [2.0, 8.0],
                "value": [3.18, 8.14],
                "reason": "Table 2 bounds the represented B2 observations at 31.8-81.4 kbar.",
            },
            {
                "path": "experimental_temperature_range_k",
                "source_value": None,
                "value": [296.15, 873.15],
                "reason": "Tables 2-3 span 23-600 degC.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.2138/am-2002-0701",
                "location": "Equation BE1, Tables 2-3, and pages 807-810",
            }
            append_correction(record, correction)


def normalized_doi(reference: object) -> str | None:
    if not isinstance(reference, dict):
        return None
    value = reference.get("doi")
    if not isinstance(value, str) or not value.strip():
        return None
    return value.lower().removeprefix("https://doi.org/").removeprefix("doi:")


def deferred_reason(record: dict[str, Any], doi: str | None) -> str:
    if record["identifier"] in FORCED_DEFERRED:
        return FORCED_DEFERRED[record["identifier"]]
    if doi in DEFERRED_BY_DOI:
        return DEFERRED_BY_DOI[doi]
    if doi is None:
        return (
            "The migrated record does not contain an unambiguous DOI or equivalent "
            "primary-source locator. The EOS equation and every stored parameter "
            "could not be traced to one citable primary table or equation."
        )
    return (
        "The cited primary publication was identified, but an accessible primary "
        "copy did not establish every stored EOS parameter, its equation, units, "
        "reference state, phase assignment, uncertainty convention, and represented "
        "data range. The record is retained for interchange but is not executable."
    )


def audit_record(record: dict[str, Any], material_file: str) -> dict[str, Any]:
    result = dict(record)
    # ``audit_corrections`` is generated entirely by this script. Rebuild it
    # from the current rules so moved or withdrawn corrections cannot remain
    # attached to an unrelated record on subsequent runs.
    result.pop("audit_corrections", None)
    restore_primary_model_inputs(result)
    doi = normalized_doi(result.get("reference"))
    previous = result.get("scientific_validation") or {}
    migration = previous.get("migration_source")

    record_evidence = VALIDATED_RECORD_SOURCES.get(result["identifier"])
    doi_evidence = VALIDATED_SOURCES.get(doi) if doi is not None else None
    if result["identifier"] not in FORCED_DEFERRED and (
        record_evidence is not None or doi_evidence is not None
    ):
        evidence = dict(record_evidence or doi_evidence or {})
        evidence["doi"] = doi
        validation: dict[str, Any] = {
            "status": "primary_source_validated",
            "note": (
                "Independently checked against the cited primary publication or "
                "its official/author manuscript; no external software catalog was "
                "used as scientific authority."
            ),
            "audit_date": AUDIT_DATE,
            "verified_fields": list(VERIFIED_FIELDS),
            "primary_source_check": evidence,
        }
    else:
        reason = deferred_reason(result, doi)
        validation = {
            "status": "deferred",
            "note": reason,
            "audit_date": AUDIT_DATE,
            "unresolved": [
                "complete primary equation/parameter trace",
                "reference-state and unit provenance",
                "published uncertainty/covariance and validity provenance",
            ],
            "primary_source_check": {
                "doi": doi,
                "access_url": f"https://doi.org/{doi}" if doi else None,
                "outcome": "insufficient_primary_evidence",
                "reason": reason,
            },
        }
    if migration is not None:
        validation["migration_source"] = migration
    result["scientific_validation"] = validation

    if result["identifier"] == "graphite_hanfland_1989_murnaghan_1":
        result["parameter_errors"] = dict(result["parameter_errors"])
        result["parameter_errors"]["V0"] = 0.02
        append_correction(
            result,
            {
                "path": "parameter_errors.V0",
                "source_value": None,
                "value": 0.02,
                "reason": "The primary paper reports V0=35.12(2) A^3.",
                "primary_reference": {
                    "doi": "10.1103/physrevb.39.12598",
                    "location": "page 12599, Murnaghan-fit paragraph",
                },
            },
        )

    if result["identifier"] == "b4c_somayazulu_2023_bm3_1":
        result["scientific_validation"]["reported_inconsistencies"] = [
            {
                "field": "eos.order",
                "abstract": "third-order Birch-Murnaghan",
                "figure_1_caption": "second-order Birch-Murnaghan",
                "resolution": (
                    "Retain BM3 because K0'=3.3(1), which is incompatible with "
                    "a conventional BM2 fit where K0'=4 is fixed."
                ),
            }
        ]

    return result


def main() -> None:
    entries: list[dict[str, Any]] = []
    for path in sorted(MATERIALS.glob("*.eosmat")):
        document = json.loads(path.read_text(encoding="utf-8"))
        if document["identifier"] == "iceviii":
            document["name"] = "Ice VIII (D2O)"
            document["formula"] = "D2O"
            document["aliases"] = ["Heavy ice VIII", "Deuterated ice VIII"]
            document["notes"] = (
                "Deuterated tetragonal ice VIII. The isotope identity is part of "
                "the EOS provenance: Besson et al. (1994) measured D2O, not H2O. "
                "The legacy powder lines are retained for interchange because a "
                "publication-verified ordered atomic-site model has not yet been curated."
            )
        records = []
        for record in document["eos_records"]:
            audited = audit_record(record, path.name)
            records.append(audited)
            check = audited["scientific_validation"]
            entries.append(
                {
                    "material": document["identifier"],
                    "file": path.name,
                    "record": audited["identifier"],
                    "label": audited["label"],
                    "doi": normalized_doi(audited.get("reference")),
                    "status": check["status"],
                    "note": check["note"],
                    "primary_source_check": check["primary_source_check"],
                }
            )
        document["eos_records"] = records
        path.write_text(
            json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )

    counts = Counter(entry["status"] for entry in entries)
    if len(entries) != 147:
        raise ValueError(f"Expected 147 EOS records, found {len(entries)}")
    if "pending_primary_source_check" in counts:
        raise ValueError("Primary-source audit left pending records")

    report = {
        "format": "peritheos.primary-source-audit",
        "format_version": 1,
        "audit_date": AUDIT_DATE,
        "policy": {
            "scientific_authority": "primary publications and official supplements",
            "external_catalogs": (
                "BurnMan and Pytheos source code were not inspected or used. "
                "Dioptas supplied migration/file provenance only."
            ),
            "validated_definition": (
                "Equation, all stored parameters, units, reference state, phase, "
                "published errors when available, and represented validity/data "
                "range traced to primary evidence."
            ),
            "deferred_behavior": "retained for interchange; rejected by Material.from_eosmat",
        },
        "summary": {"records": len(entries), **dict(sorted(counts.items()))},
        "records": entries,
    }
    REPORT.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    manifest_path = MATERIALS / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["scientific_validation"] = {
        "audit_date": AUDIT_DATE,
        "report": "../primary-source-audit.json",
        "counts": dict(sorted(counts.items())),
        "policy": (
            "Only primary_source_validated records are executable; deferred records "
            "remain available for lossless catalog interchange."
        ),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
