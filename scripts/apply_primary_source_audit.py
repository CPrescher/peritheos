"""Apply Peritheos's primary-source audit to migrated ``.eosmat`` records.

This script is intentionally separate from the mechanical Dioptas importer.
The importer establishes file provenance; this script records the independent
scientific review.  A record is executable only after the cited primary paper
or its official supplement has established the equation, every stored EOS
parameter, units/reference state, phase, and the represented data range.

Run from the repository root::

    python scripts/apply_primary_source_audit.py

The generated audit report is committed so releases do not depend on network
access. External software catalogs are not treated as scientific authority.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

AUDIT_DATE = "2026-09-01"
CATALOG_AUDIT_DATE = "2026-09-03"
REPORT_AUDIT_DATE = "2026-09-04"
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

DERIVED_REFERENCE_ISOTHERM_RECORDS = {
    "diamond_correa_2008_dewaele_anchored",
    "diamond_benedict_2014_dewaele_anchored",
}

DERIVED_REFIT_RECORDS = {
    "molybenum_carbide_mo2c_haines_2001_bm3_refit",
    "neon_fcc_hemley_1989_bm3_refit",
}

CURRENT_SOURCE_AUDIT_RECORDS = {
    "mgo_b1_luo_2023_vinet_thermal_5",
    "ca_perovskite_caracas_2005_bm3_3",
    "rbcl_b2_campbell_1994_bm3_1",
    "rbcl_b2_campbell_1994_bm3_1",
    "sio2_stv_andr_wang_2012_vinet_mgd_2",
}


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
    "10.1103/physrevb.78.024101": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.78.024101/fulltext",
        [
            "Equations 2-7 and 13-18",
            "Table I, diamond row",
            "Figures 3 and 5-8",
            "Section II.A",
        ],
        "The record implements the complete published diamond Helmholtz branch: "
        "the DFT-GGA Vinet cold curve, logarithmic-moment double-Debye weights "
        "and zero-point energy, and the volume-independent T^2 anharmonic term. "
        "Per-atom volumes and energies were converted to Peritheos molar working "
        "units; public volumes remain eight-atom conventional-cell values.",
    ),
    "10.1103/physrevb.89.224109": source(
        "https://arxiv.org/pdf/1311.4577",
        [
            "Equations 3-7",
            "Table I, diamond column",
            "Section III.A",
            "Figure 6",
        ],
        "The record implements the complete diamond Helmholtz branch, including "
        "the motionless-ion Vinet cold curve, volume-dependent double-Debye "
        "weights and zero-point energy, and the T^2 anharmonic term. Per-atom "
        "volumes and energies were converted to Peritheos molar working units; "
        "public volumes remain eight-atom conventional-cell values.",
    ),
    "10.1007/s002690000108": source(
        "https://doi.org/10.1007/s002690000108",
        [
            "page 620, Table 2, measured cubic-phase P-V-T data",
            "page 621, third-order Birch-Murnaghan equation and 300 K fit",
            "page 621, Table 3 and thermal-expansion discussion",
        ],
        "The stored record is only the 300 K reference isotherm. The paper's "
        "separate 25 GPa thermal-expansion result is documented but is not "
        "silently turned into a full thermal EOS. The notation V0=130.6(3) "
        "means an uncertainty of 0.3 A^3, not 3 A^3.",
    ),
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
    "10.1016/j.rgg.2013.01.005": source(
        "https://doi.org/10.1016/j.rgg.2013.01.005",
        [
            "Equations 1-15",
            "Table 1 input reference values",
            "Table 4 optimized parameters for all eleven markers",
            "Figures 3-11 cross-calibration tests",
        ],
        "The stored Sokolova marker coefficients reproduce the corresponding "
        "columns of Tables 1 and 4. Sokolova et al. (2016) is recorded "
        "separately as the spreadsheet implementation and correction source, "
        "not represented as a new experimental fit.",
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
    "10.1016/0022-3697(93)90106-2": source(
        "https://doi.org/10.1016/0022-3697(93)90106-2",
        [
            "page 6, ambient NiAs-type lattice parameters",
            "page 6, Table 1, measured P-a-c data",
            "page 6, third-order Birch-Murnaghan fit and Figure 1",
        ],
        "The conventional hexagonal-cell V0 and its error are independently "
        "propagated from a=3.4395(2) A and c=5.3514(7) A. The paper reports "
        "no parameter covariance.",
    ),
    "10.1016/0031-9201(92)90063-2": source(
        "https://doi.org/10.1016/0031-9201(92)90063-2",
        [
            "page 4, Figure 3 and Birch-Murnaghan fit paragraph",
            "page 4, Table 1 and V0 footnote",
            "page 5, Table 2, this-study 1 bar cell",
        ],
        "The source prints K0=161.2 GPa for fixed K0'=4 and V0=1513.1 A^3, "
        "but no fit uncertainty for K0 or V0. The migrated +/-4 GPa was not "
        "a fit error from this EOS and is removed.",
    ),
    "10.1016/s0012-821x(00)00033-9": source(
        "https://doi.org/10.1016/S0012-821X(00)00033-9",
        [
            "page 71, Table 2, AntA-D and AntB-D P-V data",
            "pages 78-79, second-order Birch-Murnaghan fit and Figure 10",
        ],
        "The paper normalizes two antigorite-derived phase-D runs with "
        "different measured ambient volumes. Peritheos therefore exposes two "
        "explicit reference-volume variants sharing the published K0=134(5) "
        "GPa, rather than inventing one shared absolute V0.",
    ),
    "10.1016/s0038-1098(99)00322-1": source(
        "https://doi.org/10.1016/S0038-1098(99)00322-1",
        [
            "page 125, Table 1, bcc and fcc P-V observations",
            "page 125, Table 2, experimental Vinet parameters",
            "pages 125-126, combined bcc-fcc fit description and Figure 3",
        ],
        "This is explicitly one empirical Vinet curve fitted across the small "
        "bcc-fcc volume discontinuity. Its bcc two-atom-cell volume convention "
        "is retained, and it is not described as a phase-specific bcc EOS.",
    ),
    "10.1029/94jb00127": source(
        "https://doi.org/10.1029/94JB00127",
        [
            "page 11767, Table 1, separate RbCl and CsCl compression blocks",
            "page 11767, Figures 2-3 and Eulerian finite-strain fits",
        ],
        "The paper fixes the CsCl reference from a0=4.123 A and the RbCl-B2 "
        "reference from the hypothetical zero-pressure density 3.3068(10) "
        "Mg m^-3. Table 1 separately prints 13 CsCl and 24 RbCl observations; "
        "the CsCl fit also uses nine corrected V/V0 values from Yagi (1978).",
    ),
    "10.1029/jb078i035p08470": source(
        "https://doi.org/10.1029/JB078i035p08470",
        [
            "page 8471, Table 2, 23+/-3 degC P-V observations",
            "page 8472, first Birch formulation and fitted coefficients",
            "page 8472, stated one-standard-deviation fit errors",
        ],
        "The first printed Birch formulation is the BM3 form represented by "
        "the record. The fixed ambient V0 error uses the paper's stated 0.3% "
        "volume accuracy. The fit spans a small reversible tetragonal "
        "distortion without a volume discontinuity.",
    ),
    "10.1029/jb079i008p01165": source(
        "https://doi.org/10.1029/JB079i008p01165",
        [
            "page 1165, ambient a=8.394(3) A",
            "page 1167, Table 3 and printed Birch-Murnaghan equation",
            "pages 1167-1168, assumed K0'=4+/-0.4 and K0 error budget",
        ],
        "V0 and its error are propagated from the paper's own ambient lattice "
        "parameter. The +/-10 GPa K0 error is the authors' sum of independent "
        "fit, pressure-scale, and assumed-K0' contributions, not a covariance "
        "matrix or a uniform one-sigma statistical error.",
    ),
    "10.1029/jb086ib12p11773": source(
        "https://doi.org/10.1029/JB086iB12p11773",
        [
            "page 11774, Table 1, B2-SrO density observations",
            "pages 11775-11776, equations (1)-(7), finite-strain derivation",
            "page 11776, second-order fit and Figure 3",
        ],
        "The source reports rho0=6.14(20) Mg/m^3 and K0=160(19) GPa with "
        "K0'=4. Peritheos converts rho0 and its uncertainty to the one-formula-"
        "unit B2 conventional-cell volume; no covariance is published.",
    ),
    "10.1029/jb094ib03p03037": source(
        "https://doi.org/10.1029/JB094iB03p03037",
        [
            "page 3037, MW60 ambient a=4.2805(4) A",
            "pages 3038-3039, Tables 1-2, MW60 BM2 fit and P-V data",
            "page 3041, equation (1) and discussion of fixing K0'=4",
        ],
        "This record is the fit-specific MW60 second-order result in Table 1, "
        "not the paper's separate global composition model. V0 error is "
        "propagated from the primary ambient lattice measurement.",
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
    "10.1029/2000gl012606": source(
        "https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2000GL012606",
        [
            "page 1875, Experimental section, ambient lattice parameters and V0",
            "page 1876, weighted least-squares EOS paragraph",
            "Figure 3 and caption, page 1877",
        ],
        "The stored record is the reported 300 K third-order Birch-Murnaghan "
        "fit. V0 is the independently measured ambient reference volume; the "
        "paper reports no fit covariance or confidence convention for the "
        "printed +/- values. The abstract and title round the compression "
        "range to 73 GPa; Figure 1 labels the highest diffraction pattern "
        "73.2 GPa, which is retained as the experimental maximum.",
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
    "10.1063/1.1726610": source(
        "https://pubs.aip.org/aip/jcp/article/44/11/4223/209674/"
        "Lattice-Parameters-of-Nine-Oxides-and-Sulfides-as",
        [
            "Equation (4), journal page 4226",
            "Tables II, III, and VI, journal pages 4224-4227",
        ],
        "The CoO record is the empirical whole-range Murnaghan fit printed in "
        "Table VI, not a Birch-Murnaghan fit. V0 is calculated from the Table II "
        "ambient cubic lattice parameter. The paper prints no parameter errors "
        "or covariance and cautions that B0 and B0' need not equal the true "
        "one-atmosphere derivatives.",
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
        [
            "Table I, journal pages 104505-2 to 104505-3",
            "Equations (1)-(3), journal page 104505-4",
            "Figure 2 caption and Table II",
        ],
        "The source uses BM2 with T0=300 K, exponential constant-alpha "
        "reference-volume scaling, and dK/dT fixed to zero. Figure 2 displays "
        "the ice-VI P-V series at 300 and 340 K; fitting those 30 endpoint-isotherm "
        "rows reproduces Table II. Molar volumes were converted with Z=10 for ice "
        "VI and Z=2 for ice VII.",
    ),
    "10.1063/1.342969": source(
        "https://doi.org/10.1063/1.342969",
        [
            "ambient density paragraph, page 1535",
            "Section III and Equation (8), pages 1536-1537",
            "Equations (20)-(29), pages 1540-1541",
            "Table V, page 1541",
        ],
        "The stored record reproduces Equation (29) and the top rows of Table V: "
        "a 300 K BM3 isotherm plus a logarithmic-volume thermal-pressure slope. "
        "The partial published error on (dKT/dT)V is retained; the article also "
        "identifies an additional unquantified contribution from K0' uncertainty.",
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
        ["Equations 1-3", "Table 1", "Ne EOS discussion and Figure 5 caption"],
        "Equation 3 uses the published variable-exponent Debye-temperature law. "
        "The room-temperature Ne fit combines Hemley et al. (1989, ref. 45), "
        "Finger et al. (1981, ref. 47), and this study. Figure 5 states that the "
        "Hemley pressures were recalculated to the Dewaele et al. (2004) ruby scale.",
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
        "The experimental EOS is fitted to relative volume and reports only K0 "
        "and K0' as fit coefficients, so the measured ambient V0 is fixed in the "
        "reproduction. The measured hexagonal subcell was converted to the "
        "equivalent four-formula-unit orthorhombic conventional-cell volume used "
        "by the record.",
    ),
    "10.1098/rsta.2022.0331": source(
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC10069115/",
        ["Abstract", "Figure 1 caption", "Table 1"],
        "The paper is internally inconsistent: the abstract says third-order "
        "Birch-Murnaghan, while the Figure 1 caption says second-order despite "
        "reporting K0'=3.3(1). Peritheos retains BM3 and records the conflict.",
    ),
    "10.1103/fxgq-96sg": source(
        "https://doi.org/10.1103/fxgq-96sg",
        [
            "Equation (4), page 8",
            "Table I and footnotes, page 4",
            "Table II, pages 9-12",
            "Section III.E, pages 8 and 12-14",
            "Experimental methods, page 2",
        ],
        "The ten stored records are the Cu-anchored reduced-300 K Vinet P-V "
        "fits. V0 is fixed; the table's quoted K0 and K0' uncertainties are "
        "preserved. The article does not state a confidence level or publish "
        "parameter covariance, so neither is inferred.",
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
        "Table 3 row. Its footnote specifies a staged fit of K0 and K0' to the "
        "room-temperature data before alpha0; the italic row is the distinct "
        "simultaneous four-parameter fit. The article explicitly declines "
        "individual parameter errors because of covariance, but reports "
        "alpha0*K0=0.0275+/-0.0009 kbar/K.",
    ),
    "10.2138/am-2019-6779": source(
        "https://doi.org/10.2138/am-2019-6779",
        [
            "Equations 1-6",
            "Table 1",
            "pages 721-723 of the final article",
            "MSA deposit AM-19-56779, Supplemental Table S1 workbook",
        ],
        "The final article's preferred Sokolova-Pt fit fixes V0=54.5 A^3 and "
        "theta0=235 K and reports K0=18.3(3) GPa, K0'=5.60(3), gamma0=2.3(2), "
        "and q=0.8(2). Equation 6 is the integrated-Gruneisen Debye-temperature "
        "law. The accepted manuscript's gamma0=0.58(5) and q=0.9(2) are "
        "superseded.",
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
        [
            "Equation 1 and goethite EOS discussion, manuscript pages 6-7",
            "Goethite thermodynamic parameters, Table 1 on manuscript page 15",
            "MSA depository item AM-08-056, Table 1",
        ],
        "The source defines a BM3 reference plus Mie-Gruneisen-Debye thermal "
        "pressure and says that all 65 deposited P-V-T rows determine K0 and "
        "K0'. The deposit's row 32 values are retained verbatim but flagged as "
        "an isolated source-data/refinement anomaly. An independent fit of the "
        "47 plotted Figure 3a marker centers also fails badly, and Equation 1 "
        "is malformed as printed; the record is retained for provenance but is "
        "not recommended for quantitative use.",
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
    "neon_fcc_hemley_1989_bm3_refit": source(
        "https://doi.org/10.1103/PhysRevB.39.11820",
        [
            "Table I, page 11823",
            "Equation-of-state discussion and Figure 4, pages 11823-11824",
        ],
        "Table I provides 21 raw 300 K observations on the Mao ruby scale. "
        "The stored coefficients are a Peritheos errors-in-variables fit to "
        "that table alone; the paper's published thermally reduced 0 K fit "
        "also includes the earlier Finger et al. low-pressure series.",
    ),
    "diamond_correa_2008_dewaele_anchored": source(
        "https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.094106",
        [
            "Dewaele et al. (2008), Table III",
            "Correa et al. (2008), equations 2-7 and 13-18 and Table I",
        ],
        "Derived composition of the primary-source-validated Dewaele 298 K "
        "Vinet isotherm and the reference-relative Correa thermal free energy.",
    ),
    "diamond_benedict_2014_dewaele_anchored": source(
        "https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.094106",
        [
            "Dewaele et al. (2008), Table III",
            "Benedict et al. (2014), equations 3-7 and Table I",
        ],
        "Derived composition of the primary-source-validated Dewaele 298 K "
        "Vinet isotherm and the reference-relative Benedict thermal free energy.",
    ),
    "aragonite_martinez_1996_bm2_2": source(
        "https://rruff.info/doclib/am/vol81/AM81_611.pdf",
        [
            "Equation (1), page 615",
            "Table 3, page 616",
            "Table 6 and EOS discussion, pages 618-619",
            "Table 7, page 620",
        ],
        "The Table 6 isotherms are BM2 fits with K0'=4 assumed. Regressing "
        "their printed V0,T values gives a mean expansion coefficient of "
        "6.484e-5 K^-1, reproducing the Table 7 value 6.5(1)e-5 K^-1. "
        "Unweighted and error-weighted regressions of the printed K0,T values "
        "give -0.01969 and -0.01702 GPa/K, bracketing the reported "
        "-0.018(2) GPa/K. The fit-specific 298 K K0 error is retained; the "
        "conflicting summary error printed in the prose and Table 7 is "
        "recorded explicitly.",
    ),
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
    "kcl_b2_dewaele_2012_vinet_3": source(
        "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.85.214105/fulltext",
        [
            "Table I, 298 K measured P-V points",
            "Table III, experimental B2 KCl Vinet fit and 2.6-165 GPa range",
            "Equation (2) and Table V, P-V-T calibration",
            "Section IV, stated 0-200 GPa and 300-7000 K calibration range",
        ],
        "The paper provides the complete B2-KCl P-V-T equation and parameters, "
        "explicitly recommends it for pressure calibration, and prints no "
        "fitted-parameter errors or covariance.",
    ),
    "kcl_b2_chidester_2021_bm3_5": source(
        "https://link.aps.org/accepted/10.1103/PhysRevB.104.094107",
        [
            "Methods, effective KCl temperature and Pt pressure basis",
            "Section III and Equations 1-4",
            "Table I and Figure 3",
            "Supplemental Table A1; author-deposited SuppTable_KCl.csv",
            "Dewaele et al. (2012), Table I",
        ],
        "The preferred BM3+MGD fit reports V0=32.0(3) cm^3/mol, "
        "K0=24(1) GPa, K0'=4.56(5), gamma0=2.9(4), q=1.0(1), and fixed "
        "thetaD=235 K. It simultaneously fits all Dewaele et al. (2012) "
        "room-temperature B2 data and the 155 new high-temperature rows. "
        "Equations 3-4 leave theta(V) implicit; the thermodynamically integrated "
        "constant-q relation reproduces all five coefficients and the reported "
        "1.6 GPa high-temperature RMSE.",
    ),
    "kcl_campbell_1991_bm2_1": source(
        "https://doi.org/10.1016/0022-3697(91)90181-X",
        [
            "Campbell and Heinz (1991), primary abstract: B2/B1 V0 ratio, K0, fixed K0', and 2-56 GPa range",
            "Dewaele et al. (2012), Table III: experimental B1-KCl V0",
        ],
        "This is an explicit two-primary-source composite. Campbell and Heinz "
        "report V0(B2)/V0(B1)=0.8483+/-0.0057, K0=28.7+/-0.6 GPa, and "
        "K0'=4 fixed. Multiplication by Dewaele et al.'s B1 V0=62.36 "
        "A^3/formula unit gives V0(B2)=52.899988 A^3. Because Dewaele et al. "
        "publish no V0 error, only the Campbell ratio error is propagated.",
    ),
    "indium_nitride_munoz_1993_murnaghan_1": source(
        "https://iopscience.iop.org/article/10.1088/0953-8984/5/33/010/pdf",
        [
            "section 3, pages 6016-6018",
            "Figure 1 caption",
            "Table 1, wurtzite (PP) row",
        ],
        "The paper fits theoretical E(V) points with the Murnaghan equation. "
        "The conventional wurtzite-cell V0 is calculated from its theoretical "
        "a0=3.483 A and c0=5.7039 A. The finite-plane-wave-cutoff sensitivity "
        "is a convergence estimate, not a covariance-derived fit error.",
    ),
    "nickel_oxide_noguchi_1999_bm3_1": source(
        "https://www.jstage.jst.go.jp/article/jshpreview1992/7/0/7_0_832/_pdf",
        [
            "Noguchi et al. (1998), page 832, ambient lattice parameter",
            "Noguchi et al. (1998), pages 833-834, Mie-Gruneisen reduction and Murnaghan-Birch fit",
            "Noguchi et al. (1999), journal abstract, final 147.6 GPa isotherm coefficients",
        ],
        "The official open primary conference paper documents the same team's "
        "sample reference lattice and shock-to-300 K reduction. The final journal "
        "article extends the data and reports K0=191 GPa and K0'=3.9. The stored "
        "Z=3 rhombohedral-cell V0 and uncertainty are propagated from a0=4.177(1) A; "
        "the journal reports no coefficient errors for K0 or K0'.",
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


# Valid DOI groups may still contain a specifically deferred record. This map
# is keyed by record ID so DOI-level evidence cannot accidentally promote one.
FORCED_DEFERRED: dict[str, str] = {}


DEFERRED_BY_DOI: dict[str, str] = {
    "10.1016/0022-3697(91)90181-x": (
        "The primary abstract reports the B2/B1 zero-pressure volume ratio, K0, "
        "and fixed K0'=4, but the stored absolute B2 V0 additionally depends on an "
        "ambient B1 volume that is not established by the accessible primary text."
    ),
    "10.1029/rf002p0029": (
        "The cited DOI is a handbook chapter on thermal expansion rather than an "
        "unambiguous original source for the stored FeO static EOS parameter set. "
        "The record is deferred until its actual primary EOS source is identified."
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
}


# Equation (1) and Table 1 of Sokolova et al. require both the atom count n
# and atomic/effective atomic number Z.  The mechanical Dioptas migration did
# not carry all of these model inputs, so the audit restores them directly
# from the primary table rather than treating a software catalog as authority.
SOKOLOVA_COMPOSITION: dict[str, tuple[float, float]] = {
    "aluminum_sokolova_2013_holzapfel_2": (1.0, 13.0),
    "copper_sokolova_2013_holzapfel_2": (1.0, 29.0),
    "diamond_sokolova_2013_holzapfel_3": (1.0, 6.0),
    "gold_sokolova_2013_holzapfel_4": (1.0, 79.0),
    "mgo_sokolova_2013_holzapfel_4": (2.0, 10.34),
    "molybdenum_sokolova_2013_holzapfel_2": (1.0, 42.0),
    "niobium_sokolova_2013_holzapfel_2": (1.0, 41.0),
    "platinum_sokolova_2013_holzapfel_3": (1.0, 78.0),
    "silver_sokolova_2013_holzapfel_2": (1.0, 47.0),
    "tantalum_sokolova_2013_holzapfel_3": (1.0, 73.0),
    "tungsten_sokolova_2013_holzapfel_4": (1.0, 74.0),
}

SOKOLOVA_2013_REFERENCE = {
    "authors": ["Sokolova", "Dorogokupets", "Litasov"],
    "year": 2013,
    "source": "Russian Geology and Geophysics",
    "volume": "54",
    "locator": "181-199",
    "details": "Tables 1 and 4",
    "doi": "10.1016/j.rgg.2013.01.005",
}

SOKOLOVA_2016_IMPLEMENTATION = {
    "role": "spreadsheet implementation, conventions, and corrections",
    "citation": "Sokolova et al. (2016), equations 1-12, Tables 1-3, and Appendix A",
    "doi": "10.1016/j.cageo.2016.06.002",
}

SUN_SILICA_RECORDS = {
    "silica_cacl2_sun_2019_bm2_1",
    "seifertite_sun_2019_bm2_1",
}

BEZACIER_ICE_RECORDS = {
    "ice_vi_bezacier_2014_bm2_1",
    "ice_vii_bezacier_2014_bm2_1",
}

REMOVED_UNSUPPORTED_RECORDS = {
    "aragonite_martinez_1996_bm3_1",
    "feo_fei_1995_bm3_1",
    "tungsten_hixson_1992_bm3_1",
}


def curate_migrated_catalog() -> None:
    """Apply record removals, consolidation, and source-required record splits."""
    majorite_path = MATERIALS / "majorite.eosmat"
    duplicate_path = MATERIALS / "mgsio3.eosmat"
    source_path = duplicate_path if duplicate_path.exists() else majorite_path
    document = json.loads(source_path.read_text(encoding="utf-8"))
    document.update(
        {
            "name": "Majorite (MgSiO3 tetragonal garnet)",
            "aliases": ["MgSiO3 majorite", "tetragonal MgSiO3 garnet"],
            "formula": "MgSiO3",
            "formula_units_per_cell": 16,
            "space_group": "I41/a",
            "space_group_number": 88,
            "phase": "tetragonal MgSiO3 garnet (majorite), I41/a",
            "identifier": "majorite",
            "notes": (
                "Tetragonal MgSiO3 majorite consolidated from the duplicate "
                "Dioptas/JCPDS majorite and mgsio3-maj entries. The richer "
                "mgsio3-maj diffraction list is retained. Majorite has the "
                "garnet formula Mg3(MgSi)Si3O12, equivalent to 4 MgSiO3, and "
                "space group I41/a with 16 MgSiO3 formula units in the "
                "conventional cell. Source file(s): JCPDS/current user/jcpds/"
                "mgsio3-maj.jcpds, JCPDS/dac_user_jcpds/mgsio3-maj.jcpds, "
                "JCPDS/current user/jcpds/ver3/majorite.jcpds."
            ),
        }
    )
    record = document["eos_records"][0]
    record["identifier"] = "majorite_yagi_1992_bm3_1"
    migration = record.get("scientific_validation", {}).get("migration_source")
    if migration is not None:
        migration["file"] = "mgsio3.json; majorite.json"
    majorite_path.write_text(
        json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    if duplicate_path.exists():
        duplicate_path.unlink()

    phase_d_path = MATERIALS / "phase_d.eosmat"
    phase_d = json.loads(phase_d_path.read_text(encoding="utf-8"))
    migrated = phase_d["eos_records"][0]
    migration = migrated.get("scientific_validation", {}).get("migration_source")

    def phase_d_variant(
        run: str,
        identifier: str,
        v0: float,
        v0_error: float,
        maximum_pressure: float,
    ) -> dict[str, Any]:
        return {
            "label": f"Shieh et al. (2000), antigorite-derived phase D ({run})",
            "reference": dict(migrated["reference"]),
            "eos": {
                "type": "BM2",
                "parameters": {"V0": v0, "K0": 134.0},
                "model": "birch_murnaghan_2",
            },
            "parameter_errors": {"V0": v0_error, "K0": 5.0},
            "fixed_parameters": ["V0"],
            "temperature_ref": 300.0,
            "experimental_pressure_range_gpa": [0.0, maximum_pressure],
            "pressure_range_status": "reported_exactly",
            "parameter_provenance": {
                "V0": (
                    f"Shieh et al. (2000), Table 2, {run}-D ambient unit-cell volume"
                ),
                "K0": (
                    "Shieh et al. (2000), pages 78-79: joint fit to the "
                    "antigorite-derived phase-D volume data"
                ),
                "equation": (
                    "Shieh et al. (2000), pages 78-79: second-order "
                    "Birch-Murnaghan (K0'=4)"
                ),
            },
            "notes": (
                f"Reference-volume variant for the {run} antigorite-derived "
                f"phase-D run. Table 2 reports V0={v0}+/-{v0_error} A^3. "
                "The paper fits the combined antigorite-derived data with "
                "a second-order Birch-Murnaghan curve and reports K0=134+/-5 "
                "GPa, but it does not define one shared absolute V0 because AntA "
                "and AntB have different ambient volumes. The two Peritheos "
                "records preserve those measured reference volumes explicitly. "
                "No parameter covariance or confidence convention is reported."
            ),
            "identifier": identifier,
            "scientific_validation": {
                "status": "pending_primary_source_check",
                **({"migration_source": dict(migration)} if migration else {}),
            },
        }

    phase_d["eos_records"] = [
        phase_d_variant(
            "AntA",
            "phase_d_ant_a_shieh_2000_bm2_1",
            88.12,
            0.32,
            24.6,
        ),
        phase_d_variant(
            "AntB",
            "phase_d_ant_b_shieh_2000_bm2_1",
            87.191,
            0.097,
            22.0,
        ),
    ]
    phase_d["phase"] = "Fe-bearing hexagonal phase D derived from antigorite"
    phase_d["notes"] = (
        "Hexagonal phase-D diffraction interchange record. The crystallographic "
        "pattern was imported from Dioptas/JCPDS, while both executable EOS "
        "records use the separately printed AntA and AntB ambient volumes in "
        "Shieh et al. (2000), Table 2. Phase D is compositionally variable; the "
        "paper's antigorite starting material contained about 3.7 wt% Fe."
    )
    phase_d_path.write_text(
        json.dumps(phase_d, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


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
    if identifier == "kcl_campbell_1991_bm2_1":
        old_v0 = record["eos"]["parameters"]["V0"]
        old_v0_error = record["parameter_errors"].get("V0")
        record["label"] = (
            "Campbell and Heinz (1991) + Dewaele et al. (2012), B2 BM2 composite"
        )
        record["eos"]["parameters"]["V0"] = 52.899988
        record["parameter_errors"]["V0"] = 0.355452
        record["fixed_parameters"] = ["V0"]
        record["temperature_ref"] = 298.0
        record["parameter_provenance"] = {
            "V0": (
                "Campbell and Heinz (1991) V0(B2)/V0(B1)=0.8483+/-0.0057 "
                "multiplied by Dewaele et al. (2012) experimental B1-KCl "
                "V0=62.36 A^3/formula unit"
            ),
            "K0": "Campbell and Heinz (1991) primary abstract: 28.7+/-0.6 GPa",
            "K0_prime": "Campbell and Heinz (1991) primary abstract: fixed at 4.0 (BM2)",
        }
        record["notes"] = (
            "Explicit composite B2-KCl reference isotherm. Campbell and Heinz "
            "report V0(B2)/V0(B1)=0.8483+/-0.0057, K0=28.7+/-0.6 GPa, "
            "and K0'=4 constrained over 2-56 GPa. The absolute B2 V0 is "
            "0.8483*62.36=52.899988 A^3/formula unit, using the independently "
            "published experimental B1-KCl V0 of Dewaele et al. (2012). "
            "Its 0.355452 A^3 uncertainty propagates the Campbell ratio error "
            "only; Dewaele et al. report no B1 V0 error or covariance. This "
            "composite is not presented as a single-paper Campbell parameter set."
        )
        for correction in (
            {
                "path": "eos.parameters.V0",
                "source_value": old_v0,
                "value": 52.899988,
                "reason": "Derive the absolute B2 reference volume from two explicitly named primary inputs.",
            },
            {
                "path": "parameter_errors.V0",
                "source_value": old_v0_error,
                "value": 0.355452,
                "reason": "Propagate 62.36*0.0057; no B1-V0 error is published to include.",
            },
            {
                "path": "fixed_parameters",
                "source_value": [],
                "value": ["V0"],
                "reason": "The absolute V0 is constructed from published inputs rather than fitted in this record.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1016/0022-3697(91)90181-X + 10.1103/PhysRevB.85.214105",
                "location": "Campbell abstract; Dewaele Table III",
            }
            append_correction(record, correction)
    if identifier == "cscl_campbell_1994_bm3_1":
        record["label"] = "Campbell and Heinz (1994), CsCl BM3"
        record["eos"] = {
            "type": "BM3",
            "parameters": {
                "V0": 70.087408867,
                "K0": 17.01,
                "K0_prime": 5.49,
            },
            "model": "birch_murnaghan_3",
        }
        record["parameter_errors"] = {
            "V0": None,
            "K0": 0.29,
            "K0_prime": 0.15,
        }
        record["fixed_parameters"] = ["V0"]
        record["temperature_ref"] = 298.0
        record["experimental_pressure_range_gpa"] = [0.0, 28.7]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": (
                "Campbell and Heinz (1994), page 11767: accepted ambient "
                "a0=4.123 A; V0=a0^3 for Pm-3m CsCl"
            ),
            "K0": "Campbell and Heinz (1994), page 11767 and Figure 2",
            "K0_prime": "Campbell and Heinz (1994), page 11767 and Figure 2",
            "equation": (
                "Campbell and Heinz (1994), page 11767: Eulerian finite-strain "
                "analysis with linear normalized stress"
            ),
        }
        record["notes"] = (
            "Room-temperature third-order Birch-Murnaghan fit to the authors' "
            "13 Table 1 observations and nine corrected Yagi (1978) Table 1 "
            "compression ratios. The primary paper explicitly "
            "uses a0=4.123 A, giving V0=70.087408867 A^3 for the one-formula-"
            "unit Pm-3m cell; it reports no uncertainty for that fixed lattice "
            "parameter. K0=17.01(29) GPa and K0'=5.49(15). Campbell and Heinz do "
            "not state the numerical weights used for their normalized-stress "
            "regression, so exact central-value reproduction remains fit-protocol limited."
        )
    if identifier == "fe3o4_mao_1974_bm3_1":
        record["label"] = "Mao et al. (1974), magnetite BM3"
        record["eos"] = {
            "type": "BM3",
            "parameters": {
                "V0": 591.434826984,
                "K0": 183.0,
                "K0_prime": 4.0,
            },
            "model": "birch_murnaghan_3",
        }
        record["parameter_errors"] = {
            "V0": 0.634133124,
            "K0": 10.0,
            "K0_prime": 0.4,
        }
        record["fixed_parameters"] = ["V0", "K0_prime"]
        record["temperature_ref"] = 296.15
        record["experimental_pressure_range_gpa"] = [0.0, 32.0]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": (
                "Mao et al. (1974), page 1165: a0=8.394+/-0.003 A; "
                "V0=a0^3 and sigma(V0)=3*a0^2*sigma(a0)"
            ),
            "K0": "Mao et al. (1974), pages 1167-1168, summed error budget",
            "K0_prime": ("Mao et al. (1974), pages 1167-1168: assumed 4+/-0.4"),
            "equation": "Mao et al. (1974), page 1167, printed Birch-Murnaghan equation",
        }
        record["notes"] = (
            "Room-temperature (23+/-3 degC) magnetite BM3. V0 and its "
            "propagated error come from the paper's a0=8.394(3) A sample, not "
            "the independent structural lattice in this file. K0'=4+/-0.4 was "
            "assumed, and the authors formed the +/-10 GPa K0 uncertainty by "
            "adding fit, NaCl pressure-scale, and assumed-K0' contributions; "
            "it is not a covariance-derived one-sigma error. The low-pressure "
            "spinel reflections persist to 32 GPa, but a high-pressure phase "
            "begins to coexist above about 25 GPa."
        )
    if identifier == "majorite_yagi_1992_bm3_1":
        record["label"] = "Yagi et al. (1992), MgSiO3 majorite BM3"
        record["eos"] = {
            "type": "BM3",
            "parameters": {"V0": 1513.1, "K0": 161.2, "K0_prime": 4.0},
            "model": "birch_murnaghan_3",
        }
        record["parameter_errors"] = {
            "V0": None,
            "K0": None,
            "K0_prime": None,
        }
        record["fixed_parameters"] = ["V0", "K0_prime"]
        record["temperature_ref"] = 298.0
        record["experimental_pressure_range_gpa"] = [0.0, 9.72]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": "Yagi et al. (1992), Table 1 footnote and Table 2 this-study row",
            "K0": "Yagi et al. (1992), page 4, Figure 3 fit paragraph",
            "K0_prime": "Yagi et al. (1992), page 4: fixed at 4",
            "equation": "Yagi et al. (1992), page 3, Birch-Murnaghan fit description",
        }
        record["notes"] = (
            "Room-temperature MgSiO3 tetragonal-garnet BM3 with the measured "
            "one-bar V0=1513.1 A^3 and K0'=4 assumed. The fit gives K0=161.2 "
            "GPa. The paper prints no fit error or covariance for V0 or K0, so "
            "the migrated +/-4 GPa value is removed rather than inferred from "
            "the separate composition-regression uncertainties. The diffraction "
            "structure in this file is an independent 1518.5 A^3 cell."
        )
    if identifier == "mgfe60o_richet_1989_bm3_1":
        record["label"] = "Richet et al. (1989), (Mg0.4Fe0.6)O BM2"
        record["identifier"] = "mgfe60o_richet_1989_bm2_1"
        record["eos"] = {
            "type": "BM2",
            "parameters": {"V0": 78.430232810125, "K0": 149.0},
            "model": "birch_murnaghan_2",
        }
        record["parameter_errors"] = {"V0": 0.0219872163, "K0": 4.0}
        record["fixed_parameters"] = ["V0"]
        record["temperature_ref"] = 298.0
        record["experimental_pressure_range_gpa"] = [0.0, 49.4]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": (
                "Richet et al. (1989), page 3037: a0=4.2805(4) A; "
                "V0=a0^3 for the four-formula-unit rocksalt cell"
            ),
            "K0": "Richet et al. (1989), Table 1, MW60 this-work static fit",
            "equation": (
                "Richet et al. (1989), Table 1 footnote and equation (1): "
                "second-order Eulerian finite strain with K0'=4"
            ),
        }
        record["notes"] = (
            "Fit-specific room-temperature BM2 for MW60, the sample containing "
            "60 mol% FeO. It is distinct from the paper's later global linear-"
            "composition model. V0 and its error are propagated from the primary "
            "a0=4.2805(4) A measurement; K0=149+/-4 GPa and K0'=4 fixed. "
            "Table 2 spans 0-49.4 GPa. No covariance or confidence convention "
            "for the K0 error is reported."
        )
    if identifier == "nis_campbell_1993_bm3_1":
        record["label"] = "Campbell and Heinz (1993), NiAs-type NiS BM3"
        record["eos"] = {
            "type": "BM3",
            "parameters": {
                "V0": 54.8262666013039,
                "K0": 156.0,
                "K0_prime": 4.4,
            },
            "model": "birch_murnaghan_3",
        }
        record["parameter_errors"] = {
            "V0": 0.009596193736601715,
            "K0": 10.0,
            "K0_prime": 1.2,
        }
        record["fixed_parameters"] = ["V0"]
        record["temperature_ref"] = 298.0
        record["experimental_pressure_range_gpa"] = [0.0, 44.9]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": (
                "Campbell and Heinz (1993), page 6: a0=3.4395(2) A and "
                "c0=5.3514(7) A; V0=sqrt(3)/2*a0^2*c0"
            ),
            "K0": "Campbell and Heinz (1993), page 6, BM3 fit paragraph",
            "K0_prime": "Campbell and Heinz (1993), page 6, BM3 fit paragraph",
            "equation": "Campbell and Heinz (1993), page 6 and Figure 1",
        }
        record["notes"] = (
            "Room-temperature third-order Birch-Murnaghan fit for metastable "
            "NiAs-type (B8) NiS. V0 and its independent propagated error use "
            "the primary paper's ambient a and c values. K0=156(10) GPa and "
            "K0'=4.4(1.2). The last tabulated readable P-V point is 44.9 GPa; "
            "above 45 GPa the diffraction lines become unreadable. No fitted-"
            "parameter covariance is published."
        )
    if identifier in {
        "indium_nitride_mu_oz_1993_bm3_1",
        "indium_nitride_munoz_1993_murnaghan_1",
    }:
        record["identifier"] = "indium_nitride_munoz_1993_murnaghan_1"
        record["label"] = "Muñoz and Kunc (1993), theoretical wurtzite Murnaghan"
        record["eos"] = {
            "type": "Murnaghan",
            "parameters": {
                "V0": 59.92519880888224,
                "K0": 166.0,
                "K0_prime": 3.8,
            },
            "model": "murnaghan",
        }
        record["parameter_errors"] = {"V0": None, "K0": None, "K0_prime": None}
        record["fixed_parameters"] = ["V0"]
        record["pressure_range_status"] = "theoretical"
        record.pop("experimental_pressure_range_gpa", None)
        record["parameter_provenance"] = {
            "V0": "Table 1 theoretical a0=3.483 A and c0=5.7039 A; V0=sqrt(3)/2*a0^2*c0",
            "K0": "Table 1 wurtzite pseudopotential row: 166 GPa",
            "K0_prime": "Table 1 wurtzite pseudopotential row: 3.8",
            "equation": "section 3 and Figure 1 caption: Murnaghan E(V)",
        }
        record["notes"] = (
            "Ab initio 0 K static Murnaghan energy-volume fit for wurtzite InN. "
            "The theoretical conventional-cell V0=59.92519880888224 A^3 is "
            "calculated from Table 1 a0=3.483 A and c0=5.7039 A; K0=166 GPa "
            "and K0'=3.8 are from the same pseudopotential row. These are not "
            "experimental EOS parameters. The paper estimates finite-cutoff "
            "sensitivity from 40/70-Ryd calculations, but does not publish "
            "statistical coefficient errors or covariance, so none are inferred."
        )
        for correction in (
            {
                "path": "identifier",
                "source_value": "indium_nitride_mu_oz_1993_bm3_1",
                "value": "indium_nitride_munoz_1993_murnaghan_1",
                "reason": "Correct the author spelling and identify the primary Murnaghan model.",
            },
            {
                "path": "eos",
                "source_value": "BM3(V0=61.7988, K0=166, K0_prime=3.8)",
                "value": "Murnaghan(V0=59.92519880888224, K0=166, K0_prime=3.8)",
                "reason": "Use the paper's Murnaghan fit and theoretical Table 1 equilibrium cell, not an experimental card V0.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1088/0953-8984/5/33/010",
                "location": "section 3, Figure 1 caption, and Table 1",
            }
            append_correction(record, correction)
    if identifier in {
        "li_bcc_hanfland_1999_bm3_1",
        "li_hanfland_1999_vinet_1",
    }:
        record["identifier"] = "li_hanfland_1999_vinet_1"
        record["label"] = "Hanfland et al. (1999), combined bcc-fcc Li Vinet"
        record["eos"] = {
            "type": "Vinet",
            "parameters": {"V0": 43.245, "K0": 11.32, "K0_prime": 3.62},
            "model": "vinet",
        }
        record["parameter_errors"] = {
            "V0": None,
            "K0": 0.10,
            "K0_prime": 0.04,
        }
        record["fixed_parameters"] = ["V0"]
        record["temperature_ref"] = 298.0
        record["experimental_pressure_range_gpa"] = [0.0, 21.1]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": (
                "Hanfland et al. (1999), Table 2: fixed atomic V0=21.6225 "
                "A^3, multiplied by two for the bcc conventional-cell convention"
            ),
            "K0": "Hanfland et al. (1999), Table 2, this-work Vinet row",
            "K0_prime": "Hanfland et al. (1999), Table 2, this-work Vinet row",
            "equation": (
                "Hanfland et al. (1999), page 125, combined bcc-fcc Vinet "
                "fit description and Figure 3"
            ),
        }
        record["notes"] = (
            "Empirical 298 K Vinet fit jointly spanning bcc and fcc lithium, "
            "exactly as defined by the primary paper. It is not a phase-specific "
            "bcc EOS: the bcc-fcc transition occurs near 7.5 GPa with a "
            "0.16(3)% volume discontinuity. Public volumes use twice the paper's "
            "atomic volume so they remain compatible with this file's two-atom "
            "bcc conventional-cell convention, including for fcc observations. "
            "V0 was fixed; no V0 error or parameter covariance is published."
        )
        for correction in (
            {
                "path": "identifier",
                "source_value": "li_bcc_hanfland_1999_bm3_1",
                "value": "li_hanfland_1999_vinet_1",
                "reason": "The cited study's combined bcc/fcc fit is Vinet, not a bcc BM3.",
            },
            {
                "path": "eos.type",
                "source_value": "BM3",
                "value": "Vinet",
                "reason": "Correct the migrated equation family to the primary Vinet fit.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1016/S0038-1098(99)00322-1",
                "location": "page 125, Tables 1-2 and Figure 3",
            }
            append_correction(record, correction)
    if identifier in {
        "sno2_cubic_27gpa_ono_2000_bm3_1",
        "sno2_pa_3_at_48gpa_ono_2000_bm3_1",
    }:
        record["label"] = "Ono et al. (2000), cubic Pa-3 SnO2 300 K BM3"
        record["eos"] = {
            "type": "BM3",
            "parameters": {"V0": 130.6, "K0": 252.0, "K0_prime": 3.5},
            "model": "birch_murnaghan_3",
        }
        record["parameter_errors"] = {
            "V0": 0.3,
            "K0": 28.0,
            "K0_prime": 2.2,
        }
        record["fixed_parameters"] = []
        record["temperature_ref"] = 300.0
        record["experimental_pressure_range_gpa"] = [16.09, 28.85]
        record["pressure_range_status"] = "reported_exactly"
        record["experimental_temperature_range_k"] = [300.0, 1400.0]
        record["validity"] = {
            "pressure_gpa": [16.09, 28.85],
            "temperature_k": [300.0, 300.0],
            "notes": [
                "This executable record is only the published 300 K reference isotherm.",
                "The paper also reports a separate 25 GPa thermal-expansion coefficient, which is not a complete thermal EOS and is not composed here.",
            ],
        }
        record["parameter_provenance"] = {
            "V0": "Ono et al. (2000), page 621, 300 K BM3 paragraph",
            "K0": "Ono et al. (2000), page 621, 300 K BM3 paragraph",
            "K0_prime": "Ono et al. (2000), page 621, 300 K BM3 paragraph",
            "equation": "Ono et al. (2000), page 621, printed third-order BM equation",
        }
        record["notes"] = (
            "The 300 K reference isotherm for cubic Pa-3 SnO2. V0=130.6(3) "
            "A^3 means 130.6+/-0.3 A^3; the migrated +/-3.0 value was a "
            "decimal-place error. K0=252(28) GPa and K0'=3.5(2.2). V0 is a "
            "fitted extrapolation because the cubic phase transforms on "
            "decompression. Table 2 spans 16.09-28.85 GPa and 300-1400 K, "
            "but this record intentionally implements only the 300 K curve. "
            "No fitted-parameter covariance is published."
        )
    if identifier == "sro_liu_1973_bm3_1":
        record["label"] = "Liu and Bassett (1973), SrO volumetric BM3"
        record["eos"] = {
            "type": "BM3",
            "parameters": {"V0": 137.388096, "K0": 91.3, "K0_prime": 4.3},
            "model": "birch_murnaghan_3",
        }
        record["parameter_errors"] = {
            "V0": 0.412164288,
            "K0": 2.7,
            "K0_prime": 0.3,
        }
        record["fixed_parameters"] = ["V0"]
        record["temperature_ref"] = 296.15
        record["experimental_pressure_range_gpa"] = [0.0, 34.05]
        record["experimental_temperature_range_k"] = [293.15, 299.15]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": (
                "Liu and Bassett (1973), Table 2: a0=5.160 A; V0=a0^3, "
                "with the stated 0.3% volume accuracy"
            ),
            "K0": "Liu and Bassett (1973), page 8472, first Birch formulation",
            "K0_prime": "Liu and Bassett (1973), page 8472, first Birch formulation",
            "equation": "Liu and Bassett (1973), page 8472, first printed Birch equation",
        }
        record["notes"] = (
            "Volumetric BM3 fit to all 23+/-3 degC SrO data through 34.05 GPa. "
            "K0=91.3+/-2.7 GPa and K0'=4.3+/-0.3 are one-standard-deviation "
            "fit errors. The fixed V0 error represents the paper's stated 0.3% "
            "volume accuracy. Between about 7 and 30.7 GPa the nominal B1 cell "
            "shows a small tetragonal distortion, but the paper reports no "
            "volume discontinuity and deliberately fits all data together. No "
            "parameter covariance is published."
        )
    if identifier in {
        "sro_b2_sato_1981_bm3_1",
        "sro_b2_sato_1981_bm2_1",
    }:
        record["identifier"] = "sro_b2_sato_1981_bm2_1"
        record["label"] = "Sato and Jeanloz (1981), B2 SrO BM2"
        record["eos"] = {
            "type": "BM2",
            "parameters": {"V0": 28.0224, "K0": 160.0},
            "model": "birch_murnaghan_2",
        }
        record["parameter_errors"] = {"V0": 0.9127817589576549, "K0": 19.0}
        record["fixed_parameters"] = []
        record["temperature_ref"] = 300.0
        record["experimental_pressure_range_gpa"] = [32.0, 59.0]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": (
                "Sato and Jeanloz (1981), page 11776: rho0=6.14+/-0.20 "
                "Mg/m^3, converted to the Z=1 Pm-3m conventional-cell volume"
            ),
            "K0": "Sato and Jeanloz (1981), page 11776, second-order fit",
            "equation": (
                "Sato and Jeanloz (1981), equations (1)-(7) and page 11776: "
                "second-order Eulerian finite strain with K0'=4"
            ),
        }
        record["notes"] = (
            "Room-temperature second-order finite-strain EOS for B2 SrO. The "
            "published rho0=6.14+/-0.20 Mg/m^3 converts to V0=28.0224+/-"
            "0.9128 A^3 for one SrO per Pm-3m cell; K0=160+/-19 GPa and "
            "K0'=4 fixed. B2 reflections were measured from 32 to 59 GPa, "
            "with B1-B2 coexistence at 32 and 37 GPa and a transition pressure "
            "of 36+/-4 GPa. No covariance is published."
        )
    if record.get("thermal", {}).get("type") == "AlphaKT":
        record["thermal"]["model"] = "thermal_reference_state"
    if identifier in {
        "coo_clendenen_1966_bm3_1",
        "coo_clendenen_1966_murnaghan_1",
    }:
        record["identifier"] = "coo_clendenen_1966_murnaghan_1"
        record["label"] = "Clendenen and Drickamer (1966), CoO Murnaghan"
        record["eos"] = {
            "type": "Murnaghan",
            "parameters": {
                "V0": 77.199941512,
                "K0": 190.5,
                "K0_prime": 3.9,
            },
            "model": "murnaghan",
        }
        record["parameter_errors"] = {
            "V0": None,
            "K0": None,
            "K0_prime": None,
        }
        record["fixed_parameters"] = ["V0"]
        record["experimental_pressure_range_gpa"] = [0.0, 30.8]
        record["pressure_range_status"] = "reported_exactly"
        record["parameter_provenance"] = {
            "V0": (
                "Table II ambient cubic a0=4.258 A; V0=a0^3 for the "
                "four-formula-unit conventional cell"
            ),
            "K0": "Table VI B0=1905 kbar; converted to 190.5 GPa",
            "K0_prime": "Table VI dimensionless B0'=3.9",
        }
        record["notes"] = (
            "Room-temperature empirical Murnaghan fit from Equation (4) and "
            "Table VI. The authors fitted B0 and B0' over the entire measured "
            "range and explicitly caution that they need not equal the true "
            "one-atmosphere bulk modulus and derivative. Table III gives the "
            "smoothed CoO compression data through 308 kbar (30.8 GPa). No "
            "parameter errors, covariance, confidence convention, or exact "
            "measurement temperature are printed."
        )
        for correction in (
            {
                "path": "identifier",
                "source_value": "coo_clendenen_1966_bm3_1",
                "value": "coo_clendenen_1966_murnaghan_1",
                "reason": "The primary paper uses the Murnaghan equation printed as Equation (4).",
            },
            {
                "path": "eos",
                "source_value": "BM3(V0=93.8242, K0=190.5, K0_prime=3.9)",
                "value": "Murnaghan(V0=77.199941512, K0=190.5, K0_prime=3.9)",
                "reason": (
                    "Use Equation (4), Table VI, and V0=a0^3 from the "
                    "Table II ambient a0=4.258 A."
                ),
            },
            {
                "path": "fixed_parameters",
                "source_value": [],
                "value": ["V0"],
                "reason": (
                    "Table VI fits B0 and B0'; V0 is the independently listed "
                    "ambient lattice reference in Table II."
                ),
            },
            {
                "path": "experimental_pressure_range_gpa",
                "source_value": None,
                "value": [0.0, 30.8],
                "reason": "Table III lists CoO compression data from 0 to 308 kbar.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1063/1.1726610",
                "location": "Equation (4) and Tables II, III, and VI",
            }
            append_correction(record, correction)
    if identifier == "aragonite_martinez_1996_bm2_2":
        record["label"] = "Martinez et al. (1996), staged P-V-T BM2"
        record["experimental_pressure_range_gpa"] = [0.0, 8.18]
        record["experimental_temperature_range_k"] = [298.0, 973.0]
        record["thermal"] = {
            "type": "AlphaKT",
            "model": "thermal_reference_state",
            "thermal_expansion_law": "constant",
            "reference_volume_law": "linear_temperature",
            "parameters": {
                "Tr": 298.0,
                "alpha0": 6.5e-5,
                "dK_dT": -0.018,
            },
            "parameter_errors": {
                "Tr": None,
                "alpha0": 0.1e-5,
                "dK_dT": 0.002,
            },
            "fixed_parameters": ["Tr"],
        }
        record["parameter_error_confidence"] = None
        record["notes"] = (
            "Staged second-order Birch-Murnaghan P-V-T parameterization from "
            "Martinez et al. (1996). Each Table 6 isotherm fixes K0'=4; Equation "
            "(2) supplies the reported linear K0(T) slope, and Equation (3) uses "
            "the reported mean expansion coefficient in the direct relation "
            "V0(T)=V0(298 K)*[1+alpha_bar*(T-298 K)]. Table 6 reports the 298 K "
            "V0=227.5+/-0.8 A^3 and sigma(K0)=3.48 GPa; the prose and Table 7 "
            "instead summarize the K0 error as approximately 4.3 GPa. Table 6 "
            "contains fitted isotherms through 973 K. The source reports no "
            "parameter covariance or confidence convention, and the stored "
            "pressure/temperature extrema are an experimental envelope rather "
            "than a rectangular validity guarantee."
        )
    if identifier == "cementite_scott_2001_bm3_1":
        old_fixed_parameters = list(record.get("fixed_parameters", []))
        record["eos"] = {
            "type": "BM3",
            "parameters": {"V0": 155.26, "K0": 175.4, "K0_prime": 5.1},
            "model": "birch_murnaghan_3",
        }
        record["parameter_errors"] = {
            "V0": 0.14,
            "K0": 3.5,
            "K0_prime": 0.3,
        }
        record["parameter_error_confidence"] = None
        record["fixed_parameters"] = ["V0"]
        record["experimental_pressure_range_gpa"] = [0.0, 73.2]
        record["experimental_temperature_range_k"] = [300.0, 300.0]
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Scott et al. report a weighted third-order Birch-Murnaghan fit to "
            "300 K compression data: V0=155.26+/-0.14 A^3 is the separately "
            "measured ambient unit-cell volume, K0T=175.4+/-3.5 GPa, and "
            "K0T'=5.1+/-0.3. The abstract and title round the compression "
            "range to 73 GPa; Figure 1 labels the highest diffraction pattern "
            "73.2 GPa, which is retained as the experimental maximum. The "
            "paper reports no parameter covariance or confidence convention "
            "for the printed +/- values. The 73 GPa envelope is not a thermal "
            "validity claim and should not be extrapolated to core pressures."
        )
        if old_fixed_parameters != ["V0"]:
            append_correction(
                record,
                {
                    "path": "fixed_parameters",
                    "source_value": old_fixed_parameters,
                    "value": ["V0"],
                    "reason": (
                        "The paper identifies V0 as the measured ambient volume "
                        "and reports only K0T and K0T' as EOS-fit results."
                    ),
                    "primary_reference": {
                        "doi": "10.1029/2000GL012606",
                        "location": "pages 1875-1876, Experimental and EOS paragraphs",
                    },
                },
            )
    if identifier in SOKOLOVA_COMPOSITION:
        n, atomic_number = SOKOLOVA_COMPOSITION[identifier]
        is_mgo = identifier == "mgo_sokolova_2013_holzapfel_4"
        record["label"] = (
            "Sokolova et al. (2013), Holzapfel thermal model (2016 workbook)"
        )
        record["reference"] = dict(SOKOLOVA_2013_REFERENCE)
        earlier_source = (
            {
                "role": "prior MgO thermodynamic EOS lineage",
                "citation": "Dorogokupets (2010)",
                "doi": "10.1007/s00269-010-0367-2",
            }
            if is_mgo
            else {
                "role": "preceding simultaneous-optimization fit",
                "citation": "Dorogokupets et al. (2012), Table 4",
                "doi": "10.5800/GT-2012-3-2-0067",
            }
        )
        record["source_lineage"] = [
            {
                "role": "reference volume, composition, and shock inputs",
                "citation": "Sokolova et al. (2013), Table 1",
                "doi": "10.1016/j.rgg.2013.01.005",
            },
            {
                "role": "final cross-calibrated Holzapfel and thermal coefficients",
                "citation": "Sokolova et al. (2013), Table 4",
                "doi": "10.1016/j.rgg.2013.01.005",
            },
            earlier_source,
            dict(SOKOLOVA_2016_IMPLEMENTATION),
        ]
        if is_mgo:
            record["source_lineage"].append(
                {
                    "role": "implemented MgO anharmonic-coefficient correction",
                    "citation": "Sokolova et al. (2016), Table 1 and correction discussion",
                    "doi": "10.1016/j.cageo.2016.06.002",
                }
            )
        if identifier == "diamond_sokolova_2013_holzapfel_3":
            record["source_lineage"].extend(
                [
                    {
                        "role": "principal room-temperature compression constraint",
                        "citation": "Occelli, Loubeyre, and LeToullec (2003)",
                        "doi": "10.1038/nmat831",
                    },
                    {
                        "role": "moderate-temperature diamond P-V-T context",
                        "citation": "Dewaele et al. (2008)",
                        "doi": "10.1103/PhysRevB.77.094106",
                    },
                    {
                        "role": "ambient heat-capacity constraint",
                        "citation": "Victor (1962)",
                        "doi": "10.1063/1.1701288",
                    },
                    {
                        "role": "ambient thermophysical-property constraint",
                        "citation": "Reeber and Wang (1996)",
                        "doi": "10.1007/BF02666175",
                    },
                    {
                        "role": "elastic-modulus constraint",
                        "citation": "Zouboulis et al. (1998)",
                        "doi": "10.1103/PhysRevB.57.2889",
                    },
                ]
            )
        record["notes"] = (
            "This record uses the self-consistent Sokolova et al. (2013) "
            "pressure-scale fit: Table 1 supplies reference volume/composition "
            "inputs and Table 4 supplies the final Holzapfel and thermal "
            "coefficients. The fit jointly considers shock-wave, ultrasonic, "
            "X-ray diffraction, dilatometric, and thermochemical measurements; "
            "it is not derived from one experimental dataset. Sokolova et al. "
            "(2016) republishes the parameters and supplies the executable "
            "Excel/VBA implementation, reference-temperature convention, and "
            "corrected equations, rather than a new fit dataset. "
            + (
                "The implemented MgO a0 coefficient follows the explicit 2016 "
                "correction, so that term has dual 2013/2016 provenance. "
                if is_mgo
                else ""
            )
            + "The publications do not provide individual parameter errors, "
            "parameter covariance, or a complete machine-readable list of fit "
            "points and weights. The *_sokolova_2013 identifier names the "
            "scientific fit year; the 2016 workbook remains explicit in "
            "source_lineage as the implementation and correction source."
        )
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
                        "restored from the Sokolova et al. (2016) model "
                        "definitions."
                    ),
                    "primary_reference": {
                        "doi": "10.1016/j.cageo.2016.06.002",
                        "location": "Equations (1)-(3) and Table 2",
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

    if identifier == "gold_anderson_1989_bm3_1":
        old_v0 = record["eos"]["parameters"]["V0"]
        old_thermal = dict(record["thermal"])
        record["label"] = "Anderson et al. (1989), Au logarithmic-volume thermal EOS"
        record["eos"]["parameters"]["V0"] = 67.79
        record["fixed_parameters"] = ["V0", "K0", "K0_prime"]
        record["temperature_ref"] = 300.0
        record["thermal"] = {
            "type": "LogVolumeThermalPressure",
            "parameters": {
                "Tr": 300.0,
                "alpha_KT_ref": 0.00714,
                "dK_dT_V": -0.0115,
            },
            "model": "log_volume_thermal_pressure",
            "parameter_errors": {
                "Tr": None,
                "alpha_KT_ref": None,
                "dK_dT_V": 0.001,
            },
            "fixed_parameters": ["Tr"],
        }
        record["validity"] = {
            "pressure_gpa": [0.0, 222.44],
            "temperature_k": [300.0, 3000.0],
            "volume_ratio": [0.66, 1.0],
            "notes": [
                "Table V calculation grid; a reference parameterization, not an independent experimental envelope."
            ],
        }
        record["parameter_provenance"] = {
            "reference_isotherm": {
                "V0": (
                    "page 1535 ambient density 19.30 g/cm^3 and atomic mass "
                    "196.967 g/mol; converted to the four-atom fcc cell"
                ),
                "K0": "Equation (29) and Table V; adopted KT0=166.65 GPa",
                "K0_prime": ("Equation (29) and Table V; adopted (dKT/dP)T=5.4823"),
            },
            "thermal_correction": {
                "Tr": "Equations (21), (27b), and (29); 300 K",
                "alpha_KT_ref": ("Equations (20), (21), (26), and (29); 7.14e-3 GPa/K"),
                "dK_dT_V": ("Section III and Equations (26)-(29); -11.5e-3 GPa/K"),
            },
        }
        record["notes"] = (
            "Anderson et al. Equation (29): the adopted 300 K BM3 curve "
            "KT0=166.65 GPa and K0'=5.4823 is combined with "
            "[0.00714-0.0115*ln(V0/V)]*(T-300) GPa. V0=67.79 A^3 is "
            "converted from the paper's ambient density 19.30 g/cm^3 and "
            "atomic mass 196.967 g/mol for the four-atom fcc cell. Table V "
            "covers V/V0=0.66-1 and 300-3000 K. The stored 0.001 GPa/K "
            "dK_dT_V error is the explicit propagated contribution reported "
            "for the near-identical K0'=5.5 calculation; the paper notes an "
            "additional unquantified contribution from K0' uncertainty. No "
            "complete covariance matrix or errors for the adopted static "
            "parameters are published."
        )
        for correction in (
            {
                "path": "eos.parameters.V0",
                "source_value": old_v0,
                "value": 67.79,
                "reason": (
                    "Convert the primary ambient density and atomic mass to the "
                    "four-atom conventional fcc cell; do not inherit another "
                    "catalog's reference volume."
                ),
            },
            {
                "path": "fixed_parameters",
                "source_value": [],
                "value": ["V0", "K0", "K0_prime"],
                "reason": (
                    "Equation (29) and Table V adopt these reference-isotherm "
                    "inputs rather than fitting them in this paper."
                ),
            },
            {
                "path": "thermal",
                "source_value": old_thermal,
                "value": dict(record["thermal"]),
                "reason": (
                    "Equation (29) is an additive linear thermal pressure whose "
                    "slope varies with ln(V0/V), not a shifted reference-state "
                    "EOS with constant expansivity and dK_dT=0."
                ),
            },
            {
                "path": "validity",
                "source_value": None,
                "value": dict(record["validity"]),
                "reason": "Use the complete calculation grid printed in Table V.",
            },
        ):
            correction["primary_reference"] = {
                "doi": "10.1063/1.342969",
                "location": "pages 1535-1541, Equations (20)-(29), and Table V",
            }
            append_correction(record, correction)

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
        record["fixed_parameters"] = ["V0"]
        record["experimental_pressure_range_gpa"] = [0.0, 46.0]
        record["temperature_ref"] = 300.0
        record["notes"] = (
            "Room-temperature BM3 fit to relative-volume X-ray diffraction data "
            "through 46 GPa: K0=307+/-5 GPa and K0'=6.2+/-0.3. Because the "
            "source fits and plots V/V0 and does not report V0 as a fitted "
            "coefficient, the refit fixes V0 at the value calculated from this "
            "sample's measured ambient hexagonal subcell a=3.0128(4) A and "
            "c=4.7357(9) A, converted to the equivalent four-formula-unit "
            "orthorhombic conventional-cell volume; its error is propagated from "
            "the reported lattice-parameter errors. The dedicated literature "
            "reproduction documents the remaining high-pressure sensitivity."
        )
        for correction in (
            {
                "path": "fixed_parameters",
                "source_value": [],
                "value": ["V0"],
                "reason": (
                    "Figure 2 and its caption express the fitted observations as "
                    "V/V0 and publish only K0 and K0'; the measured ambient "
                    "reference volume is therefore held fixed during the "
                    "reproduction."
                ),
            },
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
        record["fixed_parameters"] = ["V0"]
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
            "Walker et al.'s preferred bold Table 3 B2-KCl BE1 result is staged: "
            "the fictive V0=53.53 A^3 reference is held fixed while K0=23.7 GPa "
            "and K0'=4.4 are fitted to the room-temperature data, then Equation "
            "BE1 fits the directly reported alpha0*K0=0.0275+/-0.0009 kbar/K "
            "thermal-pressure coefficient to the P-V-T data. The separate italic "
            "Table 3 row is a simultaneous four-parameter solution (V0=55.25 A^3, "
            "K0=14.8 GPa, K0'=6.9, alpha0=0.00018 K^-1) and is not this record. "
            "The paper states that individual parameter errors are not meaningful "
            "because of correlation, so only the published error of the identifiable "
            "alpha0*K0 product is retained. The represented data span 3.18-8.14 GPa "
            "and 23-600 degC."
        )
        for correction in (
            {
                "path": "fixed_parameters",
                "source_value": [],
                "value": ["V0"],
                "reason": (
                    "The preferred bold Table 3 protocol fits only K0 and K0' to "
                    "the room-temperature data before fitting alpha0; V0 is the "
                    "held fictive reference. The italic row is the separate "
                    "simultaneous fit."
                ),
            },
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

    if identifier == "kcl_b2_tateno_2019_vinet_4":
        record["reference"]["authors"] = [
            "Tateno",
            "Komabayashi",
            "Hirose",
            "Hirao",
            "Ohishi",
        ]
        record["thermal"]["parameters"].update({"gamma0": 2.3, "q": 0.8})
        record["thermal"]["parameter_errors"].update({"gamma0": 0.2, "q": 0.2})
        record["thermal"]["debye_temperature_law"] = "integrated_gruneisen"
        record["notes"] = (
            "Preferred final-publication Tateno fit on the Sokolova Pt scale. The "
            "final article reports gamma0=2.3+/-0.2 and q=0.8+/-0.2 and defines "
            "theta(V)=theta0*exp[(gamma0-gamma(V))/q], the integrated-Gruneisen "
            "Debye law. The accepted manuscript instead printed gamma0=0.58+/-0.05 "
            "and q=0.9+/-0.2; those superseded values are not used. Table 1 also "
            "gives a Holmes-Pt alternative, which is not represented by this record. "
            "All 39 observations come from the official MSA Supplemental Table S1 "
            "workbook."
        )
        for correction in (
            {
                "path": "thermal.parameters",
                "source_value": {"gamma0": 0.58, "q": 0.9},
                "value": {"gamma0": 2.3, "q": 0.8},
                "reason": (
                    "Use the final published thermal coefficients rather than the "
                    "superseded accepted-manuscript values."
                ),
                "primary_reference": {
                    "doi": "10.2138/am-2019-6779",
                    "location": (
                        "Equation 6, Table 1, and thermal-EOS discussion on pages "
                        "721-722"
                    ),
                },
            },
            {
                "path": "thermal.debye_temperature_law",
                "source_value": "variable_exponent",
                "value": "integrated_gruneisen",
                "reason": (
                    "Final Equation 6 defines theta(V)=theta0*exp[(gamma0-gamma(V))/q]."
                ),
                "primary_reference": {
                    "doi": "10.2138/am-2019-6779",
                    "location": "Equation 6",
                },
            },
            {
                "path": "scientific_validation.primary_data_check",
                "source_value": "accepted-manuscript split table",
                "value": "official MSA Supplemental Table S1 workbook",
                "reason": (
                    "The split manuscript layout does not preserve Pt-KCl row "
                    "alignment for runs 3 and 4; the official XLSX deposit does."
                ),
                "primary_reference": {
                    "doi": "10.2138/am-2019-6779",
                    "location": ("MSA deposit AM-19-56779, 6779TableS1 revised.xlsx"),
                },
            },
        ):
            append_correction(record, correction)

    if identifier == "kcl_b2_chidester_2021_bm3_5":
        record["thermal"]["debye_temperature_law"] = "integrated_gruneisen"
        record["fit_datasets"] = [
            "kcl_dewaele_2012_table1_compression",
            "kcl_chidester_2021_supplemental_pvt",
        ]
        record["notes"] = (
            "The preferred Birch-Murnaghan reference EOS from Table I. Published "
            "molar V0=32.0(3) cm^3/mol is converted to one B2 formula unit per "
            "conventional cell with the exact Avogadro constant. Chidester et al. "
            "fitted the 155 new high-temperature observations simultaneously with "
            "all 123 Dewaele et al. (2012) room-temperature B2 observations. "
            "Pressures in the author-deposited high-temperature table were "
            "calculated from the Dorogokupets-Oganov (2007) Pt EOS. This is an "
            "experimental, effective-temperature-calibrated P-V-T pressure scale: "
            "KCl volumes and Pt-derived pressures are measured constraints, while "
            "the KCl temperature coordinate comes from the authors' gradient model "
            "for this laser-heated DAC geometry."
        )
        for correction in (
            {
                "path": "fit_datasets",
                "source_value": None,
                "value": list(record["fit_datasets"]),
                "reason": (
                    "Section III states that the new high-temperature data were "
                    "fitted together with the Dewaele et al. (2012) room-temperature "
                    "B2 data."
                ),
                "primary_reference": {
                    "doi": "10.1103/PhysRevB.104.094107",
                    "location": "Section III and Figure 3",
                },
            },
            {
                "path": "thermal.debye_temperature_law",
                "source_value": "variable_exponent",
                "value": "integrated_gruneisen",
                "reason": (
                    "Equations 3-4 leave theta(V) implicit; the thermodynamically "
                    "integrated constant-q relation reproduces all five Table I "
                    "coefficients with the complete 278-row fit scope."
                ),
                "primary_reference": {
                    "doi": "10.1103/PhysRevB.104.094107",
                    "location": "Equations 3-4, Table I, and Figure 3",
                },
            },
        ):
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
    primary_data_check = previous.get("primary_data_check")
    reproduction = previous.get("reproduction")
    audit_date = (
        REPORT_AUDIT_DATE
        if result["identifier"]
        in DERIVED_REFIT_RECORDS | CURRENT_SOURCE_AUDIT_RECORDS
        else CATALOG_AUDIT_DATE
        if result["identifier"] in DERIVED_REFERENCE_ISOTHERM_RECORDS
        else AUDIT_DATE
    )

    previous_evidence = (
        previous.get("primary_source_check")
        if previous.get("status") == "primary_source_validated"
        else None
    )
    record_evidence = VALIDATED_RECORD_SOURCES.get(result["identifier"])
    if record_evidence is None and previous_evidence is not None:
        # Native records added after the original migration audit carry their
        # completed, publication-specific evidence in the document itself.
        # Preserve that evidence when deterministically rebuilding the ledger.
        record_evidence = previous_evidence
    doi_evidence = VALIDATED_SOURCES.get(doi) if doi is not None else None
    if result["identifier"] not in FORCED_DEFERRED and (
        record_evidence is not None or doi_evidence is not None
    ):
        evidence = dict(record_evidence or doi_evidence or {})
        evidence["doi"] = doi
        validation: dict[str, Any] = {
            "status": "primary_source_validated",
            "note": (
                "Primary observations and source protocol were validated; the "
                "stored coefficients are an explicitly identified Peritheos refit."
                if result["identifier"] in DERIVED_REFIT_RECORDS
                else (
                    "Derived composition of separately primary-source-validated "
                    "reference-isotherm and simulated thermal components."
                    if result["identifier"] in DERIVED_REFERENCE_ISOTHERM_RECORDS
                    else (
                        "Independently checked against the cited primary publication "
                        "(published article or official/author copy); no external "
                        "software catalog was used as scientific authority."
                    )
                )
            ),
            "audit_date": audit_date,
            "verified_fields": list(VERIFIED_FIELDS),
            "primary_source_check": evidence,
        }
    else:
        reason = deferred_reason(result, doi)
        validation = {
            "status": "deferred",
            "note": reason,
            "audit_date": audit_date,
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
    if primary_data_check is not None:
        validation["primary_data_check"] = primary_data_check
    for extension in ("reported_parameterizations", "parameterization_resolution"):
        if extension in previous:
            validation[extension] = previous[extension]
    result["scientific_validation"] = validation

    if result["identifier"] == "ca_perovskite_caracas_2005_bm3_3":
        result["scientific_validation"]["note"] = (
            "The primary article and publisher HTML were audited directly. "
            "Exactly one cubic source parameterization is executable; all 18 "
            "Table 2 fits remain distinguished in the audit metadata."
        )
        result["scientific_validation"]["verified_fields"] = [
            "equation",
            "parameters",
            "units",
            "reference_state",
            "phase",
            "crystallography",
            "published_uncertainties",
            "validity",
            "source_parameterizations",
            "numerical_reproduction",
        ]

    if result["identifier"] == "kcl_b2_tateno_2019_vinet_4":
        result["scientific_validation"]["note"] = (
            "Equation, final published parameters, fixed quantities, uncertainties, "
            "P-T range, Pt pressure basis, and row alignment checked against the "
            "final article and official MSA supplement."
        )
        result["scientific_validation"]["audit_date"] = REPORT_AUDIT_DATE
        result["scientific_validation"]["verified_fields"].append(
            "pressure_calibration"
        )

    if result["identifier"] == "mgo_b1_luo_2023_vinet_thermal_5":
        result["scientific_validation"]["note"] = previous["note"]
        result["scientific_validation"]["verified_fields"] = previous[
            "verified_fields"
        ]
        if reproduction is not None:
            result["scientific_validation"]["reproduction"] = reproduction

    if result["identifier"] == "kcl_b2_chidester_2021_bm3_5":
        result["scientific_validation"]["note"] = (
            "Equation, parameters, units, molar-to-cell conversion, P-T range, "
            "temperature convention, Pt pressure basis, complete two-dataset fit "
            "scope, and Debye-temperature convention checked against the primary "
            "paper, Dewaele et al. (2012) Table I, and the author-deposited "
            "high-temperature table."
        )
        result["scientific_validation"]["audit_date"] = REPORT_AUDIT_DATE
        result["scientific_validation"]["verified_fields"].extend(
            ["pressure_calibration", "fit_scope", "thermal_model_convention"]
        )

    if result["identifier"] == "goethite_gleason_2008_bm3_1":
        result["scientific_validation"]["note"] = (
            "Independently checked against the cited primary publication and "
            "official deposit; bibliographic validation does not imply fit "
            "reproducibility. The EOS is retained for provenance but is not "
            "recommended for quantitative use because table and figure-marker "
            "refits fail to reproduce the published coefficients."
        )
        result["scientific_validation"]["usage_recommendation"] = (
            "not_recommended_for_quantitative_use"
        )
        result["scientific_validation"]["audit_date"] = REPORT_AUDIT_DATE
        result["scientific_validation"]["verified_fields"].append(
            "fit_reproducibility"
        )

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

    if result["identifier"] == "aragonite_martinez_1996_bm2_2":
        result["scientific_validation"]["reported_inconsistencies"] = [
            {
                "field": "parameter_errors.K0",
                "table_6": "3.48 GPa",
                "text_and_table_7": "approximately 4.3 GPa",
                "resolution": (
                    "Retain 3.48 GPa because it is the fit-specific sigma K0,T "
                    "printed in Table 6 alongside V0=227.5(8) A^3."
                ),
            }
        ]

    return result


def main() -> None:
    curate_migrated_catalog()
    entries: list[dict[str, Any]] = []
    for path in sorted(MATERIALS.glob("*.eosmat")):
        document = json.loads(path.read_text(encoding="utf-8"))
        primary_phase_labels = {
            "cscl": "B2 (CsCl-type), cubic Pm-3m",
            "fe3o4": "magnetite, cubic inverse-spinel Fd-3m",
            "nis": "metastable NiAs-type, hexagonal P63/mmc",
            "rbcl": "B2 (CsCl-type), cubic Pm-3m",
            "sno2_cubic_27gpa": "high-pressure cubic Pa-3 SnO2",
            "sno2_pa_3_at_48gpa": "high-pressure cubic Pa-3 SnO2",
            "sro": "B1 (NaCl-type), cubic Fm-3m",
            "sro_b2": "B2 (CsCl-type), cubic Pm-3m",
        }
        if document["identifier"] in primary_phase_labels:
            document["phase"] = primary_phase_labels[document["identifier"]]
        if document["identifier"] == "coo":
            document["lattice"].update({"a": 4.258, "b": 4.258, "c": 4.258})
            document["notes"] = (
                "Rocksalt CoO structure. The ambient cubic lattice parameter "
                "a0=4.258 A is from Clendenen and Drickamer (1966), Table II; "
                "the atom sites and space-group metadata are retained from the "
                "Dioptas/JCPDS interchange record."
            )
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
        if document["identifier"] == "feo":
            document["notes"] = (
                "Consolidated B1 FeO entry from equivalent cubic beamline JCPDS "
                "cards. The executable EOS is composition-specific to Fe0.94O. "
                "The legacy 'Fei 1995' static record was removed because its "
                "cited DOI is a thermal-expansion chapter and does not establish "
                "that EOS."
            )
        if document["identifier"] == "aragonite":
            document["notes"] = (
                "Aragonite structure with the reproducible staged BM2 P-V-T "
                "parameterization of Martinez et al. (1996). The paper's "
                "separate global HT-BM3 reduction is intentionally excluded: "
                "its fitted V0(298 K) is omitted, and its remaining coefficients "
                "cannot be reproduced from the printed 64-point table with the "
                "documented equations under pressure- or volume-residual least "
                "squares. Ambient coordinates are from Ye et al., American "
                "Mineralogist 97, 707-712 (2012)."
            )
        if document["identifier"] == "mgfe60o":
            document["name"] = "(Mg0.4Fe0.6)O (MW60)"
            document["formula"] = "Mg0.4Fe0.6O"
            document["phase"] = "rocksalt magnesiowustite, cubic Fm-3m"
            document["formula_units_per_cell"] = 4
            document["space_group"] = "Fm-3m"
            document["space_group_number"] = 225
            document["cell_contents"] = (
                "4 (Mg0.4Fe0.6)O formula units per conventional rocksalt cell"
            )
            document["notes"] = (
                "MW60 magnesiowustite containing 60 mol% FeO. The EOS volume "
                "uses the four-cation/four-anion conventional rocksalt cell. "
                "The legacy diffraction peaks are retained for interchange."
            )
        if document["identifier"] == "li_bcc":
            document["phase"] = "bcc structure; EOS fit spans bcc and fcc Li"
            document["cell_contents"] = (
                "2 Li atoms per bcc conventional cell; fcc observations in the "
                "combined EOS are expressed as equivalent two-atom volumes"
            )
            document["notes"] = (
                "Ambient bcc lithium crystallography. The bundled literature "
                "EOS is the primary paper's empirical bcc-fcc combined fit, not "
                "a single-phase bcc curve; its volume convention is documented "
                "at record level."
            )
        if document["identifier"] == "tungsten":
            document["notes"] = (
                "Tungsten pressure calibrant; the static-DAC Dewaele et al. "
                "(2004) Vinet fit is the default, with Cu-referenced Shen and "
                "Smith and Sokolova et al. fits as alternatives. The legacy "
                "Hixson-Fritz BM3 record was removed because the cited shock "
                "paper publishes Hugoniot-derived tabular isotherms, not the "
                "migrated standalone BM3 reduction."
            )
        records = []
        for record in document["eos_records"]:
            if record["identifier"] in REMOVED_UNSUPPORTED_RECORDS:
                continue
            audited = audit_record(record, path.name)
            records.append(audited)
            check = audited["scientific_validation"]
            entry = {
                "material": document["identifier"],
                "file": path.name,
                "record": audited["identifier"],
                "label": audited["label"],
                "doi": normalized_doi(audited.get("reference")),
                "status": check["status"],
                "note": check["note"],
                "primary_source_check": check["primary_source_check"],
            }
            if "usage_recommendation" in check:
                entry["usage_recommendation"] = check["usage_recommendation"]
            entries.append(entry)
        document["eos_records"] = records
        path.write_text(
            json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )

    counts = Counter(entry["status"] for entry in entries)
    if len(entries) != 163:
        raise ValueError(f"Expected 163 EOS records, found {len(entries)}")
    if "pending_primary_source_check" in counts:
        raise ValueError("Primary-source audit left pending records")

    report = {
        "format": "peritheos.primary-source-audit",
        "format_version": 1,
        "audit_date": REPORT_AUDIT_DATE,
        "policy": {
            "scientific_authority": "primary publications and official supplements",
            "external_catalogs": (
                "External software catalogs are not used as scientific authority. "
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
    manifest["materials"] = len(list(MATERIALS.glob("*.eosmat")))
    manifest["eos_records"] = len(entries)
    manifest["scientific_validation"] = {
        "audit_date": REPORT_AUDIT_DATE,
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
