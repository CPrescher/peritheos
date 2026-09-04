#!/usr/bin/env python3
"""Refit every bundled EOS record for which primary observations permit it.

The generated ledger is intentionally conservative.  A record without row-level
observations, or whose pressures require a calibration that Peritheos cannot yet
execute, is reported as ``not_refittable`` instead of being fitted to invented
or circularly reconstructed data.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from peritheos import get_material_document, list_material_documents
from peritheos.eos import ThermalEOS
from peritheos.eos.rt import BM2, BM3, BM4, Murnaghan, Vinet
from peritheos.eos.thermal import ThermalReferenceStateEOS
from peritheos.fitting import fit_joint_eos, fit_rt_eos
from peritheos.materials import Material

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "peritheos" / "data"
DEFAULT_JSON = ROOT / "docs" / "data" / "primary-eos-refits.json"
DEFAULT_MARKDOWN = ROOT / "docs" / "primary-eos-refits.md"

MODEL_CLASSES = {
    "BM2": BM2,
    "BM3": BM3,
    "BM4": BM4,
    "Murnaghan": Murnaghan,
    "Vinet": Vinet,
}

# These observations do not define the pressure-volume fit stored by the record.
INDIRECT_DATA = {
    "platinum_holmes_1989_vinet_1": (
        "The bundled rows are shock-Hugoniot qualification experiments; the stored "
        "equilibrium Vinet curve is a theoretical 300 K isotherm and cannot be "
        "refitted directly to those rows."
    ),
    "mgsio3_post_perovskite_mosenfelder_2009_bm3_1": (
        "The bundled rows are shock states and the source's thermal reduction cannot "
        "be reconstructed as a direct P-V-T least-squares fit because most rows do "
        "not report temperature."
    ),
    "nickel_oxide_noguchi_1999_bm3_1": (
        "The bundled rows are Hugoniot states; the stored 300 K isotherm is the "
        "source's Mie-Gruneisen reduction, not a direct fit to Hugoniot P-V pairs."
    ),
    "diamond_correa_2008_dewaele_anchored": (
        "The linked diffraction rows constrain only the Dewaele reference isotherm; "
        "the Correa thermal term is a separately published theoretical model."
    ),
    "diamond_benedict_2014_dewaele_anchored": (
        "The linked diffraction rows constrain only the Dewaele reference isotherm; "
        "the Benedict thermal term is a separately published theoretical model."
    ),
    "iridium_anzellini_2025_bm3_1": (
        "The bundled rows are all heated states. The stored coefficients are the "
        "300 K reference part of a combined thermal fit, but the record does not "
        "represent the source's thermal correction needed to refit those rows."
    ),
}

SHEN_PREFIX = "shen_smith_2026_table_s1_simultaneous_volumes"

# The Fei et al. neon points are digitized from a plot with no experimental
# error bars. The stored uncertainties describe the digitization, not the
# weights used in the paper's least-squares fit.
UNWEIGHTED_DATASETS = {"neon_fei_2007_figure5_digitized"}

COMBINED_FIT_DATASET_RECORDS = {
    "cscl_campbell_1994_bm3_1",
    "neon_fcc_fei_2007_bm3_1",
    "neon_fcc_fei_2007_vinet_2",
}

# These tables report one-sigma uncertainties on a cubic lattice parameter.
# Convert sigma_a to sigma_V with the first-order derivative of V=a^3.
CUBIC_LATTICE_SIGMA_DATASETS = {
    "cscl_campbell_1994_table1_compression",
    "rbcl_campbell_1994_table1_compression",
}

FIT_QUALIFICATIONS = {
    "cscl_campbell_1994_bm3_1": (
        "Complete source-data reproduction: Table 1 contains a distinct 13-row "
        "CsCl block, not the RbCl rows previously attached to this material. The "
        "fit now also includes all nine 25 degC V/V0 values from Yagi (1978) after "
        "the source-stated 4.118-to-4.123 A reference correction. The numerical "
        "weights used for Campbell and Heinz's normalized-stress regression are not "
        "published, but the complete unweighted 22-point refit recovers both "
        "coefficients within combined two-sigma uncertainty. See the "
        "[dedicated Campbell-Heinz reproduction]"
        "(literature-reproductions.md#campbell-heinz-1994-cscl-and-rbcl)."
    ),
    "kcl_b2_chidester_2021_bm3_5": (
        "Corrected source-scope reproduction: Chidester et al. fitted the 155 "
        "new high-temperature rows simultaneously with all 123 Dewaele et al. "
        "(2012) room-temperature B2 rows. The unweighted 278-row fit uses the "
        "thermodynamically integrated Gruneisen Debye-temperature law and "
        "recovers every published coefficient. See the [dedicated Chidester "
        "reproduction](literature-reproductions.md#kcl-chidester-2021)."
    ),
    "kcl_b2_tateno_2019_vinet_4": (
        "Corrected final-publication reproduction: the model uses the final "
        "gamma0=2.3 and q=0.8 coefficients, Equation 6's integrated-Gruneisen "
        "Debye law, and the correctly aligned official MSA Supplemental Table S1 "
        "workbook. See the [dedicated Tateno reproduction]"
        "(literature-reproductions.md#kcl-tateno-2019)."
    ),
    "kcl_walker_2002_bm3_2": (
        "Source-protocol reproduction: the preferred bold Table 3 result is "
        "staged, not a simultaneous four-parameter fit. V0 is held at the "
        "published fictive reference, K0 and K0' are fitted to the eight "
        "23-24 degC rows, and alpha_KT is then fitted to all 39 P-V-T rows."
    ),
    "neon_fcc_fei_2007_bm3_1": (
        "Conditional partial reproduction: Fei et al. fitted Hemley et al. (1989, "
        "ref. 45), Finger et al. (1981, ref. 47), and their new observations. The "
        "34 bundled points combine the 13 new Figure 5 observations with all 21 "
        "Hemley Table I observations after the source-stated Mao-to-Dewaele ruby "
        "recalculation. Finger's low-pressure rows remain unavailable. V0 is held "
        "fixed and the available subset is fitted without weights."
    ),
    "neon_fcc_fei_2007_vinet_2": (
        "Conditional partial reference-isotherm reproduction: Fei et al. fitted "
        "Hemley et al. (1989, ref. 45), Finger et al. (1981, ref. 47), and their new "
        "observations. The 34 bundled points combine the 13 new Figure 5 observations "
        "with all 21 Hemley Table I observations after the source-stated "
        "Mao-to-Dewaele ruby recalculation. Finger's low-pressure rows remain "
        "unavailable. V0 and thermal coefficients are held fixed, and the available "
        "subset is fitted without weights."
    ),
    "rbcl_b2_campbell_1994_bm3_1": (
        "Complete source-data reproduction: all 24 RbCl-B2 Table 1 rows are fitted "
        "with the paper's hypothetical zero-pressure density held fixed. The refit "
        "recovers K0 and K0' within combined two-sigma uncertainty. See the "
        "[dedicated Campbell-Heinz reproduction]"
        "(literature-reproductions.md#campbell-heinz-1994-cscl-and-rbcl)."
    ),
}

# Record-specific conclusions for cases that remain outside both the numerical
# and uncertainty criteria. These are deliberately phrased as diagnoses or
# bounded hypotheses; unresolved source-method details are not presented as fact.
INVESTIGATION_NOTES = {
    "b4c_somayazulu_2023_bm3_1": (
        "The reference isotherm and gamma0 are fixed, leaving q as the only free "
        "coefficient. The published curve and the refit have similar pressure "
        "residuals, but q moves from 2.1 to about 1.05. This points to weak q "
        "identifiability and sensitivity to the source's exact effective-variance "
        "objective, rather than an equation-evaluation error. The existing B4C "
        "literature reproduction documents that distinction in detail."
    ),
    "coo_clendenen_1966_murnaghan_1": (
        "The input is a nine-point smoothed and rounded table, not the original "
        "measurements, and the paper gives no coefficient covariance or exact "
        "regression objective. K0 remains close, while K0' is the sensitive "
        "coefficient. Reproducing the authors' pre-smoothing data and regression "
        "variable is the next evidence needed."
    ),
    "cscl_campbell_1994_bm3_1": (
        "The original extreme discrepancy was a data-assignment error: the resource "
        "contained 21 rows from Table 1's RbCl block and evaluated them with CsCl's "
        "larger reference volume. After restoring the 13-row CsCl block and adding "
        "Yagi's nine Table 1 ratios with the Campbell-Heinz reference correction, "
        "the complete 22-point refit is statistically compatible with the published "
        "fit. Only the exact normalized-stress regression weights remain unpublished."
    ),
    "goethite_gleason_2008_bm3_1": (
        "The audit uses the 27 rows nearest the reference temperature from a "
        "65-row P-V-T deposit. The source's exact temperature selection, thermal "
        "correction, weighting, and covariance are not encoded. The much softer "
        "K0/K0' trade-off in the refit is therefore most plausibly a fit-scope "
        "difference; reproducing the paper's row mask and objective is required."
    ),
    "ice_vi_bezacier_2014_bm2_1": (
        "The 42 K temperature span and roughly 1.3 GPa pressure span provide weak "
        "leverage on alpha0 once V0 and K0 are also free. The refitted curve only "
        "modestly improves the RMSE, but alpha0 shifts strongly. This is consistent "
        "with coefficient correlation; the source covariance or exact fitting "
        "weights are needed to distinguish it from a convention difference."
    ),
    "kcl_walker_2002_bm3_2": (
        "The original large discrepancy was a validation-protocol error: it compared "
        "the paper's preferred staged Table 3 row with an unconstrained simultaneous "
        "four-parameter refit. Reproducing the Table 3 footnote instead gives "
        "K0=23.775 GPa and K0'=4.416 from the eight room-temperature rows, followed "
        "by alpha_KT=0.002766 GPa/K from all 39 rows. These differ from the printed "
        "23.7, 4.4, and 0.00275 by 0.32%, 0.36%, and 0.59%, respectively. Strict "
        "uncertainty parity remains unavailable only because Walker et al. decline "
        "individual elastic-parameter errors due to covariance. The paper separately "
        "prints an italic simultaneous solution (V0=55.25 A^3, K0=14.8 GPa, "
        "K0'=6.9, alpha0=0.00018 K^-1), confirming that the two protocols must not "
        "be conflated. The corresponding unweighted Peritheos simultaneous fit also "
        "recovers that alternate row closely (V0=55.392 A^3, K0=14.189 GPa, "
        "K0'=7.157). See the [dedicated Walker reproduction]"
        "(literature-reproductions.md#kcl-walker-2002)."
    ),
    "kcl_b2_tateno_2019_vinet_4": (
        "The original failure combined two source-control errors. The record used "
        "gamma0=0.58, q=0.9, and a variable-exponent Debye law from the accepted "
        "manuscript, while the final article reports gamma0=2.3, q=0.8, and the "
        "integrated-Gruneisen law in Equation 6. In addition, the split manuscript "
        "table had mismatched the Pt and KCl halves of runs 3 and 4. The corrected "
        "official MSA Table S1 data recover all four fitted coefficients within "
        "combined two-sigma uncertainty. See the [dedicated Tateno reproduction]"
        "(literature-reproductions.md#kcl-tateno-2019)."
    ),
    "molybenum_carbide_mo2c_haines_2001_bm3_1": (
        "Only Figure 2 markers are available, split across pressure-medium regimes. "
        "K0 and V0 remain broadly consistent, but K0' is sensitive to marker "
        "placement, weighting, and the highest-pressure points. The original numeric "
        "P-V array or the authors' finite-strain regression is needed for a decisive "
        "coefficient comparison."
    ),
    "palladium_baty_2024_bm3_1": (
        "The complete supplementary table is present, yet the published curve fits "
        "the checked rows substantially less well than the refit and K0 changes by "
        "about 20%. This makes source row selection, run grouping, pressure-scale "
        "treatment, or the fitted residual variable more likely than transcription "
        "rounding. The publication's exact fitting script or row mask is needed."
    ),
    "silicon_vii_anzellini_2019_vinet_1": (
        "Si-VII is observed only at 46-94 GPa, so all three zero-pressure Vinet "
        "coefficients are extrapolated and strongly covariant. The refit drives K0 "
        "to its lower bound, a direct sign that this isolated phase subset does not "
        "identify an unconstrained zero-pressure EOS. The source's constraints or "
        "joint phase-fit protocol must be recovered."
    ),
    "sno2_cubic_27gpa_ono_2000_bm3_1": (
        "Only eleven 300 K rows over 16-29 GPa constrain an ambient reference state "
        "that the phase cannot retain on decompression. K0' reaches its lower bound "
        "and K0 compensates, demonstrating extrapolation-driven non-identifiability. "
        "The source's constraints, covariance, and any use of heated rows must be "
        "reproduced before interpreting the coefficient shift."
    ),
    "sno2_pa_3_at_48gpa_ono_2000_bm3_1": (
        "This alias uses the same eleven 300 K rows and the same parameterization as "
        "the cubic SnO2 record, so it inherits the identical extrapolation problem: "
        "K0' reaches zero and K0 compensates. It is not an independent failed "
        "experiment; resolution requires the same source constraints and fit scope."
    ),
    "wadsleyite_katsura_2009_bm3_1": (
        "The full 85-row table is present, but the published parameterization misses "
        "those rows by much more than the refit and gamma0 shifts far outside its "
        "reported error. The pressure calibration remains unresolved, and an MGD "
        "reference-state or Debye-energy convention mismatch is also possible. A "
        "row-by-row source-equation reproduction should precede any coefficient "
        "revision."
    ),
}

# Dataset choices that cannot be inferred uniquely from generic quantity metadata.
PRESSURE_COLUMNS = {
    "neon_hemley_1989_table1_compression": "pressure_gpa",
    "neon_hemley_1989_table1_fei_recalculated": "pressure_gpa_dewaele_2004",
    "aluminum_dewaele_2004_table1_compression": "ruby_pressure_revised_gpa",
    "copper_dewaele_2004_table1_compression": "ruby_pressure_revised_gpa",
    "gold_dewaele_2004_table1_compression": "ruby_pressure_revised_gpa",
    "tungsten_dewaele_2004_table1_compression": "ruby_pressure_revised_gpa",
    "silver_dewaele_2008_table2_compression": "ruby_pressure_dewaele_gpa",
    "diamond_dewaele_2008_table1_pvt": "pressure_h05_gpa",
    "gold_takemura_2008_table3_compression": "ruby_pressure_dorogokupets_gpa",
    "gold_fratanduono_2021_table_s3_compression": "isotherm_298k_pressure_gpa",
    "nickel_dewaele_2008_table2_compression": "ruby_pressure_dewaele_gpa",
    "mgsio3_post_perovskite_sakai_2016_table_s1_pvt": "pressure_bm3_gpa",
    "rhodium_rodrigo_ramon_2024_table2_compression": "pressure_gpa",
    "chromium_anzellini_2022_table2_compression": "pressure_w_gpa",
    "ringwoodite_meng_1994_table1_pvt": "observed_pressure_gpa",
    "bridgmanite_tange_2012_table1_pvt#BM3": "pressure_bm3_gpa",
    "bridgmanite_tange_2012_table1_pvt#Vinet": "pressure_vinet_gpa",
    "mgfe_perovskite_knittle_1987_table2_compression": ("pressure_after_heating_gpa"),
}

VOLUME_COLUMNS = {
    "alumina_dewaele_2013_table1_compression": "a_a",
    "alpha_quartz_angel_1997_table1_compression": "unit_cell_volume_a3",
    # Table 1 normalizes both phases to the ambient B1 volume. Using the
    # measured cubic lattice parameter avoids applying the B2 record's V0 to
    # that B1-normalized column.
    "cao_richet_1988_table1_compression": "lattice_a_angstrom",
    "bridgmanite_tange_2012_table1_pvt": "unit_cell_volume_a3",
    "chromium_anzellini_2022_table2_compression": "chromium_unit_cell_volume_a3",
    "iridium_anzellini_2025_tables_s1_s3_pvt": "iridium_lattice_a_angstrom",
    "kcl_tateno_2019_table_s1_pvt": "kcl_unit_cell_volume_a3",
    "neon_hemley_1989_table1_compression": "volume_a3_conventional_cell",
    "mgsio3_post_perovskite_sakai_2016_table_s1_pvt": ("mgsio3_unit_cell_volume_a3"),
    "nis_campbell_1993_table1_compression": "a_a",
    "rhodium_rodrigo_ramon_2024_table2_compression": "rhodium_unit_cell_volume_a3",
    "ringwoodite_meng_1994_table1_pvt": "observed_volume_a3",
    "diamond_dewaele_2008_table1_pvt": "lattice_a_angstrom",
    "feo_fischer_2011_table_s1_pvt#feo_fischer_2011_bm3_2": (
        "b1_feo_molar_volume_cm3_mol"
    ),
    "feo_fischer_2011_table_s1_pvt#feo_b8_2_fischer_2011_bm3_1": (
        "b8_feo_molar_volume_cm3_mol"
    ),
    "silicon_anzellini_2019_tables1_4_6_7_compression": ("silicon_lattice_a_angstrom"),
    "silicon_carbide_b3_miozzi_2018_data_set_s1_eos": ("sic_unit_cell_volume_a3"),
    "titanium_alpha_dewaele_2015_table4_compression": "lattice_a_angstrom",
    "titanium_omega_dewaele_2015_table4_compression": "lattice_a_angstrom",
}

PHASE_FILTERS = {
    "cao_richet_1988_bm3_1": {"phase": "B1", "used_in_eos_fit": "yes"},
    "cao_b2_richet_1988_bm3_1": {"phase": "B2", "used_in_eos_fit": "yes"},
    "phase_d_ant_a_shieh_2000_bm2_1": {"sample": "1"},
    "phase_d_ant_b_shieh_2000_bm2_1": {"sample": "2"},
    "forsterite_finkelstein_2014_bm3_1": {
        "phase": "forsterite_I",
        "used_in_forsterite_i_eos_fit": "yes",
    },
    "silicon_anzellini_2019_vinet_1": {"phase": "I"},
    "silicon_v_anzellini_2019_vinet_1": {"phase": "V"},
    "silicon_vii_anzellini_2019_vinet_1": {"phase": "VII"},
    "silicon_x_anzellini_2019_vinet_1": {"phase": "X"},
    "sno2_hazen_1981_bm3_1": {"compound": "SnO2"},
    "geo2_rutile_hazen_1981_bm3_1": {"compound": "GeO2"},
}


@dataclass
class Series:
    dataset_id: str
    pressure: np.ndarray
    volume: np.ndarray
    temperature: np.ndarray | None
    pressure_sigma: np.ndarray | None
    volume_sigma: np.ndarray | None
    temperature_sigma: np.ndarray | None
    pressure_column: str
    volume_column: str
    temperature_column: str | None
    selection: str


def _number(value: Any) -> float:
    if value is None or value == "":
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _load_rows(dataset: dict[str, Any]) -> list[dict[str, Any]]:
    names = [column["name"] for column in dataset["columns"]]
    if "rows" in dataset:
        return [dict(zip(names, row)) for row in dataset["rows"]]
    path = DATA_ROOT / dataset["resource"]["path"]
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        next(reader)
        return [dict(zip(names, row)) for row in reader]


def _column_map(dataset: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {column["name"]: column for column in dataset["columns"]}


def _value_columns(dataset: dict[str, Any]) -> list[dict[str, Any]]:
    return [column for column in dataset["columns"] if column["role"] == "value"]


def _pressure_column(dataset: dict[str, Any], model_type: str) -> str:
    key = f"{dataset['identifier']}#{model_type}"
    if key in PRESSURE_COLUMNS:
        return PRESSURE_COLUMNS[key]
    if dataset["identifier"] in PRESSURE_COLUMNS:
        return PRESSURE_COLUMNS[dataset["identifier"]]
    candidates = [
        column
        for column in _value_columns(dataset)
        if "pressure" in column["quantity"]
        and "residual" not in column["quantity"]
        and "calculated" not in column["quantity"]
        and "before" not in column["quantity"]
        and "standard" not in column["quantity"]
        and "spatial" not in column["quantity"]
    ]
    exact = [column for column in candidates if column["quantity"] == "pressure"]
    if exact:
        candidates = exact
    if len(candidates) != 1:
        raise ValueError(
            "ambiguous pressure columns: "
            + ", ".join(column["name"] for column in candidates)
        )
    return candidates[0]["name"]


def _temperature_column(dataset: dict[str, Any]) -> str | None:
    candidates = [
        column
        for column in _value_columns(dataset)
        if column["quantity"] in {"temperature", "measured_temperature"}
    ]
    if not candidates:
        return None
    exact = [column for column in candidates if column["quantity"] == "temperature"]
    return (exact or candidates)[0]["name"]


def _volume_column(dataset: dict[str, Any], record_id: str) -> str:
    key = f"{dataset['identifier']}#{record_id}"
    if key in VOLUME_COLUMNS:
        return VOLUME_COLUMNS[key]
    if dataset["identifier"] in VOLUME_COLUMNS:
        return VOLUME_COLUMNS[dataset["identifier"]]
    candidates = []
    for column in _value_columns(dataset):
        quantity = column["quantity"]
        if any(
            token in quantity
            for token in (
                "volume",
                "density",
                "lattice_parameter",
            )
        ) and not any(
            token in quantity
            for token in (
                "pressure_standard",
                "calibrant",
                "residual",
                "calculated",
                "axial",
            )
        ):
            candidates.append(column)
    priorities = (
        "conventional_unit_cell_volume",
        "unit_cell_volume",
        "cell_volume",
        "volume",
        "atomic_volume",
        "molar_volume",
        "specific_volume",
        "volume_ratio",
        "relative_volume",
        "density",
        "lattice_parameter_ratio",
        "lattice_parameter",
    )
    for priority in priorities:
        selected = [column for column in candidates if column["quantity"] == priority]
        if len(selected) == 1:
            return selected[0]["name"]
        selected = [
            column
            for column in candidates
            if column["quantity"].endswith("_" + priority)
        ]
        if len(selected) == 1:
            return selected[0]["name"]
    if len(candidates) == 1:
        return candidates[0]["name"]
    raise ValueError(
        "ambiguous volume columns: "
        + ", ".join(column["name"] for column in candidates)
    )


def _sigma_column(dataset: dict[str, Any], value_name: str) -> str | None:
    columns = dataset["columns"]
    explicit = [column for column in columns if column.get("of") == value_name]
    if explicit:
        return explicit[0]["name"]
    value = _column_map(dataset)[value_name]
    quantity = value["quantity"]
    candidates = [
        column
        for column in columns
        if column["quantity"] == quantity
        and column["role"] in {"standard_deviation", "standard_error", "uncertainty"}
    ]
    return candidates[0]["name"] if len(candidates) == 1 else None


def _pressure_factor(unit: str) -> float:
    return {"GPa": 1.0, "kbar": 0.1, "Mbar": 100.0}[unit]


def _temperature_values(values: np.ndarray, unit: str) -> np.ndarray:
    if unit == "K":
        return values
    if unit in {"degC", "C", "celsius"}:
        return values + 273.15
    raise ValueError(f"unsupported temperature unit {unit!r}")


def _formula_atom_count(formula: str) -> float:
    tokens = re.findall(r"([A-Z][a-z]?)([0-9]*\.?[0-9]*)", formula)
    if not tokens:
        return 1.0
    return sum(float(count) if count else 1.0 for _, count in tokens)


def _volume_values(
    values: np.ndarray,
    column: dict[str, Any],
    document: dict[str, Any],
    record: dict[str, Any],
    rows: list[dict[str, Any]],
) -> tuple[np.ndarray, float]:
    quantity = column["quantity"]
    name = column["name"]
    unit = column["unit"]
    v0 = float(record["eos"]["parameters"]["V0"])
    formula_units = float(document.get("formula_units_per_cell") or 1.0)

    if "/atom" in unit:
        factor = formula_units * _formula_atom_count(str(document["formula"]))
        return values * factor, factor
    if "/formula_unit" in unit or quantity == "formula_unit_volume":
        return values * formula_units, formula_units
    if quantity in {"volume_ratio", "relative_volume", "normalized_volume"} or (
        "ratio" in quantity and "volume" in quantity
    ):
        if record["identifier"] == "kcl_campbell_1991_bm2_1":
            # The source normalizes B2 measurements to the ambient B1 volume,
            # not to the extrapolated B2 V0 stored by this record.
            return values * 62.36, 62.36
        if record["identifier"] == "nacl_b2_heinz_1984_bm2_1":
            # Heinz and Jeanloz tabulate V/V01 with V01 the B1 reference volume;
            # the record stores the B2 extrapolated V0 = 0.929 * V01.
            b1_v0 = v0 / 0.929
            return values * b1_v0, b1_v0
        return values * v0, v0
    if quantity == "lattice_parameter_ratio":
        return values**3 * v0, 3.0 * v0
    if "density" in quantity:
        if "sro_b2" in record["identifier"]:
            return v0 * 6.14 / values, math.nan
        raise ValueError("density-to-cell-volume conversion is not specified")
    if "molar_volume" in quantity:
        factor = formula_units / 0.602214076
        return values * factor, factor
    if "atomic_volume" in quantity:
        factor = formula_units * _formula_atom_count(str(document["formula"]))
        return values * factor, factor
    if "specific_volume" in quantity:
        factor = formula_units
        return values * factor, factor
    if "volume" in quantity or "volume" in name:
        return values, 1.0

    # Lattice-only source tables.  The material symmetry and available a/b/c
    # columns determine the conventional-cell volume.
    by_name = {
        key: np.array([_number(row.get(key)) for row in rows]) for key in rows[0]
    }
    lower = str(document.get("symmetry", "")).lower()
    a = values
    b_name = next((key for key in by_name if "lattice_b" in key or key == "b_a"), None)
    c_name = next((key for key in by_name if "lattice_c" in key or key == "c_a"), None)
    if "cubic" in lower or (b_name is None and c_name is None):
        return a**3, math.nan
    if "hexagonal" in lower or "trigonal" in lower:
        if c_name is not None and not np.any(np.isfinite(by_name[c_name])):
            c_name = b_name
        if c_name is None:
            raise ValueError("hexagonal lattice requires c")
        return (math.sqrt(3.0) / 2.0) * a**2 * by_name[c_name], math.nan
    if b_name is not None and c_name is not None:
        return a * by_name[b_name] * by_name[c_name], math.nan
    raise ValueError("cannot derive volume from lattice columns")


def _select_rows(
    rows: list[dict[str, Any]], record_id: str
) -> tuple[list[dict[str, Any]], str]:
    filters = PHASE_FILTERS.get(record_id, {})
    selected = [
        row
        for row in rows
        if all(str(row.get(name, "")) == expected for name, expected in filters.items())
    ]
    detail = ", ".join(f"{key}={value}" for key, value in filters.items()) or "all rows"
    return selected, detail


def _series(
    document: dict[str, Any], record: dict[str, Any], dataset: dict[str, Any]
) -> Series:
    rows, selection = _select_rows(_load_rows(dataset), record["identifier"])
    if not rows:
        raise ValueError("row selection is empty")
    model_type = record["eos"]["type"]
    if dataset["identifier"] == "silicon_anzellini_2019_tables1_4_6_7_compression":
        pressure_name = "pressure_gold_gpa|pressure_tungsten_gpa|pressure_ruby_gpa"
    else:
        pressure_name = _pressure_column(dataset, model_type)
    volume_name = _volume_column(dataset, record["identifier"])
    temperature_name = _temperature_column(dataset)
    columns = _column_map(dataset)

    if dataset["identifier"] == "silicon_anzellini_2019_tables1_4_6_7_compression":
        pressure_name = "pressure_gold_gpa|pressure_tungsten_gpa|pressure_ruby_gpa"
        pressure = np.array(
            [
                next(
                    (
                        value
                        for name in (
                            "pressure_gold_gpa",
                            "pressure_tungsten_gpa",
                            "pressure_ruby_gpa",
                        )
                        if np.isfinite(value := _number(row.get(name)))
                    ),
                    math.nan,
                )
                for row in rows
            ]
        )
        pressure_unit = "GPa"
    else:
        pressure = np.array([_number(row.get(pressure_name)) for row in rows])
        pressure_unit = columns[pressure_name]["unit"]
    pressure *= _pressure_factor(pressure_unit)
    raw_volume = np.array([_number(row.get(volume_name)) for row in rows])
    volume, volume_factor = _volume_values(
        raw_volume, columns[volume_name], document, record, rows
    )
    temperature = None
    if temperature_name is not None:
        temperature = np.array([_number(row.get(temperature_name)) for row in rows])
        temperature = _temperature_values(
            temperature, columns[temperature_name]["unit"]
        )

    def uncertainties(name: str, factor: float = 1.0) -> np.ndarray | None:
        sigma_name = _sigma_column(dataset, name)
        if sigma_name is None:
            return None
        result = np.array([_number(row.get(sigma_name)) for row in rows]) * factor
        return result

    pressure_sigma = (
        None
        if "|" in pressure_name
        else uncertainties(pressure_name, _pressure_factor(pressure_unit))
    )
    volume_sigma = uncertainties(volume_name, volume_factor)
    if dataset["identifier"] in CUBIC_LATTICE_SIGMA_DATASETS:
        lattice_sigma = uncertainties(volume_name)
        if lattice_sigma is not None:
            volume_sigma = 3.0 * raw_volume**2 * lattice_sigma
    temperature_sigma = (
        uncertainties(temperature_name) if temperature_name is not None else None
    )

    finite = np.isfinite(pressure) & np.isfinite(volume) & (volume > 0.0)
    if temperature is not None:
        finite &= np.isfinite(temperature) & (temperature > 0.0)
    pressure = pressure[finite]
    volume = volume[finite]
    temperature = None if temperature is None else temperature[finite]
    pressure_sigma = _usable_sigma(pressure_sigma, finite)
    volume_sigma = _usable_sigma(volume_sigma, finite)
    temperature_sigma = _usable_sigma(temperature_sigma, finite)
    if pressure.size < 3:
        raise ValueError(f"only {pressure.size} usable observations")
    return Series(
        dataset_id=dataset["identifier"],
        pressure=pressure,
        volume=volume,
        temperature=temperature,
        pressure_sigma=pressure_sigma,
        volume_sigma=volume_sigma,
        temperature_sigma=temperature_sigma,
        pressure_column=pressure_name,
        volume_column=volume_name,
        temperature_column=temperature_name,
        selection=selection,
    )


def _usable_sigma(values: np.ndarray | None, mask: np.ndarray) -> np.ndarray | None:
    if values is None:
        return None
    values = values[mask]
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        return None
    return values


def _bounds(name: str, value: float) -> tuple[float, float]:
    if name.endswith("V0"):
        return (max(value * 0.25, 1.0e-9), value * 2.0)
    if name.endswith("K0"):
        return (max(value * 0.05, 1.0e-6), value * 5.0)
    if name.endswith("K0_prime"):
        return (0.0, 20.0)
    if name.endswith("K0_double_prime"):
        return (-2.0, 2.0)
    if name in {"gamma0", "q", "theta0"}:
        return (1.0e-9, max(value * 10.0, 10.0))
    if name.startswith("alpha"):
        scale = max(abs(value) * 20.0, 1.0e-8)
        return (-scale, scale)
    if name == "dK_dT":
        return (-1.0, 1.0)
    scale = max(abs(value) * 20.0, 1.0)
    return (-scale, scale)


def _static_parameters(
    record: dict[str, Any],
) -> tuple[dict[str, float], dict[str, float]]:
    values = {
        name: float(value)
        for name, value in record["eos"]["parameters"].items()
        if name in {"V0", "K0", "K0_prime", "K0_double_prime"}
    }
    fixed_names = set(record.get("fixed_parameters", ()))
    fixed_names.update(record["eos"].get("fixed_parameters", ()))
    fixed = {name: value for name, value in values.items() if name in fixed_names}
    initial = {name: value for name, value in values.items() if name not in fixed_names}
    return initial, fixed


def _thermal_parameters(
    record: dict[str, Any],
) -> tuple[dict[str, float], dict[str, float]]:
    thermal = record["thermal"]
    values = {
        name: float(value)
        for name, value in thermal["parameters"].items()
        if value is not None and name not in {"phi0"}
    }
    fixed_names = set(thermal.get("fixed_parameters", ()))
    # Reference temperature and atom count are definitions even in older migrated
    # records whose original interchange metadata omitted an explicit fixed list.
    fixed_names.update({"Tr", "n"} & values.keys())
    fixed = {name: value for name, value in values.items() if name in fixed_names}
    initial = {name: value for name, value in values.items() if name not in fixed_names}
    return initial, fixed


def _configuration(record: dict[str, Any]) -> dict[str, Any]:
    thermal = record.get("thermal", {})
    result = dict(thermal.get("configuration", {}))
    for name in (
        "debye_temperature_law",
        "thermal_expansion_law",
        "reference_volume_law",
    ):
        if name in thermal:
            result[name] = thermal[name]
    return result


def _published_parameters(record: dict[str, Any], thermal: bool) -> dict[str, float]:
    values = {
        name: float(value)
        for name, value in record["eos"]["parameters"].items()
        if name in {"V0", "K0", "K0_prime", "K0_double_prime"}
    }
    if thermal:
        values = {f"rt_eos.{name}": value for name, value in values.items()}
        values.update(
            {
                name: float(value)
                for name, value in record["thermal"]["parameters"].items()
                if value is not None and name != "phi0"
            }
        )
    return values


def _published_errors(record: dict[str, Any], thermal: bool) -> dict[str, float]:
    errors = {
        name: float(value)
        for name, value in record.get("parameter_errors", {}).items()
        if value is not None
    }
    if thermal:
        errors = {f"rt_eos.{name}": value for name, value in errors.items()}
        errors.update(
            {
                name: float(value)
                for name, value in record["thermal"].get("parameter_errors", {}).items()
                if value is not None
            }
        )
    return errors


def _similar(name: str, published: float, fitted: float) -> bool:
    difference = abs(fitted - published)
    relative = difference / max(abs(published), 1.0e-30)
    if name.endswith("V0"):
        return relative <= 0.05
    if name.endswith("K0"):
        return relative <= 0.15
    if name.endswith("K0_prime"):
        return difference <= 1.0 or relative <= 0.20
    if name in {"gamma0", "q", "theta0"}:
        return relative <= 0.25
    if name.startswith("alpha") or name == "dK_dT":
        return relative <= 0.30
    return relative <= 0.20


def _compare(
    record: dict[str, Any], result: Any, thermal: bool, volume_scale: float
) -> tuple[str, list[dict[str, Any]]]:
    published = _published_parameters(record, thermal)
    published_errors = _published_errors(record, thermal)
    comparisons = []
    all_close = True
    all_exact = True
    for name in result.free_parameters:
        source = published[name]
        fitted = float(result.parameters[name])
        fit_error = float(result.standard_errors[name])
        if thermal and name == "rt_eos.V0":
            fitted /= volume_scale
            fit_error /= volume_scale
        source_error = published_errors.get(name)
        combined = None
        within_uncertainty = None
        if source_error is not None and np.isfinite(fit_error):
            combined = math.hypot(source_error, fit_error)
            within_uncertainty = abs(fitted - source) <= 2.0 * combined
        similar = _similar(name, source, fitted)
        uncertainty_parity = bool(within_uncertainty) and similar
        all_close &= similar or bool(within_uncertainty)
        all_exact &= uncertainty_parity
        comparisons.append(
            {
                "parameter": name,
                "published": source,
                "published_error": source_error,
                "refit": fitted,
                "refit_error": fit_error if np.isfinite(fit_error) else None,
                "difference": fitted - source,
                "relative_difference": (
                    abs(fitted - source) / abs(source) if source != 0.0 else None
                ),
                "within_combined_2sigma": within_uncertainty,
                "similar": similar,
            }
        )
    if all_exact and comparisons:
        return "parity", comparisons
    if all_close and comparisons:
        return "similar", comparisons
    return "parity_not_achieved", comparisons


def _combined_fit_dataset(
    document: dict[str, Any], record: dict[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    """Combine explicitly declared P-V inputs after per-dataset transformations."""
    datasets = {item["identifier"]: item for item in document.get("datasets", [])}
    identifiers = list(record["fit_datasets"])
    rows: list[list[float]] = []
    for identifier in identifiers:
        series = _series(document, record, datasets[identifier])
        rows.extend(
            [float(pressure), float(volume)]
            for pressure, volume in zip(series.pressure, series.volume)
        )
    return (
        {
            "identifier": "+".join(identifiers),
            "columns": [
                {
                    "name": "pressure_gpa",
                    "quantity": "pressure",
                    "unit": "GPa",
                    "role": "value",
                },
                {
                    "name": "volume_a3_conventional_cell",
                    "quantity": "conventional_unit_cell_volume",
                    "unit": "angstrom^3/conventional_unit_cell",
                    "role": "value",
                },
            ],
            "rows": rows,
        },
        identifiers,
    )


def _fit_record(
    document: dict[str, Any], record: dict[str, Any], dataset: dict[str, Any]
) -> dict[str, Any]:
    record_id = record["identifier"]
    dataset_identifiers = [dataset["identifier"]]
    if record_id in COMBINED_FIT_DATASET_RECORDS:
        dataset, dataset_identifiers = _combined_fit_dataset(document, record)
    if record_id in INDIRECT_DATA:
        return {
            "status": "not_refittable",
            "reason": INDIRECT_DATA[record_id],
            "dataset_identifiers": dataset_identifiers,
        }
    if dataset["identifier"] == SHEN_PREFIX:
        return {
            "status": "not_refittable",
            "reason": (
                "The workbook contains simultaneous volumes but no pressures, and "
                "the record declares its Cu anchor as reference_model_not_supported."
            ),
            "dataset_identifiers": dataset_identifiers,
        }

    chidester_high_temperature = None
    source_protocol_unweighted = False
    if record_id == "kcl_b2_chidester_2021_bm3_5":
        datasets = {item["identifier"]: item for item in document["datasets"]}
        room_temperature = _series(
            document, record, datasets["kcl_dewaele_2012_table1_compression"]
        )
        chidester_high_temperature = _series(
            document, record, datasets["kcl_chidester_2021_supplemental_pvt"]
        )
        assert chidester_high_temperature.temperature is not None
        reference_temperature = float(record["thermal"]["parameters"]["Tr"])
        dataset_identifiers = [
            "kcl_dewaele_2012_table1_compression",
            "kcl_chidester_2021_supplemental_pvt",
        ]
        series = Series(
            dataset_id="+".join(dataset_identifiers),
            pressure=np.concatenate(
                (room_temperature.pressure, chidester_high_temperature.pressure)
            ),
            volume=np.concatenate(
                (room_temperature.volume, chidester_high_temperature.volume)
            ),
            temperature=np.concatenate(
                (
                    np.full(room_temperature.pressure.shape, reference_temperature),
                    chidester_high_temperature.temperature,
                )
            ),
            # The paper reports an ordinary simultaneous fit and a pressure RMSE,
            # not an effective-variance objective. Its 300 K source rows also lack
            # row-wise sigmas, so the scientifically equivalent fit is unweighted.
            pressure_sigma=None,
            volume_sigma=None,
            temperature_sigma=None,
            pressure_column="pressure_gpa + pt_derived_pressure_gpa",
            volume_column=(
                "volume_a3_per_formula_unit + kcl_molar_volume_cm3_mol"
            ),
            temperature_column="300 K assigned + kcl_temperature_k",
            selection=(
                "all 123 Dewaele 300 K rows and all 155 Chidester high-temperature "
                "rows"
            ),
        )
        source_protocol_unweighted = True
    else:
        series = _series(document, record, dataset)
    if dataset["identifier"] in UNWEIGHTED_DATASETS:
        series.pressure_sigma = None
        series.volume_sigma = None
    material = Material.from_eosmat(document, record_identifiers=[record_id])
    executable = material.eos_records[0].eos
    thermal = isinstance(executable, ThermalEOS)
    static_initial, static_fixed = _static_parameters(record)

    if record_id.startswith("b4c_somayazulu_2023_"):
        assert series.temperature is not None
        heated = series.temperature > 300.0
        series = _masked_series(series, heated, "41 heated rows")

    if record_id == "kcl_walker_2002_bm3_2":
        return _fit_walker_staged(
            record,
            series,
            executable,
            material.eos_records[0].volume_scale,
            material.eos_records[0].reference_temperature,
        )

    if record_id == "aragonite_martinez_1996_bm2_2":
        outcome = _fit_aragonite_staged(record, series)
        outcome["published_rmse_gpa"] = _published_rmse(
            executable,
            material.eos_records[0].volume_scale,
            series,
            material.eos_records[0].reference_temperature,
        )
        return outcome

    # A non-thermal record linked to a P-V-T table represents its reference
    # isotherm; use the rows at the published reference temperature only.
    if not thermal and series.temperature is not None:
        delta = np.abs(
            series.temperature - float(material.eos_records[0].reference_temperature)
        )
        minimum = float(np.min(delta))
        mask = delta <= max(5.0, minimum + 1.0e-9)
        series = _masked_series(series, mask, f"{series.selection}; reference-T rows")
        if series.pressure.size <= len(static_initial):
            return {
                "status": "not_refittable",
                "dataset_identifiers": dataset_identifiers,
                "reason": (
                    f"Only {series.pressure.size} observation(s) lie at the "
                    f"reference temperature for {len(static_initial)} free "
                    "isothermal coefficients; the other rows require a thermal "
                    "relation that this record does not represent."
                ),
            }

    published_rmse = _published_rmse(
        executable,
        material.eos_records[0].volume_scale,
        series,
        material.eos_records[0].reference_temperature,
    )

    if thermal:
        rt_eos = executable.rt_eos
        rt_class = type(rt_eos)
        thermal_initial, thermal_fixed = _thermal_parameters(record)
        if series.temperature is None:
            # A room-temperature P-V series has no information about thermal
            # coefficients. Preserve them while refitting only the reference EOS.
            thermal_fixed.update(thermal_initial)
            thermal_initial = {}
        initial = {f"rt_eos.{key}": value for key, value in static_initial.items()}
        initial.update(thermal_initial)
        fixed = {f"rt_eos.{key}": value for key, value in static_fixed.items()}
        fixed.update(thermal_fixed)
        volume_scale = material.eos_records[0].volume_scale
        if "rt_eos.V0" in initial:
            initial["rt_eos.V0"] *= volume_scale
        if "rt_eos.V0" in fixed:
            fixed["rt_eos.V0"] *= volume_scale
        if not initial:
            raise ValueError("record has no free parameters")
        result = fit_joint_eos(
            type(executable),
            rt_class,
            volume=series.volume * material.eos_records[0].volume_scale,
            temperature=(
                series.temperature
                if series.temperature is not None
                else np.full(
                    series.pressure.shape, material.eos_records[0].reference_temperature
                )
            ),
            pressure=series.pressure,
            initial=initial,
            fixed=fixed,
            configuration=_configuration(record),
            bounds={name: _bounds(name, value) for name, value in initial.items()},
            pressure_sigma=series.pressure_sigma,
            volume_sigma=(
                None
                if series.volume_sigma is None
                else series.volume_sigma * material.eos_records[0].volume_scale
            ),
            temperature_sigma=series.temperature_sigma,
            absolute_sigma=not source_protocol_unweighted,
            max_nfev=5000,
        )
    else:
        eos_class = MODEL_CLASSES[record["eos"]["type"]]
        if not static_initial:
            raise ValueError("record has no free parameters")
        result = fit_rt_eos(
            eos_class,
            volume=series.volume,
            pressure=series.pressure,
            initial=static_initial,
            fixed=static_fixed,
            bounds={
                name: _bounds(name, value) for name, value in static_initial.items()
            },
            pressure_sigma=series.pressure_sigma,
            volume_sigma=series.volume_sigma,
            absolute_sigma=True,
            max_nfev=5000,
        )
    status, comparisons = _compare(
        record, result, thermal, material.eos_records[0].volume_scale
    )
    residuals = np.asarray(result.residuals, dtype=float)
    outcome = {
        "status": status,
        "dataset_identifiers": dataset_identifiers,
        "observations": int(series.pressure.size),
        "selection": series.selection,
        "observed_pressure_range_gpa": [
            float(np.min(series.pressure)),
            float(np.max(series.pressure)),
        ],
        "observed_volume_range": [
            float(np.min(series.volume)),
            float(np.max(series.volume)),
        ],
        "observed_temperature_range_k": (
            None
            if series.temperature is None
            else [
                float(np.min(series.temperature)),
                float(np.max(series.temperature)),
            ]
        ),
        "columns": {
            "pressure": series.pressure_column,
            "volume": series.volume_column,
            "temperature": series.temperature_column,
        },
        "fit_kind": (
            "reference_isotherm_pv"
            if thermal and series.temperature is None
            else "joint_pvt"
            if thermal
            else "isothermal_pv"
        ),
        "objective": (
            "errors_in_variables"
            if series.volume_sigma is not None
            else "pressure_residuals"
        ),
        "absolute_sigma": not source_protocol_unweighted,
        "free_parameters": list(result.free_parameters),
        "parameters": comparisons,
        "published_rmse_gpa": published_rmse,
        "rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "reduced_chi_square": (
            float(result.reduced_chi_square)
            if np.isfinite(result.reduced_chi_square)
            else None
        ),
        "degrees_of_freedom": int(result.degrees_of_freedom),
        "solver_success": bool(result.success),
        "solver_message": str(result.message),
    }
    if record_id in FIT_QUALIFICATIONS:
        outcome["qualification"] = FIT_QUALIFICATIONS[record_id]
    if chidester_high_temperature is not None:
        high_temperature_refit_pressure = result.model.pressure(
            chidester_high_temperature.volume
            * material.eos_records[0].volume_scale,
            chidester_high_temperature.temperature,
        )
        outcome["source_reported_pressure_rmse_gpa"] = 1.6
        outcome["high_temperature_published_rmse_gpa"] = _published_rmse(
            executable,
            material.eos_records[0].volume_scale,
            chidester_high_temperature,
            material.eos_records[0].reference_temperature,
        )
        outcome["high_temperature_refit_rmse_gpa"] = float(
            np.sqrt(
                np.mean(
                    (
                        np.asarray(high_temperature_refit_pressure, dtype=float)
                        - chidester_high_temperature.pressure
                    )
                    ** 2
                )
            )
        )
    return outcome


def _published_rmse(
    executable: Any,
    volume_scale: float,
    series: Series,
    reference_temperature: float,
) -> float:
    if isinstance(executable, ThermalEOS):
        temperature = (
            series.temperature
            if series.temperature is not None
            else np.full(series.pressure.shape, reference_temperature)
        )
        pressure = executable.pressure(series.volume * volume_scale, temperature)
    else:
        pressure = executable.pressure(series.volume)
    return float(
        np.sqrt(np.mean((np.asarray(pressure, dtype=float) - series.pressure) ** 2))
    )


def _fit_aragonite_staged(record: dict[str, Any], series: Series) -> dict[str, Any]:
    """Reproduce Martinez et al.'s staged BM2-isotherm procedure."""
    assert series.temperature is not None
    temperatures = []
    isotherm_v0 = []
    isotherm_k0 = []
    for temperature in sorted(set(series.temperature)):
        if temperature > 973.0:
            continue
        mask = series.temperature == temperature
        if np.count_nonzero(mask) < 3:
            continue
        fit = fit_rt_eos(
            BM2,
            volume=series.volume[mask],
            pressure=series.pressure[mask],
            initial={"V0": 227.5, "K0": 64.81},
        )
        temperatures.append(float(temperature))
        isotherm_v0.append(fit.parameters["V0"])
        isotherm_k0.append(fit.parameters["K0"])
    delta_t = np.asarray(temperatures) - 298.0
    v_coefficients, v_covariance = np.polyfit(
        delta_t, np.asarray(isotherm_v0), 1, cov=True
    )
    k_coefficients, k_covariance = np.polyfit(
        delta_t, np.asarray(isotherm_k0), 1, cov=True
    )
    v_slope, v0 = map(float, v_coefficients)
    d_k_dt, k0 = map(float, k_coefficients)
    alpha0 = v_slope / v0
    alpha_gradient = np.array([1.0 / v0, -v_slope / v0**2])
    alpha_error = float(np.sqrt(alpha_gradient @ v_covariance @ alpha_gradient))
    parameters = {
        "rt_eos.V0": v0,
        "rt_eos.K0": k0,
        "alpha0": alpha0,
        "dK_dT": d_k_dt,
    }

    standard_errors = {
        "rt_eos.V0": float(np.sqrt(v_covariance[1, 1])),
        "rt_eos.K0": float(np.sqrt(k_covariance[1, 1])),
        "alpha0": alpha_error,
        "dK_dT": float(np.sqrt(k_covariance[0, 0])),
    }
    synthetic = SimpleNamespace(
        parameters=parameters,
        standard_errors=standard_errors,
        free_parameters=tuple(parameters),
    )
    status, comparisons = _compare(record, synthetic, True, 1.0)
    fitted_model = ThermalReferenceStateEOS(
        BM2(V0=v0, K0=k0),
        Tr=298.0,
        alpha0=alpha0,
        dK_dT=d_k_dt,
        thermal_expansion_law="constant",
        reference_volume_law="linear_temperature",
    )
    fitted_pressure = fitted_model.pressure(series.volume, series.temperature)
    return {
        "status": status,
        "dataset_identifiers": [series.dataset_id],
        "observations": int(series.pressure.size),
        "selection": "eight staged isotherms from 298 to 973 K",
        "columns": {
            "pressure": series.pressure_column,
            "volume": series.volume_column,
            "temperature": series.temperature_column,
        },
        "fit_kind": "staged_bm2_isotherms",
        "objective": "unweighted pressure residuals followed by linear regressions",
        "absolute_sigma": False,
        "free_parameters": list(parameters),
        "parameters": comparisons,
        "rmse_gpa": float(
            np.sqrt(np.mean((np.asarray(fitted_pressure) - series.pressure) ** 2))
        ),
        "reduced_chi_square": None,
        "solver_success": True,
        "solver_message": "eight isotherm fits and two linear regressions completed",
    }


def _fit_walker_staged(
    record: dict[str, Any],
    series: Series,
    executable: Any,
    volume_scale: float,
    reference_temperature: float,
) -> dict[str, Any]:
    """Reproduce the preferred staged Walker et al. B2-KCl fit."""
    assert series.temperature is not None
    published = record["eos"]["parameters"]
    room_temperature = (
        np.abs(series.temperature - reference_temperature) <= 1.0 + 1.0e-9
    )
    room = _masked_series(
        series,
        room_temperature,
        "eight room-temperature rows at 23-24 degC",
    )
    static_fit = fit_rt_eos(
        BM3,
        volume=room.volume,
        pressure=room.pressure,
        initial={
            "K0": float(published["K0"]),
            "K0_prime": float(published["K0_prime"]),
        },
        fixed={"V0": float(published["V0"])},
        bounds={"K0": (1.0e-6, 200.0), "K0_prime": (0.0, 20.0)},
        absolute_sigma=False,
        max_nfev=5000,
    )
    static_parameters = static_fit.parameters
    thermal_fit = fit_joint_eos(
        type(executable),
        type(executable.rt_eos),
        volume=series.volume * volume_scale,
        temperature=series.temperature,
        pressure=series.pressure,
        initial={"alpha_KT": float(record["thermal"]["parameters"]["alpha_KT"])},
        fixed={
            "rt_eos.V0": float(static_parameters["V0"]) * volume_scale,
            "rt_eos.K0": float(static_parameters["K0"]),
            "rt_eos.K0_prime": float(static_parameters["K0_prime"]),
            "Tr": reference_temperature,
        },
        configuration=_configuration(record),
        bounds={"alpha_KT": (-0.1, 0.1)},
        absolute_sigma=False,
        max_nfev=5000,
    )
    parameters = {
        "rt_eos.K0": float(static_parameters["K0"]),
        "rt_eos.K0_prime": float(static_parameters["K0_prime"]),
        "alpha_KT": float(thermal_fit.parameters["alpha_KT"]),
    }
    standard_errors = {
        "rt_eos.K0": float(static_fit.standard_errors["K0"]),
        "rt_eos.K0_prime": float(static_fit.standard_errors["K0_prime"]),
        "alpha_KT": float(thermal_fit.standard_errors["alpha_KT"]),
    }
    combined = SimpleNamespace(
        parameters=parameters,
        standard_errors=standard_errors,
        free_parameters=tuple(parameters),
    )
    status, comparisons = _compare(record, combined, True, volume_scale)

    simultaneous = fit_joint_eos(
        type(executable),
        type(executable.rt_eos),
        volume=series.volume * volume_scale,
        temperature=series.temperature,
        pressure=series.pressure,
        initial={
            "rt_eos.V0": float(published["V0"]) * volume_scale,
            "rt_eos.K0": float(published["K0"]),
            "rt_eos.K0_prime": float(published["K0_prime"]),
            "alpha_KT": float(record["thermal"]["parameters"]["alpha_KT"]),
        },
        fixed={"Tr": reference_temperature},
        configuration=_configuration(record),
        bounds={
            "rt_eos.V0": _bounds("rt_eos.V0", float(published["V0"]) * volume_scale),
            "rt_eos.K0": _bounds("rt_eos.K0", float(published["K0"])),
            "rt_eos.K0_prime": _bounds("rt_eos.K0_prime", float(published["K0_prime"])),
            "alpha_KT": _bounds(
                "alpha_KT", float(record["thermal"]["parameters"]["alpha_KT"])
            ),
        },
        absolute_sigma=False,
        max_nfev=5000,
    )
    simultaneous_parameters = dict(simultaneous.parameters)
    simultaneous_parameters["rt_eos.V0"] /= volume_scale
    residuals = np.asarray(thermal_fit.residuals, dtype=float)
    outcome = {
        "status": status,
        "dataset_identifiers": [series.dataset_id],
        "observations": int(series.pressure.size),
        "selection": ("eight 23-24 degC rows for K0 and K0'; all 39 rows for alpha_KT"),
        "observed_pressure_range_gpa": [
            float(np.min(series.pressure)),
            float(np.max(series.pressure)),
        ],
        "observed_volume_range": [
            float(np.min(series.volume)),
            float(np.max(series.volume)),
        ],
        "observed_temperature_range_k": [
            float(np.min(series.temperature)),
            float(np.max(series.temperature)),
        ],
        "columns": {
            "pressure": series.pressure_column,
            "volume": series.volume_column,
            "temperature": series.temperature_column,
        },
        "fit_kind": "staged_reference_isotherm_then_thermal_pressure",
        "objective": "unweighted pressure residuals in two source-stated stages",
        "absolute_sigma": False,
        "free_parameters": list(parameters),
        "parameters": comparisons,
        "published_rmse_gpa": _published_rmse(
            executable, volume_scale, series, reference_temperature
        ),
        "rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "reduced_chi_square": None,
        "degrees_of_freedom": int(thermal_fit.degrees_of_freedom),
        "solver_success": bool(static_fit.success and thermal_fit.success),
        "solver_message": "room-temperature BM3 and thermal-pressure stages completed",
        "stages": [
            {
                "name": "reference_isotherm",
                "observations": int(room.pressure.size),
                "temperature_range_k": [
                    float(np.min(room.temperature)),
                    float(np.max(room.temperature)),
                ],
                "fixed_parameters": {"V0": float(static_parameters["V0"])},
                "free_parameters": ["K0", "K0_prime"],
                "parameters": {
                    "K0": float(static_parameters["K0"]),
                    "K0_prime": float(static_parameters["K0_prime"]),
                },
                "rmse_gpa": float(
                    np.sqrt(np.mean(np.asarray(static_fit.residuals, dtype=float) ** 2))
                ),
            },
            {
                "name": "thermal_pressure",
                "observations": int(series.pressure.size),
                "fixed_parameters": {
                    "V0": float(static_parameters["V0"]),
                    "K0": float(static_parameters["K0"]),
                    "K0_prime": float(static_parameters["K0_prime"]),
                },
                "free_parameters": ["alpha_KT"],
                "parameters": {"alpha_KT": float(thermal_fit.parameters["alpha_KT"])},
                "rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
            },
        ],
        "simultaneous_fit_diagnostic": {
            "purpose": (
                "Demonstrates why the original validator disagreed: this is not the "
                "protocol used for the preferred bold Table 3 row."
            ),
            "published_italic_table3_alternative": {
                "V0": 55.25,
                "K0": 14.8,
                "K0_prime": 6.9,
                "alpha0": 0.00018,
                "alpha_KT_derived": 0.002664,
            },
            "peritheos_four_parameter_refit": {
                name: float(value) for name, value in simultaneous_parameters.items()
            },
            "rmse_gpa": float(
                np.sqrt(np.mean(np.asarray(simultaneous.residuals, dtype=float) ** 2))
            ),
        },
    }
    outcome["qualification"] = FIT_QUALIFICATIONS[record["identifier"]]
    return outcome


def _masked_series(series: Series, mask: np.ndarray, selection: str) -> Series:
    def masked(value: np.ndarray | None) -> np.ndarray | None:
        return None if value is None else value[mask]

    return Series(
        dataset_id=series.dataset_id,
        pressure=series.pressure[mask],
        volume=series.volume[mask],
        temperature=masked(series.temperature),
        pressure_sigma=masked(series.pressure_sigma),
        volume_sigma=masked(series.volume_sigma),
        temperature_sigma=masked(series.temperature_sigma),
        pressure_column=series.pressure_column,
        volume_column=series.volume_column,
        temperature_column=series.temperature_column,
        selection=selection,
    )


def validate_all() -> dict[str, Any]:
    results = []
    for material_id in list_material_documents():
        document = get_material_document(material_id)
        datasets = {item["identifier"]: item for item in document.get("datasets", [])}
        for record in document["eos_records"]:
            check = record["scientific_validation"]["primary_data_check"]
            identifiers = list(check.get("dataset_identifiers", ()))
            identifiers += list(check.get("digitized_dataset_identifiers", ()))
            base = {
                "material_identifier": material_id,
                "record_identifier": record["identifier"],
                "label": record["label"],
                "model": record["eos"]["type"],
                "primary_data_status": check["status"],
                "primary_source": record["scientific_validation"][
                    "primary_source_check"
                ].get("access_url"),
                "source_notes": record.get("notes"),
                "primary_data_finding": check.get("finding"),
                "experimental_pressure_range_gpa": record.get(
                    "experimental_pressure_range_gpa"
                ),
                "experimental_temperature_range_k": record.get(
                    "experimental_temperature_range_k"
                ),
                "fixed_parameters": list(record.get("fixed_parameters", ()))
                + list(record["eos"].get("fixed_parameters", ()))
                + list(record.get("thermal", {}).get("fixed_parameters", ())),
            }
            if not identifiers:
                outcome = {
                    "status": "not_refittable",
                    "dataset_identifiers": [],
                    "reason": check["finding"],
                }
            else:
                try:
                    outcome = _fit_record(document, record, datasets[identifiers[0]])
                except Exception as error:  # keep the campaign complete and inspectable
                    outcome = {
                        "status": "refit_failed",
                        "dataset_identifiers": identifiers,
                        "reason": f"{type(error).__name__}: {error}",
                    }
            results.append({**base, **outcome})

    counts = Counter(item["status"] for item in results)
    return {
        "format": "peritheos.primary-eos-refit-validation",
        "format_version": 1,
        "generated_with": "scripts/validate_primary_eos_refits.py",
        "policy": {
            "parity": (
                "Every fitted parameter agrees within combined 2-sigma and also "
                "meets the parameter-specific numerical similarity threshold."
            ),
            "similar": (
                "Every fitted parameter either agrees within combined 2-sigma or "
                "meets the documented parameter-specific similarity threshold."
            ),
            "parity_not_achieved": (
                "At least one fitted parameter is outside both the uncertainty and "
                "similarity criteria."
            ),
            "not_refittable": (
                "The primary source supplies no direct row-level observations, or "
                "the necessary reduction/calibration is not executable."
            ),
        },
        "summary": {"total": len(results), **dict(sorted(counts.items()))},
        "records": results,
    }


def _fmt(value: Any) -> str:
    if value is None:
        return "—"
    return f"{float(value):.6g}"


def _range_text(value: list[float] | None, unit: str) -> str:
    if not value:
        return "not reported"
    return f"{_fmt(value[0])}-{_fmt(value[1])} {unit}"


def _investigation_findings(item: dict[str, Any]) -> list[str]:
    parameters = item.get("parameters", [])
    findings = []
    if item["status"] == "similar":
        uncertainty_missing = [
            p["parameter"] for p in parameters if p["within_combined_2sigma"] is None
        ]
        uncertainty_only = [
            p["parameter"]
            for p in parameters
            if not p["similar"] and p["within_combined_2sigma"] is True
        ]
        numerical_only = [
            p["parameter"]
            for p in parameters
            if p["similar"] and p["within_combined_2sigma"] is False
        ]
        if uncertainty_missing:
            findings.append(
                "Strict uncertainty parity cannot be established because a source "
                "or refit uncertainty is unavailable for "
                + ", ".join(f"`{name}`" for name in uncertainty_missing)
                + ". The point estimates nevertheless meet the numerical criterion."
            )
        if uncertainty_only:
            findings.append(
                "The point estimate exceeds the numerical limit for "
                + ", ".join(f"`{name}`" for name in uncertainty_only)
                + ", but the source and refit two-sigma intervals overlap. This is "
                "evidence of weak coefficient identification rather than a resolved "
                "curve-level disagreement."
            )
        if numerical_only:
            findings.append(
                "The magnitude is similar for "
                + ", ".join(f"`{name}`" for name in numerical_only)
                + ", but the quoted two-sigma intervals do not overlap. Differences "
                "in weighting, rounding, covariance, or the fitted residual variable "
                "remain plausible."
            )
    if item.get("primary_data_status") == "plot_only":
        findings.append(
            "The observations are digitized from a plot. Marker resolution, overlap, "
            "axis calibration, and unavailable source regression weights limit the "
            "strength of any coefficient-level conclusion."
        )
    pressure_range = item.get("observed_pressure_range_gpa")
    v0_free = any(p["parameter"].endswith("V0") for p in parameters)
    if pressure_range and pressure_range[0] > 5.0 and v0_free:
        findings.append(
            f"The lowest checked pressure is {_fmt(pressure_range[0])} GPa while V0 "
            "is free, so the ambient reference volume and correlated elastic "
            "coefficients are extrapolated rather than directly anchored."
        )
    published_rmse = item.get("published_rmse_gpa")
    refit_rmse = item.get("rmse_gpa")
    if (
        published_rmse is not None
        and refit_rmse is not None
        and published_rmse > 2.0 * max(refit_rmse, 1.0e-12)
    ):
        findings.append(
            "The refit reduces pressure RMSE by more than a factor of two. That gap "
            "is too large to attribute only to solver precision and prioritizes a "
            "source row-selection, pressure-scale, weighting, or model-convention "
            "difference."
        )
    boundary = [
        p["parameter"]
        for p in parameters
        if abs(p["refit"]) <= 1.0e-8 and p["parameter"].endswith(("K0_prime", "q"))
    ]
    if boundary:
        findings.append(
            "The optimizer places "
            + ", ".join(f"`{name}`" for name in boundary)
            + " at its physical lower bound. This is direct evidence that the chosen "
            "rows do not identify the unconstrained parameterization."
        )
    if item["record_identifier"] in INVESTIGATION_NOTES:
        findings.append(INVESTIGATION_NOTES[item["record_identifier"]])
    if not findings:
        findings.append(
            "No single failure mechanism is established by the available metadata. "
            "The next reproducible step is to recover the publication's exact row "
            "mask, fixed coefficients, residual definition, and covariance treatment."
        )
    return findings


def render_markdown(ledger: dict[str, Any]) -> str:
    summary = ledger["summary"]
    lines = [
        "# Primary EOS refit validation",
        "",
        "This is the record-by-record numerical reproduction ledger for the bundled ",
        "material library. It is generated by ",
        "`scripts/validate_primary_eos_refits.py` from the checked-in primary-data ",
        "transcriptions. A source parameterization remains authoritative; these refits ",
        "are independent diagnostics and never overwrite a library record.",
        "",
        "## Outcome",
        "",
        f"The campaign covers all **{summary['total']}** EOS records. "
        f"**{summary.get('parity', 0)}** achieve uncertainty parity, "
        f"**{summary.get('similar', 0)}** are numerically similar, "
        f"**[{summary.get('parity_not_achieved', 0)}](#parity-not-achieved)** do not "
        "achieve parity, "
        f"**{summary.get('not_refittable', 0)}** cannot be directly refitted, and "
        f"**{summary.get('refit_failed', 0)}** attempts failed before comparison.",
        "",
        "`parity` means all free coefficients agree within two combined standard ",
        "uncertainties and also meet the numerical similarity limits. This prevents an ",
        "unidentifiable extrapolation with an enormous fitted error from being labeled ",
        "parity. Where a source or refit uncertainty is unavailable, `similar` ",
        "uses explicit limits: 5% for V0, 15% for K0, 1.0 absolute or 20% for K0', ",
        "25% for gamma0/q/theta0, 30% for thermal-expansion or dK/dT terms, and 20% ",
        "for other coefficients. These broad limits identify broadly reproducible ",
        "published reductions; they are not statistical confidence statements.",
        "",
        "Fits use the equation and fixed coefficients declared by each record. Published ",
        "row-wise uncertainties are used only when complete and positive; otherwise the ",
        "fit is unweighted. P-V-T records are joint fits except where the publication ",
        "declares a staged protocol, which is reproduced explicitly. Isothermal records linked to a ",
        "P-V-T table use the rows nearest the declared reference temperature. Digitized ",
        "plot data are labeled `plot_only` and should be interpreted less strictly.",
        "A record-specific qualification is shown where the checked-in rows cover ",
        "only a subset of the source fit or source-fixed coefficients must be preserved.",
        "Regenerate the ledger with `uv run python scripts/validate_primary_eos_refits.py`; ",
        "use `--check` in continuous integration to detect stale generated files.",
        "",
        "## Records",
        "",
        "| Record | Data | n | Published → Peritheos refit | Curve/refit RMSE (GPa) | Outcome |",
        "|---|---|---:|---|---:|---|",
    ]
    for item in ledger["records"]:
        params = []
        for parameter in item.get("parameters", []):
            params.append(
                f"`{parameter['parameter']}` {_fmt(parameter['published'])} → "
                f"{_fmt(parameter['refit'])}"
            )
        data = ", ".join(item["dataset_identifiers"]) or item["primary_data_status"]
        reason = item.get("reason") or item.get("qualification")
        outcome = item["status"]
        if item["status"] in {"similar", "parity_not_achieved", "refit_failed"}:
            anchor = f"investigation-{item['record_identifier']}"
            outcome = f"[{outcome}](#{anchor})"
        if reason:
            outcome += f" — {reason}"
        outcome = outcome.replace("|", "\\|").replace("\n", " ")
        record_label = f"`{item['record_identifier']}`"
        if item.get("primary_source"):
            record_label = f"[{record_label}]({item['primary_source']})"
        rmse = f"{_fmt(item.get('published_rmse_gpa'))}/{_fmt(item.get('rmse_gpa'))}"
        lines.append(
            f"| {record_label} | `{data}` | "
            f"{item.get('observations', '—')} | {'; '.join(params) or '—'} | "
            f"{rmse} | {outcome} |"
        )

    unsuccessful = [
        item
        for item in ledger["records"]
        if item["status"] in {"parity_not_achieved", "refit_failed"}
    ]
    unavailable = [
        item for item in ledger["records"] if item["status"] == "not_refittable"
    ]
    lines.extend(["", "## Parity not achieved", ""])
    if unsuccessful:
        for item in unsuccessful:
            mismatches = [
                parameter
                for parameter in item.get("parameters", [])
                if not parameter["similar"]
            ]
            detail = "; ".join(
                f"{parameter['parameter']} {_fmt(parameter['published'])} → "
                f"{_fmt(parameter['refit'])}"
                for parameter in mismatches
            )
            reason = item.get("reason") or (
                f"outside similarity limits ({detail})" if detail else item["status"]
            )
            anchor = f"investigation-{item['record_identifier']}"
            lines.append(f"- [`{item['record_identifier']}`](#{anchor}): {reason}")
    else:
        lines.append("No completed fit fell outside the parity/similarity criteria.")

    investigated = [
        item
        for item in ledger["records"]
        if item["status"] in {"similar", "parity_not_achieved", "refit_failed"}
    ]
    lines.extend(
        [
            "",
            "## Corrections found during investigation",
            "",
            "The first campaign incorrectly reported both Fei neon records as large ",
            "V0 failures. Fei et al. fixed V0 and fitted all room-temperature data ",
            "from two earlier studies plus their own observations; freeing V0 on only ",
            "the 13 new high-pressure points created an unconstrained extrapolation. ",
            "The fit inputs are now identified explicitly as Hemley et al. (1989, ref. ",
            "45), Finger et al. (1981, ref. 47), and Fei et al.'s new observations. ",
            "For the 21 bundled Hemley rows, the printed Mao-scale pressures are ",
            "recalculated to the Dewaele et al. (2004) ruby scale as stated in Fei's ",
            "Figure 5 caption, then combined with the 13 digitized observations from ",
            "this study. With V0 fixed and thermal terms held fixed for the 300 K-only ",
            "data, the 34-point partial BM3 fit gives K0=1.4775 GPa and K0'=7.8503; ",
            "the Vinet fit gives K0=1.1439 GPa and K0'=8.2581. Both retain parity. ",
            "Finger's low-pressure rows are still missing, so these results are ",
            "explicitly classified as partial reproductions rather than exact reruns ",
            "of Fei's complete regression.",
            "",
            "The first campaign also applied the B2-CaO fitted V0 to a table column ",
            "normalized to the ambient B1 volume. Reconstructing B2 volumes from the ",
            "tabulated cubic lattice parameter removes that false failure.",
            "",
            "The Walker B2-KCl failure was a fit-protocol mismatch. The preferred ",
            "bold Table 3 result is staged, while the first campaign fitted all four ",
            "coefficients simultaneously and used the printed row uncertainties. ",
            "Following the source's unweighted staged protocol recovers K0, K0', and ",
            "alpha_KT within 0.6%; an unweighted simultaneous diagnostic also recovers ",
            "the paper's separate italic Table 3 row.",
            "",
            "The Tateno B2-KCl failure combined a source-version error with a data-row ",
            "alignment error. The final article replaces the accepted manuscript's ",
            "gamma0=0.58 and q=0.9 with gamma0=2.3 and q=0.8 and uses the ",
            "integrated-Gruneisen Debye law. The accepted-manuscript table split the ",
            "Pt and KCl halves in different orders, so runs 3 and 4 were mismatched. ",
            "Using the correctly aligned official MSA Supplemental Table S1 workbook ",
            "recovers all four fitted coefficients within combined two-sigma uncertainty.",
            "",
            "The Chidester B2-KCl boundary solution was another fit-scope error. ",
            "The paper fitted its 155 new high-temperature observations together ",
            "with all 123 Dewaele et al. (2012) room-temperature B2 observations, ",
            "whereas the first campaign fitted only the new rows and applied their ",
            "row-wise uncertainties. The corrected unweighted 278-row simultaneous ",
            "fit also uses the thermodynamically integrated Gruneisen ",
            "Debye-temperature relation. It recovers all five coefficients within ",
            "combined two-sigma uncertainty, and its 1.582 GPa high-temperature ",
            "RMSE reproduces the paper's reported 1.6 GPa.",
            "",
            "The Campbell-Heinz CsCl failure was a material-assignment error. ",
            "Table 1 prints 24 RbCl rows followed by 13 CsCl rows, but the first ",
            "campaign had attached only the first 21 RbCl rows to CsCl. Restoring ",
            "the correct CsCl block and adding all nine corrected Yagi Table 1 ",
            "ratios changes the refit from K0=2.4298 GPa and K0'=17.1164 to ",
            "K0=17.6967 GPa and K0'=5.2026, which is compatible with the ",
            "combined Campbell + Yagi values. The complete 24-row RbCl ",
            "dataset is now a separate material record and refits to K0=17.8808 ",
            "GPa and K0'=5.2382, reproducing the printed 17.9(10) and 5.23(29).",
            "",
            "## Detailed non-parity investigations",
            "",
            f"The following **{len(investigated)}** sections cover every completed "
            "refit that does not meet the strict `parity` definition. `similar` means ",
            "the difference is numerically acceptable or covered by combined ",
            "uncertainty; `parity_not_achieved` means at least one coefficient is ",
            "outside both tests. Causes described as possible remain hypotheses until ",
            "the missing source fit detail is recovered.",
        ]
    )
    for item in investigated:
        anchor = f"investigation-{item['record_identifier']}"
        lines.extend(
            [
                "",
                f'<a id="{anchor}"></a>',
                "",
                f"### `{item['record_identifier']}`",
                "",
                f"**Classification:** `{item['status']}`. **Model:** `{item['model']}`. "
                f"**Data:** `{', '.join(item['dataset_identifiers'])}` with "
                f"{item.get('observations', '—')} selected observations.",
                "",
                "| Parameter | Published | Refit ± 1σ | Relative difference | "
                "Within combined 2σ | Numerical limit |",
                "|---|---:|---:|---:|:---:|:---:|",
            ]
        )
        for parameter in item.get("parameters", []):
            relative = parameter.get("relative_difference")
            relative_text = "—" if relative is None else f"{100.0 * relative:.2f}%"
            within = parameter.get("within_combined_2sigma")
            within_text = "—" if within is None else "yes" if within else "no"
            error = parameter.get("refit_error")
            refit_text = _fmt(parameter["refit"])
            if error is not None:
                refit_text += f" ± {_fmt(error)}"
            lines.append(
                f"| `{parameter['parameter']}` | {_fmt(parameter['published'])} | "
                f"{refit_text} | {relative_text} | {within_text} | "
                f"{'yes' if parameter['similar'] else 'no'} |"
            )
        lines.extend(
            [
                "",
                "**Fit diagnostics.** "
                f"Observed pressure range: {_range_text(item.get('observed_pressure_range_gpa'), 'GPa')}; "
                f"source-declared range: {_range_text(item.get('experimental_pressure_range_gpa'), 'GPa')}; "
                f"fit kind: `{item.get('fit_kind', '—')}`; objective: "
                f"`{item.get('objective', '—')}`; published/refit pressure RMSE: "
                f"{_fmt(item.get('published_rmse_gpa'))}/{_fmt(item.get('rmse_gpa'))} GPa; "
                f"reduced chi-square: {_fmt(item.get('reduced_chi_square'))}; "
                f"free parameters: `{', '.join(item.get('free_parameters', ()))}`; "
                f"source-fixed parameters: `{', '.join(item.get('fixed_parameters', ())) or 'none'}`.",
                "",
                f"**Source/data scope.** {item.get('primary_data_finding') or 'No additional source-scope note is registered.'}",
            ]
        )
        if item.get("source_notes"):
            lines.extend(
                ["", f"**Registered source-fit note.** {item['source_notes']}"]
            )
        lines.extend(["", "**Assessment and likely origin.**"])
        for finding in _investigation_findings(item):
            lines.append(f"- {finding}")

    lines.extend(["", "## Direct refit unavailable", ""])
    for item in unavailable:
        lines.append(f"- `{item['record_identifier']}`: {item['reason']}")
    lines.extend(
        [
            "",
            "The complete machine-readable diagnostics, including selected columns, ",
            "fit objective, coefficient errors, relative differences, solver status, ",
            "and reduced chi-square, are in ",
            "[`docs/data/primary-eos-refits.json`](data/primary-eos-refits.json).",
            "",
        ]
    )
    return "\n".join(line.rstrip() for line in lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if outputs differ")
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--markdown", type=Path, default=DEFAULT_MARKDOWN)
    args = parser.parse_args()
    ledger = validate_all()
    json_text = json.dumps(ledger, indent=2, sort_keys=True, allow_nan=False) + "\n"
    markdown_text = render_markdown(ledger)
    if args.check:
        current_json = (
            args.json.read_text(encoding="utf-8") if args.json.exists() else ""
        )
        current_markdown = (
            args.markdown.read_text(encoding="utf-8") if args.markdown.exists() else ""
        )
        if current_json != json_text or current_markdown != markdown_text:
            print("primary EOS refit documentation is stale", file=sys.stderr)
            return 1
        return 0
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json_text, encoding="utf-8")
    args.markdown.write_text(markdown_text, encoding="utf-8")
    print(json.dumps(ledger["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
