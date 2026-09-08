"""Audit the recoverable inputs behind Dorogokupets et al. (2015).

The paper describes a joint 14-parameter thermoelastic optimization and does
not state a point-weighting rule.  This script therefore treats unweighted
least squares as the primary reproducible interpretation and keeps an
uncertainty-weighted fit only as a sensitivity check.  The pressure coordinates
are reconstructed from the printed Au or MgO calibrant observations wherever
the primary tables retain them.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from peritheos import Material, get_material_document
from peritheos.materials import AU_SOKOLOVA_2013, MGO_SOKOLOVA_2013

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "peritheos/data/datasets"
RECONSTRUCTION = DATASETS / "dorogokupets-2015-298k-literature-reconstruction.csv"
REPORT = ROOT / "docs/data/dorogokupets-2015-298k-refit.json"

KATSURA = DATASETS / "bridgmanite-katsura-2009-corrected-table1-pvt.csv"
TANGE = DATASETS / "bridgmanite-tange-2012-table1-pvt.csv"
REYNARD = DATASETS / "akimotoite-reynard-1996-table1-compression.csv"
GUIGNOT = DATASETS / "mgsio3-post-perovskite-guignot-2007-table1-300k-compression.csv"
WANG = DATASETS / "akimotoite-wang-2004-tables1-2-pvt.csv"
KOMABAYASHI = DATASETS / "mgsio3-post-perovskite-komabayashi-2008-table1-pvt.csv"
ZHOU = DATASETS / "akimotoite-zhou-2014-table1-elasticity.csv"

MGO_SOKOLOVA_V0_A3 = 74.71096452981054
KATSURA_BRIDGMANITE_V0_A3 = 162.3855988669


@dataclass(frozen=True)
class Observation:
    phase: str
    source: str
    source_doi: str
    source_order: str
    temperature_k: float
    volume_a3: float
    volume_uncertainty_a3: float | None
    pressure_gpa: float
    pressure_uncertainty_gpa: float | None
    pressure_coordinate: str
    selection_note: str


PUBLISHED = {
    "bridgmanite": dict(V0=162.4007207696022, K0=252.0, K0_prime=4.38),
    "akimotoite": dict(V0=262.53122652018516, K0=215.3, K0_prime=4.91),
    "post_perovskite": dict(V0=160.74018170242837, K0=253.7, K0_prime=4.03),
}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _optional(value: str | None) -> float | None:
    return None if value in (None, "") else float(value)


def _mgo_pressure(volume_a3: float, temperature_k: float) -> float:
    return float(MGO_SOKOLOVA_2013.pressure(volume_a3, temperature_k))


def _wang_gold_pressure(original_pressure_gpa: float, temperature_k: float) -> float:
    """Map Wang's Anderson-Au pressure to the 2015 Sokolova-Au coordinate."""
    material = Material.from_eosmat(
        get_material_document("gold"),
        record_identifiers=["gold_anderson_1989_bm3_1"],
    )
    anderson = material.eos_records[0]
    gold_volume = anderson.volume(original_pressure_gpa, temperature_k)
    return float(AU_SOKOLOVA_2013.pressure(gold_volume, temperature_k))


def build_observations() -> list[Observation]:
    observations: list[Observation] = []

    # Dorogokupets et al. identify Katsura and Tange as the Prv fit inputs.
    # Katsura reports only MgO V/V0.  Applying that rounded ratio to the
    # Sokolova reference volume is an explicit scale-mapping approximation.
    for row in _rows(KATSURA):
        temperature = float(row["temperature_k"])
        if temperature > 308.0:
            continue
        mgo_ratio = float(row["mgo_v_v0"])
        observations.append(
            Observation(
                phase="bridgmanite",
                source="Katsura et al. (2009), corrected Table 1",
                source_doi="10.1029/2009GL039318",
                source_order=row["source_order"],
                temperature_k=temperature,
                volume_a3=(float(row["bridgmanite_v_v0"]) * KATSURA_BRIDGMANITE_V0_A3),
                volume_uncertainty_a3=(
                    float(row["bridgmanite_v_v0_uncertainty"])
                    * KATSURA_BRIDGMANITE_V0_A3
                ),
                pressure_gpa=_mgo_pressure(mgo_ratio * MGO_SOKOLOVA_V0_A3, temperature),
                pressure_uncertainty_gpa=float(row["pressure_uncertainty_gpa"]),
                pressure_coordinate="Sokolova-2013 MgO from rounded source V/V0",
                selection_note=(
                    "298 K slice diagnostic; exact 2015 recalculation cannot be "
                    "recovered because absolute MgO volumes are not printed"
                ),
            )
        )

    for row in _rows(TANGE):
        temperature = float(row["temperature_k"])
        if temperature != 300.0:
            continue
        observations.append(
            Observation(
                phase="bridgmanite",
                source="Tange et al. (2012), Table 1",
                source_doi="10.1029/2011JB008988",
                source_order=row["run"],
                temperature_k=temperature,
                volume_a3=float(row["unit_cell_volume_a3"]),
                volume_uncertainty_a3=float(row["unit_cell_volume_sd_a3"]),
                pressure_gpa=_mgo_pressure(
                    float(row["mgo_unit_cell_volume_a3"]), temperature
                ),
                pressure_uncertainty_gpa=float(row["pressure_bm3_sd_gpa"]),
                pressure_coordinate="Sokolova-2013 MgO from measured cell volume",
                selection_note="Source explicitly selected by Dorogokupets et al.",
            )
        )

    # The 2015 Figure 6 caption says Reynard pressures are retained on their
    # original coordinate.  The primary source publishes two alternatives, so
    # both are kept as separate, non-independent diagnostics.
    for source_order, row in enumerate(_rows(REYNARD), start=1):
        observations.append(
            Observation(
                phase="akimotoite_ruby",
                source="Reynard et al. (1996), Table 1",
                source_doi="10.2138/am-1996-1-206",
                source_order=str(source_order),
                temperature_k=298.0,
                volume_a3=float(row["unit_cell_volume_a3"]),
                volume_uncertainty_a3=float(row["unit_cell_volume_sd_a3"]),
                pressure_gpa=float(row["ruby_pressure_gpa"]),
                pressure_uncertainty_gpa=_optional(
                    row["ruby_linewidth_uncertainty_gpa"]
                ),
                pressure_coordinate="published ruby pressure",
                selection_note=(
                    "Alternative assignment of the same experiment; exact ruby "
                    "calibration is not stated"
                ),
            )
        )
        ice_pressure = _optional(row["ice_vii_pressure_gpa"])
        if ice_pressure is not None or float(row["ruby_pressure_gpa"]) == 0.0:
            observations.append(
                Observation(
                    phase="akimotoite_ice_vii",
                    source="Reynard et al. (1996), Table 1/Figure 3f",
                    source_doi="10.2138/am-1996-1-206",
                    source_order=str(source_order),
                    temperature_k=298.0,
                    volume_a3=float(row["unit_cell_volume_a3"]),
                    volume_uncertainty_a3=float(row["unit_cell_volume_sd_a3"]),
                    pressure_gpa=0.0 if ice_pressure is None else ice_pressure,
                    pressure_uncertainty_gpa=_optional(row["ice_vii_max_gradient_gpa"]),
                    pressure_coordinate="published ice-VII pressure",
                    selection_note=(
                        "Preferred Reynard scale; bracketed values are maximum "
                        "gradients, not stated one-sigma errors"
                    ),
                )
            )

    # Figure 6 assigns Wang run T0133 (Au marker, 300-773 K) to the
    # authors' Au scale and leaves the other Wang observations on their
    # original coordinate.  The paper prints 298 K; Dorogokupets label the
    # corresponding isotherm 300 K.
    for row in _rows(WANG):
        temperature = float(row["temperature_k"])
        if temperature != 298.0:
            continue
        original_pressure = float(row["pressure_original_gpa"])
        if row["run"] == "T0133":
            pressure = _wang_gold_pressure(original_pressure, temperature)
            pressure_coordinate = (
                "Sokolova-2013 Au from Anderson-1989 pressure inversion"
            )
            note = (
                "Figure 6 recalculation: inferred Au volume is the unique state "
                "that reproduces Wang's printed Anderson pressure"
            )
        else:
            pressure = original_pressure
            pressure_coordinate = "published Decker-1971 NaCl pressure"
            note = "Figure 6 retains the non-Au Wang observations as original"
        observations.append(
            Observation(
                phase="akimotoite_wang",
                source=f"Wang et al. (2004), Table {row['table']}",
                source_doi="10.1016/j.pepi.2003.08.007",
                source_order=f"{row['run']}-{row['source_order']}",
                temperature_k=temperature,
                volume_a3=float(row["unit_cell_volume_a3"]),
                volume_uncertainty_a3=float(row["unit_cell_volume_uncertainty_a3"]),
                pressure_gpa=pressure,
                pressure_uncertainty_gpa=None,
                pressure_coordinate=pressure_coordinate,
                selection_note=note,
            )
        )

    for row in _rows(GUIGNOT):
        temperature = float(row["temperature_k"])
        observations.append(
            Observation(
                phase="post_perovskite",
                source="Guignot et al. (2007), Table 1",
                source_doi="10.1016/j.epsl.2007.01.025",
                source_order=row["source_order"],
                temperature_k=temperature,
                volume_a3=float(row["mgsio3_unit_cell_volume_a3"]),
                volume_uncertainty_a3=float(
                    row["mgsio3_unit_cell_volume_uncertainty_a3"]
                ),
                pressure_gpa=_mgo_pressure(
                    float(row["mgo_lattice_a_angstrom"]) ** 3, temperature
                ),
                pressure_uncertainty_gpa=float(
                    row["pressure_original_uncertainty_gpa"]
                ),
                pressure_coordinate="Sokolova-2013 MgO from measured lattice a",
                selection_note=("Source explicitly selected by Dorogokupets et al."),
            )
        )

    for row in _rows(KOMABAYASHI):
        temperature = float(row["temperature_k"])
        if temperature != 300.0 or row["ppv_fit_included"] != "true":
            continue
        observations.append(
            Observation(
                phase="post_perovskite",
                source="Komabayashi et al. (2008), Table 1",
                source_doi="10.1016/j.epsl.2007.10.036",
                source_order=row["source_order"],
                temperature_k=temperature,
                volume_a3=float(row["ppv_volume_a3"]),
                volume_uncertainty_a3=float(row["ppv_volume_uncertainty_a3"]),
                pressure_gpa=_mgo_pressure(
                    float(row["mgo_a_angstrom"]) ** 3, temperature
                ),
                pressure_uncertainty_gpa=_optional(row["pressure_mgo_uncertainty_gpa"]),
                pressure_coordinate="Sokolova-2013 MgO from measured lattice a",
                selection_note="Source explicitly selected by Dorogokupets et al.",
            )
        )

    return observations


def rydberg_stacey(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    x = (volume / v0) ** (1.0 / 3.0)
    k = 5.0
    eta = 1.5 * k0_prime - k + 0.5
    return 3.0 * k0 * x ** (-k) * (1.0 - x) * np.exp(eta * (1.0 - x))


def _effective_sigma(
    observations: list[Observation], published: dict[str, float]
) -> np.ndarray:
    volume = np.asarray([item.volume_a3 for item in observations])
    step = np.maximum(1.0e-5 * volume, 1.0e-5)
    derivative = (
        rydberg_stacey(
            volume + step,
            published["V0"],
            published["K0"],
            published["K0_prime"],
        )
        - rydberg_stacey(
            volume - step,
            published["V0"],
            published["K0"],
            published["K0_prime"],
        )
    ) / (2.0 * step)
    pressure_sigma = np.asarray(
        [item.pressure_uncertainty_gpa or math.nan for item in observations]
    )
    volume_sigma = np.asarray(
        [item.volume_uncertainty_a3 or 0.0 for item in observations]
    )
    combined = np.sqrt(pressure_sigma**2 + (derivative * volume_sigma) ** 2)
    fallback = float(np.nanmedian(combined)) if np.any(np.isfinite(combined)) else 1.0
    return np.where(np.isfinite(combined) & (combined > 0.0), combined, fallback)


def _fit(
    observations: list[Observation], published: dict[str, float], weighted: bool
) -> dict[str, Any]:
    volume = np.asarray([item.volume_a3 for item in observations])
    pressure = np.asarray([item.pressure_gpa for item in observations])
    sigma = _effective_sigma(observations, published)

    def residual(parameters: np.ndarray) -> np.ndarray:
        calculated = rydberg_stacey(
            volume, published["V0"], parameters[0], parameters[1]
        )
        delta = calculated - pressure
        return delta / sigma if weighted else delta

    result = least_squares(
        residual,
        np.asarray([published["K0"], published["K0_prime"]]),
        bounds=([1.0, 1.0], [1000.0, 12.0]),
    )
    calculated = rydberg_stacey(
        volume, published["V0"], float(result.x[0]), float(result.x[1])
    )
    published_calculated = rydberg_stacey(
        volume,
        published["V0"],
        published["K0"],
        published["K0_prime"],
    )
    return {
        "weighting": (
            "source-reported pressure error/spread proxy combined in quadrature "
            "with published-curve volume-uncertainty projection; not the "
            "unpublished Dorogokupets objective weights"
            if weighted
            else "unweighted pressure residuals"
        ),
        "K0_gpa": float(result.x[0]),
        "K0_prime": float(result.x[1]),
        "rmse_gpa": float(np.sqrt(np.mean((calculated - pressure) ** 2))),
        "published_curve_rmse_gpa": float(
            np.sqrt(np.mean((published_calculated - pressure) ** 2))
        ),
        "maximum_absolute_residual_gpa": float(np.max(np.abs(calculated - pressure))),
        "success": bool(result.success),
    }


def reproduce() -> dict[str, Any]:
    observations = build_observations()
    wang = [o for o in observations if o.phase == "akimotoite_wang"]
    reynard_ruby = [o for o in observations if o.phase == "akimotoite_ruby"]
    reynard_ice = [o for o in observations if o.phase == "akimotoite_ice_vii"]
    groups = {
        "bridgmanite": [o for o in observations if o.phase == "bridgmanite"],
        "akimotoite_ruby": reynard_ruby,
        "akimotoite_ice_vii": reynard_ice,
        "akimotoite_wang": wang,
        "akimotoite_wang_plus_reynard_ruby": wang + reynard_ruby,
        "akimotoite_wang_plus_reynard_ice_vii": wang + reynard_ice,
        "post_perovskite": [o for o in observations if o.phase == "post_perovskite"],
    }
    output: dict[str, Any] = {
        "method": (
            "298 K slice diagnostics with V0 and k=5 fixed to Dorogokupets "
            "et al. (2015); K0 and K0_prime optimized"
        ),
        "limitations": [
            "Dorogokupets et al. state no point-weighting rule or covariance; unweighted least squares is the primary reconstruction.",
            "The joint thermal/acoustic objective is not reduced to a 298 K P-V fit.",
            "The Wang, Zhou, and Komabayashi source tables are complete, but the 14-parameter joint optimizer and its phase-specific fixed/free mask are not published.",
            "Katsura MgO V/V0 is rounded and lacks the absolute cell volume needed for exact rescaling.",
        ],
        "source_tables": {
            "Wang_2004": {"rows": len(_rows(WANG)), "sha256": sha256(WANG)},
            "Komabayashi_2008": {
                "rows": len(_rows(KOMABAYASHI)),
                "usable_post_perovskite_rows": sum(
                    row["ppv_fit_included"] == "true" for row in _rows(KOMABAYASHI)
                ),
                "sha256": sha256(KOMABAYASHI),
            },
            "Zhou_2014": {
                "rows": len(_rows(ZHOU)),
                "rows_with_KS": sum(bool(row["ks_gpa"]) for row in _rows(ZHOU)),
                "sha256": sha256(ZHOU),
            },
        },
        "zhou_2014_elasticity": zhou_elasticity_diagnostics(),
        "fits": {},
    }
    for name, rows in groups.items():
        published_key = "akimotoite" if name.startswith("akimotoite") else name
        published = PUBLISHED[published_key]
        output["fits"][name] = {
            "observations": len(rows),
            "pressure_range_gpa": [
                min(item.pressure_gpa for item in rows),
                max(item.pressure_gpa for item in rows),
            ],
            "published": published,
            "unweighted": _fit(rows, published, False),
            "reported_uncertainty_diagnostic": _fit(rows, published, True),
        }
    return output


def zhou_elasticity_diagnostics() -> dict[str, Any]:
    """Reproduce Zhou's linear KS(P,T) regression with and without sigma weights."""
    rows = [row for row in _rows(ZHOU) if row["ks_gpa"]]
    design = np.asarray(
        [
            [1.0, float(row["pressure_gpa"]), float(row["temperature_k"]) - 300.0]
            for row in rows
        ]
    )
    values = np.asarray([float(row["ks_gpa"]) for row in rows])
    sigma = np.asarray([float(row["ks_uncertainty_gpa"]) for row in rows])

    def regress(weights: np.ndarray) -> dict[str, float]:
        fitted = np.linalg.lstsq(
            design * weights[:, None], values * weights, rcond=None
        )[0]
        residuals = design @ fitted - values
        return {
            "K0S_gpa": float(fitted[0]),
            "dKS_dP": float(fitted[1]),
            "dKS_dT_gpa_per_k": float(fitted[2]),
            "rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        }

    return {
        "observations": len(rows),
        "published_linear": {
            "K0S_gpa": 219.4,
            "dKS_dP": 4.62,
            "dKS_dT_gpa_per_k": -0.0228,
        },
        "unweighted": regress(np.ones(len(rows))),
        "reported_KS_uncertainty_sensitivity": regress(1.0 / sigma),
        "interpretation": (
            "Both regressions recover the published rounded coefficients closely; "
            "the source does not identify uncertainty weighting as its objective."
        ),
    }


def write_reconstruction(path: Path = RECONSTRUCTION) -> None:
    fieldnames = list(Observation.__dataclass_fields__)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for item in build_observations():
            writer.writerow(item.__dict__)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-dataset", action="store_true")
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()
    if args.write_dataset:
        write_reconstruction()
    result = reproduce()
    if RECONSTRUCTION.exists():
        result["reconstruction_sha256"] = sha256(RECONSTRUCTION)
    if args.write_report:
        REPORT.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
