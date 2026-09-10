#!/usr/bin/env python3
"""Audit the Noguchi et al. (2013) CaSiO3 Table 1 source fit.

The 54-row table is not redistributed because the subscription article does
not grant a reusable data license. This script therefore has two modes:

* without arguments, verify and report the shipped non-row-level audit artifact;
* with ``--source-table``, rerun all four staged Table 2 fits from a lawful local
  transcription and compare its canonical digest with the audited transcription.

The local CSV must contain ``source_order``, ``pressure_gpa`` (or
``pressure_fei_gpa``), ``pressure_holmes_gpa``, ``temperature_k``, and
``volume_a3``. ``fit_included`` is checked when present. ``pt_lattice_a`` (or
``a_pt_a``) is optional; when supplied it also reruns both Pt-scale checks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from peritheos import get_material_document
from peritheos.eos.rt import BM2, BM3
from peritheos.eos.thermal import LogVolumeThermalPressure, MieGruneisenDebye
from peritheos.fitting import FitResult, fit_rt_eos, fit_thermal_eos
from peritheos.materials import Material
from peritheos.units import cell_volume_to_molar_volume

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_PATH = ROOT / "docs" / "data" / "noguchi-2013-casio3-refit.json"
RECORD_ID = "ca_perovskite_noguchi_2013_bm2_mgd_1"
HOLMES_RECORD_ID = "platinum_holmes_1989_vinet_1"
EXPECTED_EXCLUDED_ORDERS = (12, 14, 20)
STATIC_SOURCE_ORDERS = tuple(
    order for order in range(9, 21) if order not in EXPECTED_EXCLUDED_ORDERS
)
THERMAL_SOURCE_ORDERS = tuple((*range(1, 9), *range(21, 55)))
CELL_TO_MOLAR = float(cell_volume_to_molar_volume(1.0, 1.0))


def isothermal_bulk_modulus(record, volume: float, temperature: float) -> float:
    """Evaluate ``-V(dP/dV)_T`` in the record's public cell-volume unit."""
    step = volume * 1.0e-6
    derivative = (
        record.pressure(volume + step, temperature)
        - record.pressure(volume - step, temperature)
    ) / (2.0 * step)
    return -volume * derivative


def load_artifact(path: Path = ARTIFACT_PATH) -> dict[str, Any]:
    """Load the row-free committed fit summary."""
    return json.loads(path.read_text(encoding="utf-8"))


def _field(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return value
    raise ValueError(f"local source table is missing one of {names!r}")


def _boolean(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"cannot parse boolean value {value!r}")


def load_source_table(path: Path) -> list[dict[str, float | int | bool | None]]:
    """Read and validate a private local transcription without retaining it."""
    rows: list[dict[str, float | int | bool | None]] = []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for raw in csv.DictReader(stream):
            order = int(_field(raw, "source_order"))
            expected_included = order not in EXPECTED_EXCLUDED_ORDERS
            included_text = raw.get("fit_included")
            included = (
                expected_included
                if included_text in (None, "")
                else _boolean(included_text)
            )
            if included != expected_included:
                raise ValueError(
                    f"source row {order} fit_included disagrees with Table 1 footnote e"
                )
            pt_text = raw.get("pt_lattice_a") or raw.get("a_pt_a")
            rows.append(
                {
                    "source_order": order,
                    "pressure_fei_gpa": float(
                        _field(raw, "pressure_fei_gpa", "pressure_gpa")
                    ),
                    "pressure_holmes_gpa": float(_field(raw, "pressure_holmes_gpa")),
                    "temperature_k": float(_field(raw, "temperature_k")),
                    "volume_a3": float(_field(raw, "volume_a3")),
                    "pt_lattice_a": None if not pt_text else float(pt_text),
                    "fit_included": included,
                }
            )
    rows.sort(key=lambda item: int(item["source_order"]))
    orders = [int(item["source_order"]) for item in rows]
    if orders != list(range(1, 55)):
        raise ValueError(
            "local source table must contain source_order 1 through 54 once"
        )
    return rows


def canonical_digest(
    rows: list[dict[str, float | int | bool | None]], *, include_pt: bool = False
) -> str:
    """Hash normalized reported precision, not CSV layout or extra columns."""
    lines = []
    if include_pt and any(item["pt_lattice_a"] is None for item in rows):
        raise ValueError(
            "all 54 Pt lattice parameters are required for the full digest"
        )
    for item in rows:
        fields = [
            str(int(item["source_order"])),
            f"{float(item['pressure_fei_gpa']):.1f}",
            f"{float(item['pressure_holmes_gpa']):.1f}",
            f"{float(item['temperature_k']):.0f}",
            f"{float(item['volume_a3']):.3f}",
        ]
        if include_pt:
            fields.append(f"{float(item['pt_lattice_a']):.3f}")
        fields.append("1" if item["fit_included"] else "0")
        lines.append(",".join(fields))
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def _rmse(residuals: Any) -> float:
    values = np.asarray(residuals, dtype=float)
    return float(np.sqrt(np.mean(values**2)))


def _fit_summary(result: FitResult, names: tuple[str, ...]) -> dict[str, Any]:
    return {
        "parameters": {name: float(result.parameters[name]) for name in names},
        "standard_errors": {
            name: float(result.standard_errors[name]) for name in names
        },
        "pressure_rmse_gpa": _rmse(result.residuals),
        "solver_success": bool(result.success),
    }


def _thermal_arrays(
    rows: list[dict[str, float | int | bool | None]], pressure_name: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = [
        item for item in rows if int(item["source_order"]) in THERMAL_SOURCE_ORDERS
    ]
    return (
        np.asarray([float(item["volume_a3"]) for item in selected]) * CELL_TO_MOLAR,
        np.asarray([float(item["temperature_k"]) for item in selected]),
        np.asarray([float(item[pressure_name]) for item in selected]),
    )


def _static_fit(
    rows: list[dict[str, float | int | bool | None]], pressure_name: str
) -> FitResult:
    selected = [
        item for item in rows if int(item["source_order"]) in STATIC_SOURCE_ORDERS
    ]
    return fit_rt_eos(
        BM2,
        [float(item["volume_a3"]) for item in selected],
        [float(item[pressure_name]) for item in selected],
        initial={"V0": 46.0, "K0": 220.0},
        bounds={"V0": (40.0, 55.0), "K0": (100.0, 400.0)},
    )


def _mgd_fit(
    rows: list[dict[str, float | int | bool | None]],
    pressure_name: str,
    *,
    V0: float,
    K0: float,
    theta0: float | None,
) -> FitResult:
    volume, temperature, pressure = _thermal_arrays(rows, pressure_name)
    if theta0 is None:
        initial = {"theta0": 1200.0, "gamma0": 2.7, "q": 1.0}
        fixed = {"Tr": 700.0, "n": 5.0}
        bounds = {
            "theta0": (100.0, 5000.0),
            "gamma0": (0.1, 8.0),
            "q": (-5.0, 8.0),
        }
    else:
        initial = {"gamma0": 2.7, "q": 1.0}
        fixed = {"Tr": 700.0, "theta0": theta0, "n": 5.0}
        bounds = {"gamma0": (0.1, 8.0), "q": (-5.0, 8.0)}
    return fit_thermal_eos(
        MieGruneisenDebye,
        BM2(V0 * CELL_TO_MOLAR, K0),
        volume,
        temperature,
        pressure,
        initial=initial,
        fixed=fixed,
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds=bounds,
    )


def _published_curve_summary(
    rows: list[dict[str, float | int | bool | None]],
    pressure_name: str,
    model: Any,
) -> dict[str, float]:
    volume, temperature, pressure = _thermal_arrays(rows, pressure_name)
    residuals = np.asarray(model.pressure(volume, temperature)) - pressure
    return {
        "pressure_rmse_gpa": _rmse(residuals),
        "maximum_absolute_residual_gpa": float(np.max(np.abs(residuals))),
    }


def _calibration_summary(
    rows: list[dict[str, float | int | bool | None]],
    pressure_name: str,
    model: Any,
    *,
    public_pt_cell_volume: bool,
    apparent_source_anomalies: set[int],
) -> dict[str, Any]:
    predicted = []
    observed = []
    orders = []
    for item in rows:
        a = item["pt_lattice_a"]
        if a is None:
            raise ValueError("Pt lattice parameter is required for calibration checks")
        cell_volume = float(a) ** 3
        volume = (
            cell_volume if public_pt_cell_volume else cell_volume * CELL_TO_MOLAR / 4
        )
        predicted.append(float(model.pressure(volume, float(item["temperature_k"]))))
        observed.append(float(item[pressure_name]))
        orders.append(int(item["source_order"]))
    residuals = np.asarray(predicted) - np.asarray(observed)
    retained = np.asarray(
        [order not in apparent_source_anomalies for order in orders], dtype=bool
    )
    return {
        "all_rows": {
            "observations": len(rows),
            "pressure_rmse_gpa": _rmse(residuals),
            "maximum_absolute_residual_gpa": float(np.max(np.abs(residuals))),
            "within_0_5_gpa": int(np.sum(np.abs(residuals) <= 0.5)),
        },
        "excluding_apparent_source_anomalies": {
            "excluded_source_orders": sorted(apparent_source_anomalies),
            "observations": int(np.sum(retained)),
            "pressure_rmse_gpa": _rmse(residuals[retained]),
            "maximum_absolute_residual_gpa": float(np.max(np.abs(residuals[retained]))),
            "within_0_5_gpa": int(np.sum(np.abs(residuals[retained]) <= 0.5)),
        },
    }


def refit_source_table(
    rows: list[dict[str, float | int | bool | None]],
) -> dict[str, Any]:
    """Run the published two-stage protocol and both pressure-scale branches."""
    fei_static = _static_fit(rows, "pressure_fei_gpa")
    holmes_static = _static_fit(rows, "pressure_holmes_gpa")
    model1 = _mgd_fit(rows, "pressure_fei_gpa", V0=46.5, K0=207.0, theta0=None)
    model2 = _mgd_fit(rows, "pressure_fei_gpa", V0=46.5, K0=207.0, theta0=1000.0)
    model3 = _mgd_fit(rows, "pressure_holmes_gpa", V0=45.8, K0=238.0, theta0=1300.0)
    volume, temperature, pressure = _thermal_arrays(rows, "pressure_fei_gpa")
    model4 = fit_thermal_eos(
        LogVolumeThermalPressure,
        BM2(46.5 * CELL_TO_MOLAR, 207.0),
        volume,
        temperature,
        pressure,
        initial={"alpha_KT_ref": 0.012, "dK_dT_V": -0.01},
        fixed={"Tr": 700.0},
        bounds={"alpha_KT_ref": (-0.1, 0.1), "dK_dT_V": (-0.2, 0.2)},
    )
    models = {
        "model_1_fei_mgd": {
            "reference_isotherm_refit": _fit_summary(fei_static, ("V0", "K0")),
            "thermal_refit": _fit_summary(model1, ("theta0", "gamma0", "q")),
            "published_curve": _published_curve_summary(
                rows,
                "pressure_fei_gpa",
                MieGruneisenDebye(
                    BM2(46.5 * CELL_TO_MOLAR, 207.0), 700.0, 1300.0, 2.7, 1.2, 5.0
                ),
            ),
        },
        "model_2_fei_fixed_theta_mgd": {
            "thermal_refit": _fit_summary(model2, ("gamma0", "q")),
            "published_curve": _published_curve_summary(
                rows,
                "pressure_fei_gpa",
                MieGruneisenDebye(
                    BM2(46.5 * CELL_TO_MOLAR, 207.0), 700.0, 1000.0, 2.7, 1.6, 5.0
                ),
            ),
        },
        "model_3_holmes_fixed_theta_mgd": {
            "reference_isotherm_refit": _fit_summary(holmes_static, ("V0", "K0")),
            "thermal_refit": _fit_summary(model3, ("gamma0", "q")),
            "published_curve": _published_curve_summary(
                rows,
                "pressure_holmes_gpa",
                MieGruneisenDebye(
                    BM2(45.8 * CELL_TO_MOLAR, 238.0), 700.0, 1300.0, 2.8, 2.1, 5.0
                ),
            ),
        },
        "model_4_fei_log_volume": {
            "thermal_refit": {
                **_fit_summary(model4, ("alpha_KT_ref", "dK_dT_V")),
                "derived_alpha_ref_per_k": float(
                    model4.parameters["alpha_KT_ref"] / 207.0
                ),
                "derived_alpha_ref_standard_error_per_k": float(
                    model4.standard_errors["alpha_KT_ref"] / 207.0
                ),
            },
            "published_curve": _published_curve_summary(
                rows,
                "pressure_fei_gpa",
                LogVolumeThermalPressure(
                    BM2(46.5 * CELL_TO_MOLAR, 207.0),
                    700.0,
                    5.7e-5 * 207.0,
                    -0.010,
                ),
            ),
        },
    }
    result: dict[str, Any] = {
        "canonical_transcription_sha256": canonical_digest(rows),
        "models": models,
    }
    if all(item["pt_lattice_a"] is not None for item in rows):
        fei_scale = MieGruneisenDebye(
            BM3(60.38 * CELL_TO_MOLAR / 4.0, 273.0, 4.8),
            Tr=300.0,
            theta0=230.0,
            gamma0=2.69,
            q=0.5,
            n=1.0,
            debye_temperature_law="integrated_gruneisen",
        )
        platinum = get_material_document("platinum")
        holmes_scale = Material.from_eosmat(
            platinum, record_identifiers=[HOLMES_RECORD_ID]
        ).eos_records[0]
        result["canonical_transcription_with_pt_sha256"] = canonical_digest(
            rows, include_pt=True
        )
        result["pressure_calibration_checks"] = {
            "fei_2004": _calibration_summary(
                rows,
                "pressure_fei_gpa",
                fei_scale,
                public_pt_cell_volume=False,
                apparent_source_anomalies={4, 54},
            ),
            "holmes_1989": _calibration_summary(
                rows,
                "pressure_holmes_gpa",
                holmes_scale,
                public_pt_cell_volume=True,
                apparent_source_anomalies={4, 48, 54},
            ),
        }
    return result


def verify_against_artifact(result: dict[str, Any], artifact: dict[str, Any]) -> None:
    """Fail on a different transcription or numerical reproduction drift."""
    expected_source = artifact["source"]
    if result["canonical_transcription_sha256"] != expected_source["pvt_sha256"]:
        raise AssertionError(
            "local P-V-T transcription digest does not match the audit"
        )
    if "canonical_transcription_with_pt_sha256" in result and (
        result["canonical_transcription_with_pt_sha256"]
        != expected_source["pvt_and_pt_sha256"]
    ):
        raise AssertionError("local Pt-augmented transcription digest does not match")

    def compare(actual: Any, expected: Any, path: str = "") -> None:
        if isinstance(expected, dict):
            for key, value in expected.items():
                if key in {
                    "published",
                    "method",
                    "notes",
                    "source_anomalies",
                    "qualification",
                }:
                    continue
                if key not in actual:
                    raise AssertionError(f"refit output is missing {path + key}")
                compare(actual[key], value, f"{path}{key}.")
        elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
            if not np.isclose(float(actual), float(expected), rtol=2e-7, atol=2e-10):
                raise AssertionError(
                    f"refit drift at {path[:-1]}: {actual!r} != {expected!r}"
                )

    compare({"models": result["models"]}, {"models": artifact["models"]})
    if "pressure_calibration_checks" in result:
        compare(
            {"pressure_calibration_checks": result["pressure_calibration_checks"]},
            {"pressure_calibration_checks": artifact["pressure_calibration_checks"]},
        )


def _print_shipped_checks(artifact: dict[str, Any]) -> None:
    document = get_material_document("ca_perovskite")
    record = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]
    reference_pressure = record.pressure(46.5, 700.0)
    extrapolated_volume = record.volume(0.0, 300.0)
    extrapolated_modulus = isothermal_bulk_modulus(record, extrapolated_volume, 300.0)
    print(f"record={RECORD_ID}")
    print(f"audit artifact={ARTIFACT_PATH.relative_to(ROOT)}")
    print(f"source rows={artifact['selection']['source_rows']}")
    print(f"excluded source orders={artifact['selection']['excluded_source_orders']}")
    print(f"reference P(46.5 A^3, 700 K)={reference_pressure:.12g} GPa")
    print(
        "300 K zero-pressure extrapolation: "
        f"V={extrapolated_volume:.8f} A^3, K={extrapolated_modulus:.8f} GPa "
        "(paper: 45.8 A^3, 225 GPa)"
    )
    print("source-fit audit status=parity")
    print(artifact["finding"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-table",
        type=Path,
        help="lawfully obtained local CSV transcription of source Table 1",
    )
    parser.add_argument(
        "--json", action="store_true", help="print the fresh local refit result as JSON"
    )
    args = parser.parse_args()
    artifact = load_artifact()
    if args.source_table is None:
        _print_shipped_checks(artifact)
        return 0
    result = refit_source_table(load_source_table(args.source_table))
    verify_against_artifact(result, artifact)
    print("local transcription and all rerun diagnostics match the audit artifact")
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
