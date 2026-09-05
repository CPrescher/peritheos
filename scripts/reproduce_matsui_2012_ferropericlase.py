#!/usr/bin/env python3
"""Reproduce the two high-spin ferropericlase EOS curves of Matsui et al.

Primary source: Equations (1)-(6) and Tables 1-3,
doi:10.2138/am.2012.3937.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from importlib import resources

import numpy as np

from peritheos import get_material_document
from peritheos.materials import Material

CASES = {
    "mg083fe017o": "mg083fe017o_matsui_2012_bm3_mgd_1",
    "mg075fe025o": "mg075fe025o_matsui_2012_bm3_mgd_1",
}


@dataclass(frozen=True)
class Diagnostics:
    rows: int
    included_rows: int
    included_observed_rmse_gpa: float
    included_observed_max_abs_gpa: float
    table_calculated_rmse_gpa: float
    table_calculated_max_abs_gpa: float
    all_observed_rmse_gpa: float
    all_observed_max_abs_gpa: float


def load_data(material_identifier: str) -> np.ndarray:
    """Load every printed source-table row, including high-P exclusions."""
    document = get_material_document(material_identifier)
    dataset = document["datasets"][0]
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return np.array(
        [
            (
                int(row["source_row"]),
                float(row["temperature_k"]),
                float(row["cell_volume_a3"]),
                float(row["pressure_gpa"]),
                float(row["published_calculated_pressure_gpa"]),
                row["fit_included"] == "1",
            )
            for row in rows
        ],
        dtype=[
            ("source_row", int),
            ("temperature", float),
            ("volume", float),
            ("observed", float),
            ("table_calculated", float),
            ("included", bool),
        ],
    )


def diagnostics(material_identifier: str) -> Diagnostics:
    """Compare the published rounded curve with Table 1 or Table 2."""
    document = get_material_document(material_identifier)
    loaded = Material.from_eosmat(
        document,
        record_identifiers=[CASES[material_identifier]],
    ).eos_records[0]
    data = load_data(material_identifier)
    # Scalar evaluation also exercises the native backend when it is installed.
    # It avoids coupling this scientific audit to backend-specific array layout.
    predicted = np.array(
        [
            loaded.pressure(volume, temperature, check_validity=False)
            for volume, temperature in zip(data["volume"], data["temperature"])
        ],
        dtype=float,
    )
    included = data["included"]
    included_residual = predicted[included] - data["observed"][included]
    table_residual = predicted - data["table_calculated"]
    all_residual = predicted - data["observed"]
    return Diagnostics(
        rows=len(data),
        included_rows=int(np.count_nonzero(included)),
        included_observed_rmse_gpa=float(np.sqrt(np.mean(included_residual**2))),
        included_observed_max_abs_gpa=float(np.max(np.abs(included_residual))),
        table_calculated_rmse_gpa=float(np.sqrt(np.mean(table_residual**2))),
        table_calculated_max_abs_gpa=float(np.max(np.abs(table_residual))),
        all_observed_rmse_gpa=float(np.sqrt(np.mean(all_residual**2))),
        all_observed_max_abs_gpa=float(np.max(np.abs(all_residual))),
    )


def main() -> None:
    for material_identifier in CASES:
        result = diagnostics(material_identifier)
        print(
            f"{material_identifier}: {result.rows} rows ({result.included_rows} fitted)"
        )
        print(
            "  fitted rows vs observed: "
            f"RMSE={result.included_observed_rmse_gpa:.6f} GPa, "
            f"max_abs={result.included_observed_max_abs_gpa:.6f} GPa"
        )
        print(
            "  all rows vs published Pcalc: "
            f"RMSE={result.table_calculated_rmse_gpa:.6f} GPa, "
            f"max_abs={result.table_calculated_max_abs_gpa:.6f} GPa"
        )
        print(
            "  all rows vs observed (spin crossover diagnostic): "
            f"RMSE={result.all_observed_rmse_gpa:.6f} GPa, "
            f"max_abs={result.all_observed_max_abs_gpa:.6f} GPa"
        )


if __name__ == "__main__":
    main()
