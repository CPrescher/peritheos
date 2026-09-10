#!/usr/bin/env python3
"""Normalize the CC BY 4.0 Zhu et al. Mendeley v3 fit inputs to CSV."""

from __future__ import annotations

import argparse
import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).parents[1]
OUTPUT_ROOT = ROOT / "peritheos" / "data" / "datasets"


@dataclass(frozen=True)
class SourceTable:
    relative_path: str
    output_name: str
    columns: tuple[str, ...]
    rows: int
    sha256: str


TABLES = (
    SourceTable(
        "Au/Au_shock.dat",
        "zhu-2025-au-shock.csv",
        ("volume_a3_conventional_cell", "shock_pressure_gpa", "particle_velocity_m_s"),
        12,
        "fab474e9615ba8708d72d187672d93ad55ae6c15a79bc0c31e8f554beed4c2f5",
    ),
    SourceTable(
        "Au/Au_high_T.dat",
        "zhu-2025-au-zero-pressure-thermal-expansion.csv",
        ("temperature_k", "volume_a3_conventional_cell"),
        10,
        "a803b6f3805dec09e08b2312af466c80250b91239288924d061d029e35e56d5d",
    ),
    SourceTable(
        "Pt/Pt_shock_700.dat",
        "zhu-2025-pt-shock.csv",
        ("volume_a3_conventional_cell", "shock_pressure_gpa", "particle_velocity_m_s"),
        68,
        "5591f904a2c0e7200f5e9a6f0a3151607247b187d8e5473ee900990ea6526238",
    ),
    SourceTable(
        "Pt/Pt_high_T.dat",
        "zhu-2025-pt-zero-pressure-thermal-expansion.csv",
        ("temperature_k", "volume_a3_conventional_cell"),
        17,
        "986ce36b663c065f639fce79ffcd6f3f896d76f2740bad10a0e1d5eeb1647aa8",
    ),
    SourceTable(
        "MgO/MgO_pvt_40_new.dat",
        "zhu-2025-mgo-pvt.csv",
        ("pressure_gpa", "temperature_k", "volume_a3_conventional_cell"),
        213,
        "7fce0f2c195d0d3a2bb7f59e22064c283a6da64c259823466bb8616207ae34a0",
    ),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(source_root: Path, output_root: Path = OUTPUT_ROOT) -> dict[str, str]:
    """Write stable CSV normalizations and return their SHA-256 digests."""
    output_root.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    for table in TABLES:
        source = source_root / table.relative_path
        digest = _sha256(source)
        if digest != table.sha256:
            raise ValueError(
                f"Unexpected SHA-256 for {table.relative_path}: {digest}"
            )
        rows = []
        for line in source.read_text(encoding="utf-8-sig").splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("%"):
                continue
            values = stripped.split()
            if len(values) != len(table.columns):
                raise ValueError(
                    f"{table.relative_path} row has {len(values)} values; "
                    f"expected {len(table.columns)}"
                )
            rows.append(values)
        if len(rows) != table.rows:
            raise ValueError(
                f"{table.relative_path} has {len(rows)} active rows; "
                f"expected {table.rows}"
            )

        output = output_root / table.output_name
        with output.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(table.columns)
            writer.writerows(rows)
        outputs[table.output_name] = _sha256(output)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source_root",
        type=Path,
        help="Extracted Mendeley v3 root containing the Au, Pt, and MgO folders",
    )
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    args = parser.parse_args()
    for name, digest in normalize(args.source_root, args.output_root).items():
        print(f"{digest}  {name}")


if __name__ == "__main__":
    main()
