"""Compare selected records through BurnMan and Pytheos public APIs.

This optional maintainer script is deliberately isolated from Peritheos's test
baselines. It neither reads external package data nor treats agreement with
another implementation as scientific validation. Install BurnMan and/or
Pytheos separately, then run this file from the repository root.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any

from scipy.constants import Avogadro

from peritheos.materials import (
    AU_FEI_2007,
    NACL_B2_FEI_2007,
    NE_FEI_2007,
    PT_FEI_2007,
    EOSRecord,
)


@dataclass(frozen=True)
class ComparisonCase:
    label: str
    record: EOSRecord
    formula_units_per_cell: int
    pytheos_factory: Any
    burnman_factory: Any | None
    states: tuple[tuple[float, float], ...]


def _optional_libraries() -> tuple[Any, Any, Any]:
    try:
        import burnman
        import pytheos
        from uncertainties import unumpy
    except ImportError as error:  # pragma: no cover - optional maintainer tool
        raise SystemExit(
            "Install burnman, pytheos, and uncertainties in an isolated "
            "environment before running this optional comparison."
        ) from error
    return burnman, pytheos, unumpy


def main() -> None:
    burnman, pytheos, unumpy = _optional_libraries()
    from burnman import calibrants

    warnings.filterwarnings(
        "ignore",
        message="Using UFloat objects with std_dev==0 may give unexpected results",
    )

    cases = (
        ComparisonCase(
            "Au",
            AU_FEI_2007,
            4,
            pytheos.gold.Fei2007vinet,
            calibrants.Fei_2007.Au,
            ((0.9, 300.0), (0.8, 1000.0)),
        ),
        ComparisonCase(
            "Pt",
            PT_FEI_2007,
            4,
            pytheos.platinum.Fei2007vinet,
            calibrants.Fei_2007.Pt,
            ((0.9, 300.0), (0.8, 1000.0)),
        ),
        ComparisonCase(
            "NaCl-B2",
            NACL_B2_FEI_2007,
            1,
            pytheos.sodium_chloride_b2.Fei2007vinet,
            None,
            ((0.75, 300.0), (0.65, 1000.0)),
        ),
        ComparisonCase(
            "Ne",
            NE_FEI_2007,
            4,
            pytheos.neon.Fei2007vinet,
            None,
            ((0.65, 300.0), (0.5, 1000.0)),
        ),
    )

    print(
        f"# BurnMan {burnman.__version__}; "
        f"Pytheos {getattr(pytheos, '__version__', 'unknown')}"
    )
    print("material,V/V0,T_K,peritheos_GPa,pytheos_GPa,burnman_GPa")
    for case in cases:
        pytheos_scale = case.pytheos_factory()
        burnman_scale = None if case.burnman_factory is None else case.burnman_factory()
        for ratio, temperature in case.states:
            volume = case.record.reference_volume * ratio
            pressure = float(
                case.record.pressure(volume, temperature, check_validity=False)
            )
            pytheos_pressure = float(
                unumpy.nominal_values(pytheos_scale.cal_p(volume, temperature))
            )
            burnman_text = ""
            if burnman_scale is not None:
                molar_volume = volume / case.formula_units_per_cell * Avogadro * 1.0e-30
                burnman_pressure = (
                    float(burnman_scale.pressure(molar_volume, temperature)) / 1.0e9
                )
                burnman_text = f"{burnman_pressure:.9f}"
            print(
                f"{case.label},{ratio:g},{temperature:g},{pressure:.9f},"
                f"{pytheos_pressure:.9f},{burnman_text}"
            )


if __name__ == "__main__":
    main()
