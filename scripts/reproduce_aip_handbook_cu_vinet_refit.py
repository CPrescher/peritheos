"""Reproduce a diagnostic Vinet refit for the Cu handbook isotherm.

The 25 finite-pressure observations are a manually verified transcription of
the Cu column in the American Institute of Physics Handbook, Table 4d-12.  The
15 kbar entry is printed as ``09900`` and is interpreted as ``0.9900``.
"""

from __future__ import annotations

import csv
import hashlib
from importlib import resources

import numpy as np
from scipy.optimize import minimize

DATASET = "aip-handbook-table4d12-cu.csv"
DATASET_SHA256 = "67adb43b73806d5d1b3e31d4a9d9d645a2fa718a21f83641d4f433a568fb24ce"
V0_ANGSTROM3 = 11.81473546294192
PUBLISHED_K0_GPA = 140.95
PUBLISHED_K0_PRIME = 4.798
STARTS = (
    (PUBLISHED_K0_GPA, PUBLISHED_K0_PRIME),
    (100.0, 4.0),
    (180.0, 6.0),
    (250.0, 3.0),
    (60.0, 8.0),
)


def load_observations() -> tuple[np.ndarray, np.ndarray]:
    """Load the checked pressure and relative-volume cells."""
    resource = resources.files("peritheos").joinpath("data", "datasets", DATASET)
    payload = resource.read_bytes()
    actual_sha256 = hashlib.sha256(payload).hexdigest()
    if actual_sha256 != DATASET_SHA256:
        raise ValueError(f"unexpected checksum for {DATASET}: {actual_sha256}")
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    pressure = np.asarray([float(row["pressure_gpa"]) for row in rows])
    relative_volume = np.asarray([float(row["v_over_v0"]) for row in rows])
    return pressure, relative_volume


def vinet_pressure(
    relative_volume: np.ndarray, k0_gpa: float, k0_prime: float
) -> np.ndarray:
    """Evaluate Vinet pressure directly in terms of V/V0."""
    x = np.cbrt(relative_volume)
    return 3.0 * k0_gpa * (1.0 - x) / x**2 * np.exp(1.5 * (k0_prime - 1.0) * (1.0 - x))


def residual_statistics(
    pressure_gpa: np.ndarray,
    relative_volume: np.ndarray,
    k0_gpa: float,
    k0_prime: float,
) -> dict[str, float]:
    residual = vinet_pressure(relative_volume, k0_gpa, k0_prime) - pressure_gpa
    return {
        "objective8_gpa8": float(np.sum(residual**8)),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "mae_gpa": float(np.mean(np.abs(residual))),
        "max_abs_residual_gpa": float(np.max(np.abs(residual))),
    }


def refit() -> dict[str, object]:
    """Fit K0 and K0' with V0 fixed and the published m=8 objective."""
    pressure, relative_volume = load_observations()

    def objective(parameters: np.ndarray) -> float:
        residual = (
            vinet_pressure(relative_volume, parameters[0], parameters[1]) - pressure
        )
        return float(np.sum(residual**8))

    candidates = [
        minimize(
            objective,
            start,
            method="Nelder-Mead",
            bounds=((50.0, 300.0), (1.0, 10.0)),
            options={"maxiter": 100_000, "xatol": 1.0e-13, "fatol": 1.0e-20},
        )
        for start in STARTS
    ]
    result = min(candidates, key=lambda candidate: float(candidate.fun))
    if not result.success:
        raise RuntimeError(str(result.message))

    k0_gpa, k0_prime = (float(value) for value in result.x)
    return {
        "observations": len(pressure),
        "v0_angstrom3": V0_ANGSTROM3,
        "objective": "sum((P_model - P_table)**8)",
        "refit": {
            "k0_gpa": k0_gpa,
            "k0_prime": k0_prime,
            **residual_statistics(pressure, relative_volume, k0_gpa, k0_prime),
        },
        "published": {
            "k0_gpa": PUBLISHED_K0_GPA,
            "k0_prime": PUBLISHED_K0_PRIME,
            **residual_statistics(
                pressure,
                relative_volume,
                PUBLISHED_K0_GPA,
                PUBLISHED_K0_PRIME,
            ),
        },
    }


def main() -> None:
    result = refit()
    for label in ("refit", "published"):
        values = result[label]
        print(
            f"{label}: K0={values['k0_gpa']:.8f} GPa, "
            f"K0'={values['k0_prime']:.8f}, "
            f"RMSE={values['rmse_gpa']:.8f} GPa, "
            f"L8={values['objective8_gpa8']:.12g} GPa^8"
        )


if __name__ == "__main__":
    main()
