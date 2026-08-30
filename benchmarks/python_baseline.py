"""Record reproducible operation-level performance for the Python backend."""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import scipy

import peritheos
from peritheos.eos.rt import BM2, BM3, Holzapfel
from peritheos.eos.thermal import MieGruneisenDebye, Sokolova2016
from peritheos.fitting import fit_rt_eos
from peritheos.uncertainty import EOSUncertainty


def _consume(value: Any) -> float:
    """Reduce a benchmark result so every eager result is inspected."""
    if hasattr(value, "parameters"):
        return float(sum(value.parameters.values()))
    if hasattr(value, "standard_error"):
        return float(np.sum(np.asarray(value.standard_error, dtype=float)))
    return float(np.sum(np.asarray(value, dtype=float)))


def _measure(
    operation: Callable[[], Any], *, repeats: int, warmups: int = 1
) -> dict[str, Any]:
    checksum = 0.0
    for _ in range(warmups):
        checksum = _consume(operation())
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        result = operation()
        samples.append(time.perf_counter() - started)
        checksum = _consume(result)
    return {
        "median_seconds": statistics.median(samples),
        "minimum_seconds": min(samples),
        "maximum_seconds": max(samples),
        "samples_seconds": samples,
        "checksum": checksum,
    }


def _workloads(*, quick: bool) -> dict[str, tuple[Callable[[], Any], str]]:
    bm3 = BM3(50.0, 130.0, 4.3)
    scalar_iterations = 5_000 if quick else 50_000
    array_size = 20_000 if quick else 1_000_000
    inverse_size = 100 if quick else 10_000
    debye_size = 100 if quick else 10_000
    sokolova_size = 10 if quick else 100

    array_volumes = np.linspace(30.0, 60.0, array_size)
    inverse_pressures = np.linspace(-1.0, 300.0, inverse_size)

    debye = MieGruneisenDebye(
        BM3(1.0, 160.0, 4.0),
        Tr=300.0,
        theta0=800.0,
        gamma0=1.5,
        q=1.0,
        n=2.0,
    )
    debye_volumes = np.linspace(0.7, 1.0, debye_size)
    debye_temperatures = np.linspace(300.0, 3_000.0, debye_size)

    sokolova = Sokolova2016(
        Holzapfel(0.3414, 441.5, 3.9, 1.0, 6),
        Tr=298.15,
        QE1o=684.0,
        mE1=0.564,
        QE2o=1561.0,
        mE2=2.436,
        delta=-0.506,
        t=1.085,
        a_0=5.2,
        m=1.3,
        g=2.7,
        e_0=0.8,
        beta=0.35,
        QBo=480.0,
        d=2.4,
        mb=0.75,
        QB1o=1120.0,
        d1=1.6,
        mb1=0.4,
    )
    sokolova_volumes = np.linspace(0.24, 0.3414, sokolova_size)
    sokolova_temperatures = np.linspace(300.0, 4_000.0, sokolova_size)

    fit_points = 20 if quick else 100
    true_model = BM3(10.0, 120.0, 4.3)
    true_volumes = np.linspace(8.0, 10.5, fit_points)
    measured_volumes = true_volumes + 0.008 * np.sin(np.arange(fit_points))
    observed_pressures = true_model.pressure(true_volumes)

    uncertainty_points = 20 if quick else 1_000
    uncertainty_volumes = np.linspace(8.0, 10.0, uncertainty_points)
    uncertainty = EOSUncertainty(
        BM2(10.0, 120.0), parameter_errors={"V0": 0.02, "K0": 2.0}
    )
    monte_carlo_samples = 100 if quick else 2_000

    def scalar_pressure() -> float:
        result = 0.0
        for _ in range(scalar_iterations):
            result += bm3.pressure(40.0)
        return result

    def ordinary_fit():
        return fit_rt_eos(
            BM3,
            true_volumes,
            observed_pressures,
            initial={"K0": 110.0, "K0_prime": 4.0},
            fixed={"V0": 10.0},
            pressure_sigma=0.01,
            absolute_sigma=True,
        )

    def errors_in_variables_fit():
        return fit_rt_eos(
            BM3,
            measured_volumes,
            observed_pressures,
            initial={"K0": 110.0, "K0_prime": 4.0},
            fixed={"V0": 10.0},
            pressure_sigma=0.002,
            volume_sigma=0.01,
            absolute_sigma=True,
        )

    return {
        "bm3_pressure_scalar": (
            scalar_pressure,
            f"{scalar_iterations} scalar public method calls",
        ),
        "bm3_pressure_array": (
            lambda: bm3.pressure(array_volumes),
            f"{array_size} contiguous volumes",
        ),
        "bm3_volume_array": (
            lambda: bm3.volume(inverse_pressures),
            f"{inverse_size} independently bracketed roots",
        ),
        "debye_pressure_array": (
            lambda: debye.pressure(debye_volumes, debye_temperatures),
            f"{debye_size} volume-temperature states",
        ),
        "sokolova_pressure_array": (
            lambda: sokolova.pressure(sokolova_volumes, sokolova_temperatures),
            f"{sokolova_size} volume-temperature states",
        ),
        "fit_rt_ordinary": (ordinary_fit, f"{fit_points} observations"),
        "fit_rt_errors_in_variables": (
            errors_in_variables_fit,
            f"{fit_points} observations with latent volumes",
        ),
        "uncertainty_linear_pressure": (
            lambda: uncertainty.pressure(
                uncertainty_volumes, volume_sigma=0.005, full_covariance=False
            ),
            f"{uncertainty_points} states and two uncertain parameters",
        ),
        "uncertainty_monte_carlo_pressure": (
            lambda: uncertainty.pressure(
                np.array([8.5, 9.0, 9.5]),
                volume_sigma=0.005,
                method="monte_carlo",
                sample_count=monte_carlo_samples,
                random_state=42,
            ),
            f"{monte_carlo_samples} samples at three states",
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--repeats", type=int, default=None)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    repeats = options.repeats or (2 if options.quick else 5)
    if repeats < 1:
        parser.error("--repeats must be positive")

    results = {}
    for name, (operation, description) in _workloads(quick=options.quick).items():
        result = _measure(operation, repeats=repeats)
        result["description"] = description
        results[name] = result
        print(f"{name}: {result['median_seconds']:.6f} s", file=sys.stderr)

    payload = {
        "schema_version": 1,
        "backend": "python",
        "quick": options.quick,
        "repeats": repeats,
        "environment": {
            "peritheos": peritheos.__version__,
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
        },
        "benchmarks": results,
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if options.output is not None:
        options.output.parent.mkdir(parents=True, exist_ok=True)
        options.output.write_text(serialized, encoding="utf-8")
    print(serialized, end="")


if __name__ == "__main__":
    main()
