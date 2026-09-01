"""Measure Python-callback overhead against end-to-end native EOS fitting."""

from __future__ import annotations

import argparse
import json
import statistics
import time
from collections.abc import Callable
from typing import Any

import numpy as np

from peritheos import _rust
from peritheos.eos.rt import BM3


def _measure(operation: Callable[[], Any], iterations: int, repeats: int):
    operation()
    samples = []
    checksum = 0.0
    for _ in range(repeats):
        started = time.perf_counter()
        for _ in range(iterations):
            result = operation()
            checksum = float(np.sum(result.x))
        samples.append((time.perf_counter() - started) / iterations)
    return {
        "median_seconds": statistics.median(samples),
        "samples_seconds": samples,
        "checksum": checksum,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=7)
    options = parser.parse_args()
    if min(options.points, options.iterations, options.repeats) < 1:
        parser.error("points, iterations, and repeats must be positive")

    volumes = np.linspace(8.0, 10.5, options.points)
    expected = BM3(10.0, 120.0, 4.3)
    pressures = np.asarray(expected.pressure(volumes))
    pressure_sigma = np.full(options.points, 0.01)
    volume_sigma = np.full(options.points, 0.005)
    measured_volumes = volumes + 0.002 * np.sin(np.arange(options.points))
    initial = np.array([110.0, 4.0])
    lower = np.array([-np.inf, -np.inf])
    upper = np.array([np.inf, np.inf])
    prototype = BM3(10.0, 110.0, 4.0)._native

    def residual(parameters):
        model = BM3(10.0, float(parameters[0]), float(parameters[1]))
        return (np.asarray(model.pressure(volumes)) - pressures) / pressure_sigma

    def callback_ordinary():
        return _rust.fit_least_squares(residual, initial, lower, upper, loss="linear")

    def native_ordinary():
        return _rust.fit_rt_eos_native(
            prototype,
            ("K0", "K0_prime"),
            initial,
            lower,
            upper,
            pressures,
            volumes,
            pressure_sigma,
            loss="linear",
        )

    latent_initial = np.concatenate([initial, measured_volumes])
    latent_lower = np.concatenate(
        [lower, np.full(options.points, np.finfo(float).tiny)]
    )
    latent_upper = np.concatenate([upper, np.full(options.points, np.inf)])

    def latent_residual(parameters):
        model = BM3(10.0, float(parameters[0]), float(parameters[1]))
        adjusted = parameters[2:]
        return np.concatenate(
            [
                (np.asarray(model.pressure(adjusted)) - pressures) / pressure_sigma,
                (adjusted - measured_volumes) / volume_sigma,
            ]
        )

    def callback_latent():
        return _rust.fit_least_squares(
            latent_residual,
            latent_initial,
            latent_lower,
            latent_upper,
            loss="linear",
            global_parameter_count=2,
            point_count=options.points,
            latent_coordinate_count=1,
        )

    def native_latent():
        return _rust.fit_rt_eos_native(
            prototype,
            ("K0", "K0_prime"),
            initial,
            lower,
            upper,
            pressures,
            measured_volumes,
            pressure_sigma,
            volume_sigma,
            loss="linear",
        )

    operations = {
        "ordinary_callback": callback_ordinary,
        "ordinary_native": native_ordinary,
        "latent_callback": callback_latent,
        "latent_native": native_latent,
    }
    results = {
        name: _measure(operation, options.iterations, options.repeats)
        for name, operation in operations.items()
    }
    results["ordinary_speedup"] = (
        results["ordinary_callback"]["median_seconds"]
        / results["ordinary_native"]["median_seconds"]
    )
    results["latent_speedup"] = (
        results["latent_callback"]["median_seconds"]
        / results["latent_native"]["median_seconds"]
    )
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
