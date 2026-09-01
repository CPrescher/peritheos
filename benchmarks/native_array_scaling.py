"""Measure native array scaling around the parallel execution thresholds."""

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

import peritheos
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye


def _measure(operation: Callable[[], np.ndarray], repeats: int) -> dict[str, Any]:
    checksum = float(np.sum(operation()))
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        result = operation()
        samples.append(time.perf_counter() - started)
        checksum = float(np.sum(result))
    return {
        "median_seconds": statistics.median(samples),
        "minimum_seconds": min(samples),
        "maximum_seconds": max(samples),
        "samples_seconds": samples,
        "checksum": checksum,
    }


def _repeat_count(size: int) -> int:
    if size < 1_000:
        return 100
    if size < 10_000:
        return 50
    if size < 200_000:
        return 20
    return 7


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    bm3 = BM3(10.0, 120.0, 4.3)
    thermal = MieGruneisenDebye(
        BM3(1.0, 160.0, 4.0),
        Tr=300.0,
        theta0=800.0,
        gamma0=1.5,
        q=1.0,
        n=2.0,
    )
    results: dict[str, dict[str, Any]] = {
        "bm3_pressure": {},
        "bm3_volume": {},
        "debye_pressure": {},
    }

    for size in (64, 512, 2_048, 8_192, 32_768, 65_536, 131_072, 1_000_000):
        volumes = np.linspace(7.0, 11.0, size)
        results["bm3_pressure"][str(size)] = _measure(
            lambda volumes=volumes: bm3.pressure(volumes),
            _repeat_count(size),
        )

    for size in (16, 64, 256, 512, 2_048, 10_000):
        pressures = np.linspace(-1.0, 300.0, size)
        results["bm3_volume"][str(size)] = _measure(
            lambda pressures=pressures: bm3.volume(pressures),
            _repeat_count(size),
        )

    for size in (16, 128, 512, 2_048, 10_000):
        volumes = np.linspace(0.7, 1.0, size)
        temperatures = np.linspace(300.0, 3_000.0, size)
        results["debye_pressure"][str(size)] = _measure(
            lambda volumes=volumes, temperatures=temperatures: thermal.pressure(
                volumes, temperatures
            ),
            _repeat_count(size),
        )

    report = {
        "schema_version": 1,
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "numpy": np.__version__,
            "peritheos": peritheos.__version__,
            "executable": sys.executable,
        },
        "benchmarks": results,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if arguments.output is not None:
        arguments.output.write_text(f"{rendered}\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
