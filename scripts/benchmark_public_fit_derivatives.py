"""Repository-only matched benchmark of public EOS fits.

Run with the development environment and a freshly built native extension:
    python scripts/benchmark_public_fit_derivatives.py --output /tmp/fit-timings.json

Synthetic observations, fit options, model reconstruction and result processing
are identical. Only the private native Jacobian switch is changed. This switch
is a validation aid, not part of the supported public Python fitting API.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import statistics
from pathlib import Path
from time import perf_counter

import numpy as np

from peritheos import _rust
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import QCompromiseDebye
from peritheos.fitting import fit_joint_eos, fit_rt_eos


def cases():
    rng = np.random.default_rng(20261114)
    volume = np.linspace(8.5, 11.3, 80)
    pressure = BM3(11.3, 171.0, 4.7).pressure(volume) + rng.normal(
        0, 0.025, volume.size
    )
    for latent in (False, True):
        yield (
            "BM3 PV" + (" with latent V" if latent else ""),
            fit_rt_eos,
            (BM3, volume, pressure),
            dict(
                initial={"V0": 11.25, "K0": 163.0, "K0_prime": 4.1},
                pressure_sigma=0.04,
                volume_sigma=0.008 if latent else None,
                absolute_sigma=True,
                max_nfev=600,
            ),
        )
    volume = np.tile(np.linspace(8.8, 11.3, 12), 8)
    temperature = np.repeat(np.linspace(250, 2100, 8), 12)
    truth = QCompromiseDebye(BM3(11.3, 171.0, 4.7), 298.0, 820.0, 1.1, 20.0)
    pressure = truth.pressure(volume, temperature) + rng.normal(0, 0.02, volume.size)
    for errors in ("P", "PVT", "correlated PVT"):
        options = dict(
            initial={
                "rt_eos.V0": 11.25,
                "rt_eos.K0": 163.0,
                "rt_eos.K0_prime": 4.3,
                "theta0": 750.0,
                "gamma0": 1.0,
            },
            fixed={"Tr": 298.0, "n": 20.0},
            bounds={"theta0": (100, 1800), "gamma0": (0.2, 2.0)},
            absolute_sigma=True,
            max_nfev=600,
        )
        if errors == "correlated PVT":
            sigma = np.array([0.04, 0.008, 4.0])
            correlation = np.array([[1, 0.25, -0.1], [0.25, 1, 0.2], [-0.1, 0.2, 1]])
            options["observation_covariance"] = correlation * np.outer(sigma, sigma)
        else:
            options["pressure_sigma"] = 0.04
            if errors == "PVT":
                options.update(volume_sigma=0.008, temperature_sigma=4.0)
        yield (
            f"QCompromiseDebye/BM3 joint {errors}",
            fit_joint_eos,
            (QCompromiseDebye, BM3, volume, temperature, pressure),
            options,
        )


def run(repeats):
    rows = []
    original_rt = _rust.fit_rt_eos_native
    original_thermal = _rust.fit_thermal_eos_native
    try:
        for name, function, args, options in cases():
            times = {False: [], True: []}
            results = {}
            # Warm both paths, then alternate order to reduce timing-order bias.
            for iteration in range(repeats + 1):
                for analytic in (False, True) if iteration % 2 == 0 else (True, False):
                    native_calls = []

                    def rt(*args, **kwargs):
                        native_calls.append("isothermal")
                        return original_rt(
                            *args, **kwargs, analytic_derivatives=analytic
                        )

                    def thermal(*args, **kwargs):
                        native_calls.append("thermal")
                        return original_thermal(
                            *args, **kwargs, analytic_derivatives=analytic
                        )

                    _rust.fit_rt_eos_native = rt
                    _rust.fit_thermal_eos_native = thermal
                    start = perf_counter()
                    result = function(*args, **options)
                    elapsed = perf_counter() - start
                    assert len(native_calls) == 1, (
                        "benchmark must use the native EOS fitter"
                    )
                    assert result.success, result.message
                    results[analytic] = result
                    if iteration:
                        times[analytic].append(elapsed)
            analytic, numerical = results[True], results[False]
            np.testing.assert_allclose(
                list(analytic.parameters.values()),
                list(numerical.parameters.values()),
                rtol=2e-5,
            )
            np.testing.assert_allclose(
                analytic.covariance, numerical.covariance, rtol=5e-5, atol=1e-9
            )
            np.testing.assert_allclose(
                analytic.adjusted_volume, numerical.adjusted_volume, atol=2e-7
            )
            if analytic.adjusted_temperature is not None:
                np.testing.assert_allclose(
                    analytic.adjusted_temperature,
                    numerical.adjusted_temperature,
                    atol=2e-4,
                )
            assert abs(analytic.chi_square / numerical.chi_square - 1) < 1e-8
            rows.append(
                dict(
                    case=name,
                    point_count=args[-1].size,
                    analytic_seconds=times[True],
                    finite_difference_seconds=times[False],
                    median_speedup=statistics.median(times[False])
                    / statistics.median(times[True]),
                    results={
                        label: {
                            "parameters": fit.parameters,
                            "standard_errors": fit.standard_errors,
                            "covariance": fit.covariance.tolist(),
                            "chi_square": fit.chi_square,
                            "nfev": fit.nfev,
                            "njev": fit.njev,
                            "adjusted_volume": fit.adjusted_volume.tolist(),
                            "adjusted_temperature": None
                            if fit.adjusted_temperature is None
                            else fit.adjusted_temperature.tolist(),
                        }
                        for label, fit in [
                            ("analytic", analytic),
                            ("finite_difference", numerical),
                        ]
                    },
                )
            )
    finally:
        _rust.fit_rt_eos_native = original_rt
        _rust.fit_thermal_eos_native = original_thermal
    root = Path(__file__).resolve().parents[1]
    files = [
        "crates/peritheos/src/fit/mod.rs",
        "crates/peritheos/src/fit/eos.rs",
        "crates/peritheos/src/fit/derivatives.rs",
        "crates/peritheos-python/src/native_fit.rs",
        "peritheos/fitting.py",
        "scripts/benchmark_public_fit_derivatives.py",
    ]
    return dict(
        scope="Synthetic public P-V/P-V-T fits only; same compiled solver and unchanged objective/uncertainty policy",
        python=platform.python_version(),
        platform=platform.platform(),
        numpy=np.__version__,
        native_sha256=hashlib.sha256(Path(_rust.__file__).read_bytes()).hexdigest(),
        source_sha256={
            name: hashlib.sha256((root / name).read_bytes()).hexdigest()
            for name in files
        },
        repeats=repeats,
        warmups=1,
        cases=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=9)
    arguments = parser.parse_args()
    if arguments.repeats < 1:
        parser.error("--repeats must be positive")
    report = run(arguments.repeats)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    for case in report["cases"]:
        print(f"{case['case']}: {case['median_speedup']:.2f}x")
