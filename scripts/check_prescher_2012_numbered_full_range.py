"""Repository-only BM3 comparison for all numbered Prescher p1-p19 points.

Run ``python -m scripts.check_prescher_2012_numbered_full_range``.
The final Table S1 coordinates and published ruby widths are preserved.
"""

from __future__ import annotations

import csv
import hashlib
import json

import numpy as np
import scipy
from scipy import odr
from scipy.optimize import least_squares

from scripts.check_prescher_2012_coordinate_errors import coordinate_fit, odr_check
from scripts.recover_prescher_2012_fe3c import (
    ROOT,
    SOURCE,
    TABLE,
    bm3,
    fit,
    polynomial_pressure_fit,
    volume_at_pressure,
)
from scripts.recover_prescher_2012_fe3c_splits import matches_printed

FIXED_V0 = 154.5678863704812
AUDIT = (
    ROOT / "peritheos/data/research/cementite/prescher-2012-numbered-pressure-audit.csv"
)
REPORT = ROOT / "docs/data/prescher-2012-numbered-full-range-fit.json"
NOTE = ROOT / "docs/literature-reproductions/prescher-2012-numbered-full-range-fit.md"


def fixed_coordinate_fit(p, v, sp, sv, v0=FIXED_V0):
    """Fit K0, K0-prime and latent volumes, conditioning on exact fixed V0."""
    p, v, sp, sv = (np.asarray(x, dtype=float) for x in (p, v, sp, sv))

    def residual(b):
        return np.r_[(bm3(b[2:], v0, *b[:2]) - p) / sp, (b[2:] - v) / sv]

    opt = least_squares(
        residual,
        np.r_[[160.0, 6.0], v],
        jac="cs",
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
        max_nfev=5000,
    )
    if not opt.success or np.linalg.matrix_rank(opt.jac) != len(opt.x):
        raise RuntimeError("Fit did not converge with full-rank covariance")
    cov = np.linalg.inv(opt.jac.T @ opt.jac)[:2, :2]
    errors = np.sqrt(np.diag(cov))
    dof = len(p) - 2
    sse = float(opt.fun @ opt.fun)
    return {
        "parameters": dict(zip(("V0", "K0", "K0_prime"), [v0, *opt.x[:2]])),
        "free_parameters": ["K0", "K0_prime"],
        "fixed_parameters": {"V0": v0},
        "absolute_errors": errors.tolist(),
        "absolute_covariance": cov.tolist(),
        "residual_scaled_errors": (errors * np.sqrt(sse / dof)).tolist(),
        "residual_scaled_covariance": (cov * sse / dof).tolist(),
        "observations": len(p),
        "degrees_of_freedom": dof,
        "weighted_sum_squares": sse,
        "reduced_chi_square": sse / dof,
        "latent_table_volumes_a3": opt.x[2:].tolist(),
        "success": bool(opt.success),
    }


def fixed_odr_check(p, v, sp, sv, v0=FIXED_V0):
    result = odr.ODR(
        odr.RealData(p, v, sx=sp, sy=sv),
        odr.Model(lambda b, pressure: volume_at_pressure(pressure, *b)),
        beta0=[v0, 160, 6],
        ifixb=[0, 1, 1],
        maxit=300,
        sstol=1e-10,
        partol=1e-10,
    ).run()
    if result.info not in (1, 2, 3):
        raise RuntimeError(result.stopreason)
    return {
        "parameters_V0_K0_K0prime": result.beta.tolist(),
        "absolute_covariance_K0_K0prime": result.cov_beta[1:, 1:].tolist(),
        "weighted_sum_squares": float(result.sum_square),
        "stop_reason": result.stopreason,
    }


def replay():
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    with AUDIT.open() as stream:
        audit = list(csv.DictReader(stream))
    assert len(table) == 19
    assert [row["run"] for row in audit[:19]] == [f"p{n}" for n in range(1, 20)]
    p, v, sp, sv = (
        table[name]
        for name in (
            "pressure_gpa",
            "volume_a3",
            "pressure_error_gpa",
            "volume_error_a3",
        )
    )
    for source, column in (
        (p, "published_ruby_pressure_gpa"),
        (v, "published_volume_a3"),
        (sp, "published_pressure_width_gpa"),
        (sv, "published_volume_width_a3"),
    ):
        np.testing.assert_array_equal(
            source, [float(row[column]) for row in audit[:19]]
        )

    free = coordinate_fit(p, v, sp, sv)
    free["free_parameters"] = ["V0", "K0", "K0_prime"]
    free["fixed_parameters"] = {}
    free["absolute_errors"] = free.pop("absolute_errors_V0_K0_K0prime")
    free["residual_scaled_errors"] = free.pop("residual_scaled_errors_V0_K0_K0prime")
    check = odr_check(p, v, sp, sv, None)
    np.testing.assert_allclose(
        list(free["parameters"].values()),
        check["parameters_V0_K0_K0prime"],
        rtol=1e-5,
        atol=5e-5,
    )
    np.testing.assert_allclose(
        free["absolute_covariance"], check["absolute_covariance"], rtol=5e-4
    )
    free["independent_ODRPACK_check"] = check

    fixed = fixed_coordinate_fit(p, v, sp, sv)
    check = fixed_odr_check(p, v, sp, sv)
    np.testing.assert_allclose(
        list(fixed["parameters"].values()),
        check["parameters_V0_K0_K0prime"],
        rtol=1e-5,
        atol=5e-5,
    )
    np.testing.assert_allclose(
        fixed["absolute_covariance"],
        check["absolute_covariance_K0_K0prime"],
        rtol=5e-4,
    )
    fixed["independent_ODRPACK_check"] = check
    for case in (free, fixed):
        np.testing.assert_allclose(
            case["weighted_sum_squares"],
            case["independent_ODRPACK_check"]["weighted_sum_squares"],
            rtol=1e-7,
            atol=1e-7,
        )

    equal_free = fit(p, v, initial=[154.5, 160, 6])
    equal_fixed = fit(p, v, fixed_v0=FIXED_V0, initial=[154.5, 160, 6])
    np.testing.assert_allclose(
        list(equal_free["parameters"].values()),
        polynomial_pressure_fit(p, v),
        rtol=2e-5,
        atol=2e-5,
    )
    eta = (FIXED_V0 / v) ** (1 / 3)
    a = 1.5 * (eta**7 - eta**5)
    b = 0.75 * a * (eta**2 - 1)
    k0, k0_times_kp_minus_four = np.linalg.lstsq(
        np.column_stack([a, b]), p, rcond=None
    )[0]
    np.testing.assert_allclose(
        [equal_fixed["parameters"]["K0"], equal_fixed["parameters"]["K0_prime"]],
        [k0, 4 + k0_times_kp_minus_four / k0],
        rtol=1e-6,
    )

    cases = {
        "coordinate_widths_free_V0_primary": free,
        "coordinate_widths_fixed_adopted_V0_control": fixed,
        "equal_pressure_free_V0_control": equal_free,
        "equal_pressure_fixed_adopted_V0_control": equal_fixed,
    }
    target_values, target_errors = SOURCE["full_range"]
    for name, case in cases.items():
        error_key = (
            "absolute_errors"
            if name.startswith("coordinate")
            else "residual_scaled_errors"
        )
        values = [case["parameters"][key] for key in ("K0", "K0_prime")]
        errors = case[error_key][-2:]
        case["comparison_error_convention"] = error_key
        case["matches_published_printed_K0_Kp_errorK0_errorKp"] = matches_printed(
            values + errors, target_values + target_errors
        )
    return {
        "format": "peritheos.prescher-2012-numbered-full-range-fit",
        "date": "2026-10-30",
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "input_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (TABLE, AUDIT)
        },
        "protocol": {
            "model": "Third-order Birch-Murnaghan; pressure in GPa, conventional cell volume in A^3; nominal 293 K",
            "selection": "All numbered p1-p19, 1.8-47.4 GPa, final Table S1 volumes and published ruby pressures. No ambient, graph-only, neon-recalibrated or p20-p22 points.",
            "p16": "133.25 A^3 from final Table S1, corroborated by the separate red.p4p export; not 133.58 A^3 from the alternate SHELX cell.",
            "primary_objective": "Sum[((P_BM3(V_latent)-P_observed)/width_P)^2 + ((V_latent-V_observed)/width_V)^2]; V0,K0,K0_prime and latent volumes free",
            "primary_errors": "Conditional formal 1-sigma errors, treating printed coordinate widths as independent absolute standard deviations; coefficient block of full inverse normal matrix including latent-volume nuisance parameters, no residual scaling",
            "secondary_errors": "Same covariance multiplied by weighted SSE/(19-number_of_free_EOS_coefficients); retained separately",
            "confidence_qualification": "Published parameter confidence, calibration covariance and coordinate covariance are unrecovered. Ruby widths are before/after differences, not established standard deviations. Neither formal convention recovers author uncertainty semantics.",
            "fixed_V0_control": "154.5678863704812 A^3 adopted from the drawn calculated full-range curve, held exact; not the measured starting/ambient volume or a published printed coefficient",
            "equal_pressure_controls": "Same 19 observations, equal pressure weights and exact input volumes; errors use residual covariance SSE/(19-number_of_free_EOS_coefficients). These isolate weighting effects and compare with the existing 22-marker protocol.",
            "verification": "Independent inverse V(P) ODRPACK checks both coordinate fits and covariances; linear BM3 polynomial/design solutions verify equal-pressure coefficients",
        },
        "observations": [
            {
                "run": f"p{n + 1}",
                "pressure_gpa": float(p[n]),
                "volume_a3": float(v[n]),
                "reported_pressure_width_gpa": float(sp[n]),
                "reported_volume_width_a3": float(sv[n]),
            }
            for n in range(19)
        ],
        "published": {
            "K0_gpa": 161,
            "K0_error_gpa": 2,
            "K0_prime": 5.9,
            "K0_prime_error": 0.2,
            "V0": None,
            "parameter_error_confidence": None,
            "doi": "10.1103/PhysRevB.85.140402",
        },
        "cases": cases,
        "author_workflow_reproduced": False,
        "catalog_action": "Independent requested refit only; no published coefficients, source observations, catalog classifications or Studio records changed.",
    }


def write_note(report):
    lines = [
        "# Prescher p1-p19: full-range BM3 refit",
        "",
        "Source: [Prescher et al. (2012), PRB 85, 140402(R)](https://doi.org/10.1103/PhysRevB.85.140402), Figure 1 and final supplementary Table S1.",
        "",
        "All 19 numbered points use their final Table S1 volumes and published ruby pressures (1.8-47.4 GPa, nominal 293 K). The source p16 volume is 133.25 A^3. The ambient candidate is outside this selection.",
        "",
        "The primary fit uses both printed pressure and volume widths, conditionally as independent absolute 1-sigma inputs. V0, K0 and K0-prime are free. Confidence/covariance of the source widths is unknown; ruby widths are before/after reading differences.",
        "",
        "| Fit | V0 (A^3) | K0 (GPa) | K0-prime |",
        "| --- | ---: | ---: | ---: |",
        "| Published full range | not printed | 161 +/- 2 | 5.9 +/- 0.2 |",
    ]
    labels = [
        "Numbered points, coordinate widths, V0 free (primary)",
        "Numbered points, coordinate widths, adopted V0 fixed (control)",
        "Numbered points, equal pressure weights, V0 free (control)",
        "Numbered points, equal pressure weights, adopted V0 fixed (control)",
    ]
    for label, case in zip(labels, report["cases"].values()):
        errors = dict(
            zip(case["free_parameters"], case[case["comparison_error_convention"]])
        )
        cells = []
        for key in ("V0", "K0", "K0_prime"):
            value = case["parameters"][key]
            cells.append(
                f"{value:.5f} +/- {errors[key]:.5f}"
                if key in errors
                else f"{value:.8f} fixed"
            )
        lines.append("| " + label + " | " + " | ".join(cells) + " |")
    primary = next(iter(report["cases"].values()))
    scaled = primary["residual_scaled_errors"]
    lines += [
        "",
        "At the source's printed digits, the primary fit is K0=161 +/- 3 GPa and K0-prime=6.0 +/- 0.2. The fixed-V0 equal-pressure control gives 161 +/- 3 GPa and 5.9 +/- 0.2: both central coefficients match the paper's digits, but error(K0)=2.53219 GPa rounds directly to 3, not the source's 2. Do not double-round intermediate error displays.",
        "",
        f"Primary weighted SSE/(19-3) = {primary['reduced_chi_square']:.5f}. Residual-scaled formal errors are V0 +/- {scaled[0]:.5f} A^3, K0 +/- {scaled[1]:.5f} GPa and K0-prime +/- {scaled[2]:.5f}; these use a different error convention from the primary comparison.",
        "",
        "Fixed-V0 errors condition on exact adopted curve V0 and exclude its uncertainty. The source does not print V0 or recoverably specify whether it was fixed. The old 22-marker parity comparison uses digitized graph coordinates and equal pressure weights; it remains a separate conditional comparison.",
        "",
        "Independent ODRPACK and linear-pressure checks pass. Numerical comparison does not reproduce the author's input set, weighting, covariance or confidence convention. Published coefficients, observations and database/Studio assessments remain unchanged.",
        "",
        "Complete numeric evidence: `docs/data/prescher-2012-numbered-full-range-fit.json`.",
        "",
        "Repository-only replay: `.venv/bin/python -m scripts.check_prescher_2012_numbered_full_range`.",
        "",
    ]
    NOTE.write_text("\n".join(lines))


if __name__ == "__main__":
    report = replay()
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    write_note(report)
    print(REPORT)
    for name, case in report["cases"].items():
        print(name, case["parameters"], case[case["comparison_error_convention"]])
