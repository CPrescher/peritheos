#!/usr/bin/env python3
"""Rounding/convention audit and grouped validation of a research-only refit.

No source coefficients, observations or executable EOSMAT records are replaced.
Run from the repository root with PYTHONPATH=. Optional --plot needs matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares, minimize

from peritheos.eos.experimental.maltby_2024 import (
    N_A_ANGSTROM,
    Maltby2024Parameters,
    Maltby2024Published,
    Maltby2024Trial,
    R,
    fcc_shells,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/argon-fcc-dewaele-2021-supplement.csv"
OUTPUT = ROOT / "docs/data/argon-maltby-2024-refit.json"
PUBLISHED = Maltby2024Parameters()
# Half of the last explicitly printed decimal place; these are rounding
# intervals, NOT parameter standard errors. Trailing-zero ambiguity is separate.
HALF_LAST_DIGIT = dict(
    epsilon_k=0.05,
    alpha_r=0.005,
    rmin_a=0.0005,
    lambda_k_a9=50,
    z1_k=0.5,
    z2=0.005,
    z3=0.0005,
    z4_cm3=0.05,
    vref_cm3=0.005,
    theta_d0_k=0.5,
    a_d_k=0.05,
    b_d=0.000005,
    c_d=0.00000005,
    gamma_d0=0.0005,
    q_d=0.00005,
    b1=0.00000005,
    b2=0.0000000005,
    b3=0.0000000005,
    c1_k=0.00005,
    c2=0.0005,
    c3_k=0.000005,
    c4=0.005,
    a0=0.00005,
    a1=0.000005,
    a2=0.000005,
    theta0_k=0.005,
    theta1_k=0.5,
    theta2_k=0.005,
    gamma0=0.0005,
    gamma1=0.0000000005,
    gamma2=0.0005,
)
# Visually checked transcription of official SI Table SI.1, PDF page 2.
SI_COUNTS = [
    12,
    6,
    24,
    12,
    24,
    8,
    48,
    6,
    36,
    24,
    24,
    24,
    72,
    48,
    48,
    12,
    48,
    30,
    72,
    24,
    48,
    24,
    48,
    8,
    84,
    24,
    96,
    48,
    24,
    0,
    96,
    6,
    96,
    48,
    48,
    36,
    120,
    24,
    48,
    24,
    48,
    48,
    120,
    24,
    120,
    0,
    96,
    24,
    108,
    30,
    48,
    72,
    72,
    32,
    144,
    0,
    96,
    72,
    72,
    48,
    120,
    0,
    144,
    12,
]


def rounding_audit():
    base = Maltby2024Published(64)
    individual = []
    slopes = []
    for name, half in HALF_LAST_DIGIT.items():
        endpoints = [
            Maltby2024Trial(
                64, replace(PUBLISHED, **{name: getattr(PUBLISHED, name) + s * half})
            ).pressure(2.397, 70)
            for s in (-1, 1)
        ]
        slopes.append((endpoints[1] - endpoints[0]) / 2)
        individual.append(
            dict(
                parameter=name,
                half_last_digit=half,
                endpoint_pressure_mpa=np.array(endpoints) * 1000,
            )
        )
    intervals = []
    for trailing_zero in (False, True):
        widths = dict(HALF_LAST_DIGIT)
        if trailing_zero:
            widths.update(z1_k=5.0, theta1_k=5.0)
        names = list(widths)

        def trial(x):
            return Maltby2024Trial(
                64,
                replace(
                    PUBLISHED,
                    **{
                        n: getattr(PUBLISHED, n) + widths[n] * x[i]
                        for i, n in enumerate(names)
                    },
                ),
            )

        extremal = []
        for direction in (1, -1):
            # Multistart bounded optimization; numerical extrema, not interval proofs.
            solutions = []
            for start in (np.zeros(len(names)), -direction * np.sign(slopes)):
                result = minimize(
                    lambda x: direction * trial(x).pressure(2.397, 70) * 1000,
                    start,
                    method="L-BFGS-B",
                    bounds=[(-1, 1)] * len(names),
                    options={"ftol": 1e-13, "gtol": 1e-8, "maxiter": 1000},
                )
                if not result.success:
                    raise RuntimeError(result.message)
                solutions.append(result)
            result = min(solutions, key=lambda s: s.fun)
            model = trial(result.x)
            extremal.append(
                dict(
                    pressure_mpa=model.pressure(2.397, 70) * 1000,
                    volume_at_1mpa_cm3_mol=model.volume(0.001, 70) * 10,
                    pressure_over_rounded_volume_endpoints_mpa=[
                        model.pressure(v / 10, 70) * 1000 for v in (23.965, 23.975)
                    ],
                    parameters=asdict(model.parameters),
                )
            )
        intervals.append(
            dict(trailing_zero_integers_half_width_5=trailing_zero, extrema=extremal)
        )
    # A one-point match is a diagnostic only, never a fitting observation.
    rmin = brentq(
        lambda x: Maltby2024Trial(64, replace(PUBLISHED, rmin_a=x)).pressure(2.397, 70)
        - 0.001,
        3.78,
        3.83,
    )
    return dict(
        baseline_pressure_mpa=base.pressure(2.397, 70) * 1000,
        one_at_a_time=individual,
        numerical_joint_extrema=intervals,
        one_point_rmin_diagnostic=dict(
            rmin_a=rmin,
            change_in_printed_last_digits=(rmin - PUBLISHED.rmin_a) / 0.001,
            used_in_refit=False,
        ),
    )


def convention_audit():
    model = Maltby2024Published(64)
    actual = dict(zip(*fcc_shells(64)))
    differences = [
        dict(squared_distance=i, printed_count=n, geometric_count=int(actual.get(i, 0)))
        for i, n in enumerate(SI_COUNTS, 1)
        if actual.get(i, 0) != n
    ]
    literal = Maltby2024Trial(64, PUBLISHED)
    literal._shells = np.array(
        [i for i, n in enumerate(SI_COUNTS, 1) if n], dtype=float
    )
    literal._populations = np.array([n for n in SI_COUNTS if n], dtype=float)
    v, t = 23.97, 70
    rnn = (np.sqrt(2) * v / N_A_ANGSTROM) ** (1 / 3)
    rc = 8 * rnn
    eps, steepness, rmin = PUBLISHED.epsilon_k, PUBLISHED.alpha_r, PUBLISHED.rmin_a
    decay = steepness / rmin
    radii = np.sqrt(model._shells) * rnn
    rep = eps * 6 / (steepness - 6) * np.exp(steepness - decay * radii)
    att = eps * steepness / (steepness - 6) * (rmin / radii) ** 6
    b = eps * steepness / (steepness - 6) * rmin**6
    rep_rc = eps * 6 / (steepness - 6) * np.exp(steepness - decay * rc)
    prefactor = 2 * np.pi * N_A_ANGSTROM / v
    tail = prefactor * (
        rep_rc * (rc**2 / decay + 2 * rc / decay**2 + 2 / decay**3) - b / (3 * rc**3)
    )
    pair = np.dot(model._populations, rep - att) / 2 + tail
    correction = PUBLISHED.lambda_k_a9 * N_A_ANGSTROM / (eps * model.sigma_angstrom**6)
    baseline = model.pressure(v / 10, t)
    omitted_density = baseline + R * pair * correction / v**2 * 0.001
    fixed_cutoff = (
        baseline
        - R
        * (1 - correction / v)
        * prefactor
        * (rep_rc - b / rc**6)
        * rc**3
        / (3 * v)
        * 0.001
    )
    scan = [
        dict(
            cutoff_squared=int(c),
            pressure_mpa=Maltby2024Published(int(c)).pressure(v / 10, t) * 1000,
        )
        for c in fcc_shells(256)[0]
    ]
    return dict(
        shell_table_discrepancies=differences,
        geometric_pressure_mpa=baseline * 1000,
        literal_si_table_pressure_mpa=literal.pressure(v / 10, t) * 1000,
        literal_si_table_volume_cm3_mol=literal.volume(0.001, t) * 10,
        omitted_explicit_density_derivative_pressure_mpa=omitted_density * 1000,
        fixed_instantaneous_cutoff_derivative_pressure_mpa=fixed_cutoff * 1000,
        all_occupied_cutoffs_1_to_256=scan,
        interpretation="Only the geometric fcc sum and full derivative of Eqs.12-16 are used for fitting; alternatives are diagnostics, not adopted corrections.",
    )


def load_observations():
    selected, excluded = [], []
    with DATA.open() as stream:
        for line, row in enumerate(csv.DictReader(stream), 2):
            row_id = f"dewaele2021_line_{line}"
            if not row["pressure_gpa"]:
                excluded.append(
                    dict(
                        row_id=row_id,
                        sample=row["sample"],
                        csv_line=line,
                        reason="missing pressure",
                    )
                )
                continue
            pressure, temperature = (
                float(row["pressure_gpa"]),
                float(row["temperature_k"]),
            )
            if not 0 < pressure <= 16 or not 0 <= temperature <= 300:
                excluded.append(
                    dict(
                        row_id=row_id,
                        sample=row["sample"],
                        csv_line=line,
                        reason="outside diagnostic P,T limits",
                    )
                )
                continue
            selected.append(
                dict(
                    row_id=row_id,
                    sample=row["sample"],
                    csv_line=line,
                    run=int(row["run"]),
                    pressure_gpa=pressure,
                    temperature_k=temperature,
                    volume_j_bar_mol=float(row["argon_a_angstrom"]) ** 3
                    * N_A_ANGSTROM
                    / 40,
                )
            )
    return selected, excluded


def prediction(model, rows):
    return np.array(
        [model.pressure(r["volume_j_bar_mol"], r["temperature_k"]) for r in rows]
    )


def metrics(predicted, rows):
    observed = np.array([r["pressure_gpa"] for r in rows])
    residual = predicted - observed
    runs = np.array([r["run"] for r in rows])
    return dict(
        count=len(rows),
        rms_pressure_gpa=np.sqrt(np.mean(residual**2)),
        bias_gpa=np.mean(residual),
        max_absolute_gpa=np.max(abs(residual)),
        mean_absolute_relative_percent=np.mean(abs(residual / observed)) * 100,
        equal_run_rms_gpa=np.sqrt(
            np.mean([np.mean(residual[runs == run] ** 2) for run in np.unique(runs)])
        ),
    )


def fit_model(rows, weighting="equal_run", epsilon_fraction=0.05):
    """Only epsilon and rmin vary; bounds are research choices, not errors."""
    runs = np.array([r["run"] for r in rows])
    observed = np.array([r["pressure_gpa"] for r in rows])
    weights = (
        np.ones(len(rows))
        if weighting == "equal_row"
        else np.array([1 / np.sqrt(sum(runs == run)) for run in runs])
    )

    def model(x):
        return Maltby2024Trial(
            64,
            replace(
                PUBLISHED,
                epsilon_k=PUBLISHED.epsilon_k * (1 + epsilon_fraction * x[0]),
                rmin_a=PUBLISHED.rmin_a * (1 + 0.01 * x[1]),
            ),
        )

    results = [
        least_squares(
            lambda x: (prediction(model(x), rows) - observed) * weights,
            start,
            bounds=(-np.ones(2), np.ones(2)),
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
        )
        for start in ([0, 0], [-0.5, 0.5], [0.5, -0.5])
    ]
    if not all(r.success for r in results):
        raise RuntimeError("A bounded fit did not converge")
    result = min(results, key=lambda r: r.cost)
    candidate = model(result.x)
    singular = np.linalg.svd(result.jac, compute_uv=False)
    return candidate, dict(
        fitted_parameters={
            n: getattr(candidate.parameters, n) for n in ("epsilon_k", "rmin_a")
        },
        normalized_parameters=result.x,
        active_bounds=result.active_mask,
        scaled_jacobian_singular_values=singular,
        scaled_jacobian_condition=singular[0] / singular[-1],
        optimizer="scipy.optimize.least_squares, bounded TRF, three starts",
        multistart_costs=[r.cost for r in results],
        weighting=weighting,
        epsilon_relative_bound=epsilon_fraction,
        rmin_relative_bound=0.01,
    )


def cross_validate(rows, weighting="equal_run", epsilon_fraction=0.05):
    baseline = Maltby2024Published(64)
    out_of_fold = np.empty(len(rows))
    folds = []
    for held_run in sorted({r["run"] for r in rows}):
        train = [r for r in rows if r["run"] != held_run]
        indices = [i for i, r in enumerate(rows) if r["run"] == held_run]
        test = [rows[i] for i in indices]
        model, info = fit_model(train, weighting, epsilon_fraction)
        out_of_fold[indices] = prediction(model, test)
        folds.append(
            dict(
                held_out_run=held_run,
                train_rows=[r["row_id"] for r in train],
                test_rows=[r["row_id"] for r in test],
                fit=info,
                train=metrics(prediction(model, train), train),
                test=metrics(out_of_fold[indices], test),
                baseline_test=metrics(prediction(baseline, test), test),
            )
        )
    model, fit = fit_model(rows, weighting, epsilon_fraction)
    return model, dict(
        final_fit=fit,
        baseline=metrics(prediction(baseline, rows), rows),
        training=metrics(prediction(model, rows), rows),
        grouped_validation=metrics(out_of_fold, rows),
        folds=folds,
        predictions=[
            dict(
                **r,
                baseline_pressure_gpa=prediction(baseline, [r])[0],
                fitted_pressure_gpa=prediction(model, [r])[0],
                out_of_fold_pressure_gpa=out_of_fold[i],
            )
            for i, r in enumerate(rows)
        ],
    )


def thermodynamic_properties(model, volume, temperature):
    h = 0.01
    cv = (
        -temperature
        * (
            model.molar_helmholtz_energy(volume, temperature + h)
            + model.molar_helmholtz_energy(volume, temperature - h)
            - 2 * model.molar_helmholtz_energy(volume, temperature)
        )
        / h**2
    )
    bulk = model.bulk_modulus(volume, temperature)
    alpha = (
        model.pressure(volume, temperature + h)
        - model.pressure(volume, temperature - h)
    ) / (2 * h * bulk)
    cp = cv + alpha**2 * bulk * volume * temperature * 1e4
    return dict(
        pressure_gpa=model.pressure(volume, temperature),
        bulk_modulus_gpa=bulk,
        cv_j_mol_k=cv,
        cp_j_mol_k=cp,
        alpha_mk_inverse=alpha * 1000,
        gamma=alpha * bulk * volume * 1e4 / cv,
        chi_t_gpa_inverse=1 / bulk,
        chi_s_gpa_inverse=cv / cp / bulk,
    )


def property_audit(model, observations):
    fixed = thermodynamic_properties(model, 2.397, 70)
    volume = model.volume(0.001, 70)
    at_pressure = thermodynamic_properties(model, volume, 70)
    rows = [
        thermodynamic_properties(model, r["volume_j_bar_mol"], r["temperature_k"])
        for r in observations
    ]
    return dict(
        table8_fixed_volume=fixed,
        table8_at_1mpa=dict(volume_cm3_mol=volume * 10, **at_pressure),
        observed_states_stability=dict(
            count=len(rows),
            positive_bulk_modulus=all(r["bulk_modulus_gpa"] > 0 for r in rows),
            positive_cv=all(r["cv_j_mol_k"] > 0 for r in rows),
            cp_at_least_cv=all(r["cp_j_mol_k"] >= r["cv_j_mol_k"] for r in rows),
            minimum_bulk_modulus_gpa=min(r["bulk_modulus_gpa"] for r in rows),
            minimum_cv_j_mol_k=min(r["cv_j_mol_k"] for r in rows),
        ),
        caveat="Numerical stability and Table8 theory checks, not independent caloric measurements or phase-equilibrium validation.",
    )


def external_ono(model):
    path = ROOT / "peritheos/data/datasets/argon-fcc-ono-2020-table1.csv"
    rows = []
    with path.open() as stream:
        for line, r in enumerate(csv.DictReader(stream), 2):
            if float(r["pressure_gpa"]) > 16:
                continue
            v = float(r["volume_a3"]) * N_A_ANGSTROM / 40
            p = model.pressure(v, 300)
            sp, sv = float(r["pressure_sigma_gpa"]), float(r["volume_sigma_a3"])
            # One-sigma diagonal propagation is diagnostic: no covariance given.
            combined = np.hypot(
                sp, model.bulk_modulus(v, 300) * sv / float(r["volume_a3"])
            )
            rows.append(
                dict(
                    csv_line=line,
                    pressure_gpa=float(r["pressure_gpa"]),
                    calculated_pressure_gpa=p,
                    residual_gpa=p - float(r["pressure_gpa"]),
                    combined_diagonal_sigma_gpa=combined,
                    normalized_residual=(p - float(r["pressure_gpa"])) / combined,
                )
            )
    return dict(
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        rows=rows,
        rms_pressure_gpa=np.sqrt(np.mean([r["residual_gpa"] ** 2 for r in rows])),
        caveat="Three external 300K points only; different pressure calibration and stress history. Not used to select parameters or fit options.",
    )


def plot_report(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), layout="constrained")
    primary = report["experiments"]["source_informed_equal_run"]
    folds = primary["folds"]
    x = np.arange(len(folds))
    axes[0].bar(
        x - 0.18,
        [f["baseline_test"]["rms_pressure_gpa"] for f in folds],
        0.36,
        label="Printed coefficients",
        color="#596e87",
    )
    axes[0].bar(
        x + 0.18,
        [f["test"]["rms_pressure_gpa"] for f in folds],
        0.36,
        label="Refit; run held out",
        color="#d18b25",
    )
    axes[0].set(
        xticks=x,
        xticklabels=[f"Run {f['held_out_run']}" for f in folds],
        ylabel="Pressure RMS error (GPa)",
        title="Run-held-out validation",
    )
    axes[0].legend(fontsize=8)
    rows = primary["predictions"]
    for run in range(1, 6):
        subset = [r for r in rows if r["run"] == run]
        axes[1].scatter(
            [r["pressure_gpa"] for r in subset],
            [r["out_of_fold_pressure_gpa"] - r["pressure_gpa"] for r in subset],
            s=18,
            label=f"Run {run}",
        )
    axes[1].axhline(0, color="gray", lw=0.8)
    axes[1].set(
        xlabel="Observed pressure (GPa)",
        ylabel="Held-out prediction − observation (GPa)",
        title="Residuals at original observations",
    )
    axes[1].legend(fontsize=8)
    checks = report["property_checks"]
    volumes = [
        checks[n]["table8_at_1mpa"]["volume_cm3_mol"]
        for n in ("published", "source_informed_refit")
    ]
    axes[2].bar(
        [0, 1], (np.array(volumes) / 23.97 - 1) * 100, color=["#596e87", "#d18b25"]
    )
    axes[2].axhline(0, color="gray", lw=0.8)
    axes[2].set(
        xticks=[0, 1],
        xticklabels=["Printed\ncoefficients", "95-row\nrefit"],
        ylabel="Volume deviation from Table 8 (%)",
        title="70 K, 1 MPa: unfitted theory check",
    )
    fig.suptitle(
        "Maltby-based diagnostic refit: pressure improvement has a low-pressure tradeoff",
        fontsize=13,
    )
    fig.savefig(ROOT / "docs/data/argon-maltby-2024-refit.png", dpi=160)
    plt.close(fig)


def json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    observations, excluded = load_observations()
    # Dewaele's source-informed exclusion of stressed pure-Ar high-P runs.
    # This is a new diagnostic subset, not Maltby's undocumented 38-row selection.
    screened = [
        r for r in observations if r["run"] not in (1, 4) or r["pressure_gpa"] <= 5
    ]
    print("Auditing rounding and lattice conventions...", flush=True)
    report = dict(
        status="diagnostic_refit_not_promoted",
        identifier="argon_fcc_maltby_2024_diagnostic_refit",
        source_data_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),
        observation_count=len(observations),
        excluded_rows=excluded,
        source_informed_excluded_rows=[
            r["row_id"] for r in observations if r not in screened
        ],
        rounding=rounding_audit(),
        conventions=convention_audit(),
        published_parameters=asdict(PUBLISHED),
        fixed_parameters=[
            n for n in asdict(PUBLISHED) if n not in ("epsilon_k", "rmin_a")
        ],
        cutoff_squared=64,
        experiments={},
        property_checks={},
        external_ono={},
        protocol=dict(
            row_identity="Source sample labels are not unique; row_id combines dataset and original CSV line. Repeated labels and observations remain separate rows.",
            objective="Sum of squared pressure residuals; equal_run divides each row's squared residual by its training-run row count. Equal_row is a sensitivity check.",
            primary="source_informed_equal_run",
            grouping="Leave one complete run out; preprocessing and weights use training runs only; no row shuffling or tuning on held-out outcomes.",
            candidate_selection="Two static coefficients selected on physical roles; all thermal and zero-point terms fixed. The four experiments are sensitivity analyses, not a model-selection competition.",
            primary_selection="All finite 0<P<=16 GPa, 0<=T<=300 K observations, then exclude run1/run4 P>5 GPa following source hydrostatic selection rationale. Repeated rows retained.",
            uncertainty="No Dewaele coordinate uncertainties/covariance provided; weights and parameter bounds are declared numerical choices, not statistical confidence intervals.",
            independence="Maltby's published coefficients used some Dewaele data. Grouped validation measures held-out performance of this new refit; it is not a completely independent validation of the original paper.",
        ),
    )
    primary_model = None
    for name, rows, weighting, bound in [
        ("source_informed_equal_run", screened, "equal_run", 0.05),
        ("source_informed_equal_row", screened, "equal_row", 0.05),
        ("all_130_equal_run", observations, "equal_run", 0.05),
        ("source_informed_wider_epsilon", screened, "equal_run", 0.10),
    ]:
        print(f"Fitting {name}: {len(rows)} rows, five held-out runs...", flush=True)
        model, results = cross_validate(rows, weighting, bound)
        report["experiments"][name] = results
        if name == "source_informed_equal_run":
            primary_model = model
    for name, model in [
        ("published", Maltby2024Published(64)),
        ("source_informed_refit", primary_model),
    ]:
        report["property_checks"][name] = property_audit(model, observations)
        report["external_ono"][name] = external_ono(model)
    report["candidate_parameters"] = asdict(primary_model.parameters)
    report["limitations"] = [
        "Pressure improvement alone is insufficient for promotion; inspect active bounds, grouped errors and deterioration of unfitted thermodynamic checks.",
        "Full primary multi-property data, full-precision source coefficients and author implementation were not obtained; no author outreach sent.",
        "Geometric shell14 is absent, unlike official SI.1. The literal table worsens the sample discrepancy. No convention alternative is silently adopted.",
        "No phase-equilibrium calculation, complete multiproperty refit, calibrated covariance or native fitted EOS is claimed.",
    ]
    OUTPUT.write_text(
        json.dumps(report, indent=2, default=json_default, allow_nan=False) + "\n"
    )
    if args.plot:
        plot_report(report)
    print(f"Wrote {OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
