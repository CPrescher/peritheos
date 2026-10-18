"""Diagnostic replay of all ten Shen--Smith Table II EOSs; no catalog edits.

Run with the repository virtualenv Python. Pressure and Vinet formulae are
transcribed independently; compare their evaluations with the current engine.
Fit parameters are diagnostic, not recommended pressure-scale replacements.
"""

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import Avogadro
from scipy.optimize import least_squares, minimize_scalar

from peritheos import get_eos_record

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / "peritheos/data/datasets/shen-smith-2026-table-s1-volumes.csv"
# Phase, record, experiments, cutoff, literal Table II V0/K0/K0-prime.
SPECS = [
    ("Pt", "platinum_shen_2026_vinet_2", {"DAC-2"}, None, 60.364, 261.2, 5.75),
    ("Au", "gold_shen_2026_vinet_3", {"DAC-1"}, None, 67.792, 167.5, 5.85),
    ("Ta", "tantalum_shen_2026_vinet_2", {"DAC-1"}, None, 36.14, 195.2, 3.62),
    ("W", "tungsten_shen_2026_vinet_3", {"DAC-1", "DAC-2"}, None, 31.704, 307.0, 4.02),
    (
        "Mo",
        "molybdenum_shen_2026_vinet_1",
        {"DAC-1", "DAC-2"},
        None,
        31.14,
        260.6,
        4.06,
    ),
    ("MgO", "mgo_shen_2026_vinet_3", {"DAC-2"}, None, 74.636, 161.9, 4.08),
    ("NaCl-B1", "nacl_b1_shen_2026_vinet_1", {"DAC-2"}, "below30", 179.41, 23.5, 5.23),
    ("NaCl-B2", "nacl_b2_shen_2026_vinet_2", {"DAC-2"}, "above35", 41.35, 24.92, 5.63),
    ("bcc-Fe", "fe_shen_2026_vinet_1", {"DAC-1", "DAC-2"}, None, 23.506, 162.1, 5.4),
    (
        "hcp-Fe",
        "iron_shen_2026_vinet_2",
        {"DAC-1", "DAC-2"},
        None,
        22.352,
        168.55,
        5.53,
    ),
]


def copper(v):
    v0 = 4 * 63.546 / (8.939 * Avogadro) * 1e24
    y = 1 - (np.asarray(v) / v0) ** (1 / 3)
    return 3 * 133.6 * y / (1 - y) ** 2 * np.exp(6.29 * y + 2.06 * y * y + 1.65 * y**3)


def vinet(v, v0, k0, kp):
    x = (np.asarray(v) / v0) ** (1 / 3)
    return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def stats(r):
    return dict(
        mean_gpa=float(r.mean()),
        rms_gpa=float(np.sqrt(np.mean(r * r))),
        min_gpa=float(r.min()),
        max_gpa=float(r.max()),
        positive=int(np.sum(r > 0)),
        negative=int(np.sum(r < 0)),
    )


rows = list(csv.DictReader(SOURCE.open()))
cu = {(r["experiment"], r["run_number"]): r for r in rows if r["phase"] == "Cu"}
assert len(cu) == 422
results = []
fig, axes = plt.subplots(5, 2, figsize=(12, 14), constrained_layout=True)
for ax, (phase, rid, experiments, cutoff, v0, k0, kp) in zip(axes.flat, SPECS):
    paired = []
    for source_row, r in enumerate(rows, 2):
        if r["phase"] != phase or r["experiment"] not in experiments:
            continue
        calibrant = cu[(r["experiment"], r["run_number"])]
        p = float(copper(float(calibrant["unit_cell_volume_a3"])))
        if cutoff == "below30" and p >= 30:
            continue
        if cutoff == "above35" and p <= 35:
            continue
        paired.append(
            dict(source_csv_line=source_row, sample=r, copper=calibrant, pressure_gpa=p)
        )
    p = np.array([r["pressure_gpa"] for r in paired])
    v = np.array([float(r["sample"]["unit_cell_volume_a3"]) for r in paired])
    pub = vinet(v, v0, k0, kp)
    model = get_eos_record(rid).eos
    assert model.parameter_values() == dict(V0=v0, K0=k0, K0_prime=kp)
    engine_error = float(np.max(abs(pub - model.pressure(v))))
    engine_cu_error = float(
        np.max(
            abs(
                p
                - get_eos_record("copper_fratanduono_2020_vinet3_298k").eos.pressure(
                    np.array(
                        [float(r["copper"]["unit_cell_volume_a3"]) for r in paired]
                    )
                )
            )
        )
    )
    assert max(engine_error, engine_cu_error) < 1e-9

    def residual(pars):
        return vinet(v, v0, *pars) - p

    fit = least_squares(
        residual,
        [k0, kp],
        bounds=([1, 0], [1000, 15]),
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-11,
        max_nfev=5000,
    )
    assert fit.success

    # Independent one-dimensional profile: eliminate K0 analytically.
    def profile(kprime):
        f = vinet(v, v0, 1.0, kprime)
        kval = float(f @ p / (f @ f))
        return float(np.sum((kval * f - p) ** 2))

    prof = minimize_scalar(
        profile, bounds=(0, 15), method="bounded", options={"xatol": 1e-12}
    )
    f = vinet(v, v0, 1.0, prof.x)
    profile_pars = [float(f @ p / (f @ f)), float(prof.x)]
    assert np.allclose(fit.x, profile_pars, rtol=2e-6, atol=2e-6)
    refit = vinet(v, v0, *fit.x)
    # Upper decimal-rounding envelope for printed K0/K0-prime; V0 stays fixed.
    dk = 0.005 if phase in {"NaCl-B2", "hcp-Fe"} else 0.05
    dkp = 0.05 if phase == "bcc-Fe" else 0.005
    entry = dict(
        phase=phase,
        record_identifier=rid,
        observations=len(p),
        selected_experiments=sorted(experiments),
        pressure_cutoff=cutoff,
        pressure_range_gpa=[float(p.min()), float(p.max())],
        fixed_V0_a3=v0,
        published=dict(K0_gpa=k0, K0_prime=kp),
        refit=dict(K0_gpa=float(fit.x[0]), K0_prime=float(fit.x[1])),
        published_residuals=stats(pub - p),
        refit_residuals=stats(refit - p),
        independent_profile_parameters=profile_pars,
        engine_pressure_max_difference_gpa=engine_error,
        engine_copper_max_difference_gpa=engine_cu_error,
        maximum_K_Kprime_decimal_rounding_pressure_change_gpa=float(
            np.max(abs(vinet(v, v0, k0 + dk, kp + dkp) - pub))
        ),
        pairs=paired,
    )
    entry["published_residuals_by_experiment"] = {
        exp: stats(
            (pub - p)[np.array([r["sample"]["experiment"] == exp for r in paired])]
        )
        for exp in sorted(experiments)
    }
    results.append(entry)
    ax.scatter(
        p, pub - p, s=12, alpha=0.65, color="#c14f36", label="Published Table II"
    )
    ax.scatter(p, refit - p, s=9, alpha=0.65, color="#227eaa", label="Diagnostic refit")
    ax.axhline(0, color="#777777", linewidth=0.7)
    ax.set_title(f"{phase} · n={len(p)}", loc="left", fontweight="bold")
    ax.text(
        0.03,
        0.96,
        f"RMS: {stats(pub - p)['rms_gpa']:.3f} → {stats(refit - p)['rms_gpa']:.3f} GPa",
        transform=ax.transAxes,
        va="top",
        fontsize=10,
    )
    ax.set_ylim(-3.2, 3.2)
    ax.set_xlim(0, 145)
    ax.set_xlabel("Paired Cu pressure (GPa)")
    ax.set_ylabel("P(model) − P(Cu) (GPa)")
    ax.grid(alpha=0.15)
axes[0, 0].legend(loc="lower left", fontsize=9)
fig.suptitle(
    "Shen & Smith (2026): direct replay of all ten Table II EOSs\n"
    "Fixed published V₀; unweighted pressure residuals; adopted Cu scale held fixed",
    fontsize=14,
)
fig.savefig(OUT / "residuals.png", dpi=160)
fig.savefig(OUT / "residuals.svg")
payload = dict(
    date="2026-10-17",
    source_doi="10.1103/fxgq-96sg",
    source_article_sha256="fb101fe78237c91077586d8999103a10c7ee75813a13ada441ee81a880c83f72",
    source_csv=str(SOURCE),
    source_csv_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    source_csv_rows=len(rows),
    author_workflow_reproduced=False,
    qualification="Diagnostic refits only; no validated replacement standards. "
    "Source weighting, covariance and exact NaCl coexistence mask unavailable. "
    "Uses bundled Table S1 transcription, not a newly downloaded workbook. "
    "Both Pt first/last readings retained; shared Cu scale and observations are correlated.",
    source_reporting_findings={
        "Au_Figure_7_page_8": {
            "fixed_cell_V0_a3": 67.792,
            "PV_K0_gpa": 168.76,
            "PV_K0_prime": 5.82,
            "KV_K0_gpa": 167.67,
            "KV_K0_prime": 5.93,
        },
        "Au_Table_II_page_9": {
            "fixed_cell_V0_a3": 67.792,
            "PV_K0_gpa": 167.5,
            "PV_K0_prime": 5.85,
        },
        "MgO_Table_II_page_10": "Literal printed coefficients and V0 visually verified in final PDF.",
        "NaCl_Table_I_page_4_vs_Table_II_page_11": "Approximately 30-35 GPa coexistence excluded vs B2 above 33 GPa; exact author row mask unavailable.",
    },
    interpretation_limits=[
        "A nearby unweighted coefficient optimum does not establish the author workflow or calibration accuracy.",
        "W retains large structured residuals with different DAC trends; hcp Fe and other phases also retain residual structure.",
        "No formal experimental significance is assigned using these unweighted residuals: both Cu and phase volumes have errors and unknown covariance.",
        "MgO and NaCl-B2 direct unweighted refits differ materially from the printed curves; cause unresolved.",
    ],
    records=results,
)
(OUT / "results.json").write_text(json.dumps(payload, indent=2) + "\n")
for e in results:
    print(
        f"{e['phase']:8s} n={e['observations']:3d} "
        f"RMS {e['published_residuals']['rms_gpa']:.6f} -> {e['refit_residuals']['rms_gpa']:.6f} "
        f"K0={e['refit']['K0_gpa']:.6f} Kp={e['refit']['K0_prime']:.6f}"
    )
