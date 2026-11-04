"""Plot-ready replay of non-P–V audit objectives, with native residual units.

Only stored coefficients and source readers/equations are used. No regression
is performed here. A joint fit gets separate panels for its different outputs.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


@dataclass
class DiagnosticPanel:
    title: str
    x: np.ndarray
    observed: np.ndarray
    source: np.ndarray
    refit: np.ndarray
    xlabel: str
    ylabel: str
    unit: str
    swap_axes: bool = False
    sigma: np.ndarray | None = None
    x_sigma: np.ndarray | None = None
    # (label, independent coordinate, source prediction, refit prediction)
    curves: list = field(default_factory=list)
    observed_label: str = "Audit input rows"
    residual_label: str = "Model − input"
    statistic_label: str = "RMS error"

    def rms(self, side):
        return float(np.sqrt(np.mean((getattr(self, side) - self.observed) ** 2)))


@dataclass
class DiagnosticComparison:
    panels: list[DiagnosticPanel]
    note: str
    source_title: str = "Original / source fit"
    independent_refit: bool = False
    partial_replay: bool = False


@lru_cache(maxsize=1)
def ledger():
    return {
        r["record_identifier"]: r
        for r in json.loads(
            (ROOT / "docs/data/primary-eos-refits.json").read_text(encoding="utf-8")
        )["records"]
    }


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def values(rows, key):
    return np.array([float(row[key]) for row in rows])


def coefficients(record, side):
    return {
        row["parameter"]: row[side]
        for row in record["reproduction"]["assessment"]["parameters"]
    }


def panel(title, x, observed, source, refit, xlabel, ylabel, unit, **kwargs):
    arrays = [np.asarray(a, dtype=float) for a in (x, observed, source, refit)]
    if not all(a.shape == arrays[0].shape for a in arrays) or not np.all(
        np.isfinite(arrays)
    ):
        raise ValueError(f"{title}: nonfinite or mismatched replay arrays")
    return DiagnosticPanel(title, *arrays, xlabel, ylabel, unit, **kwargs)


def diagnostic_comparison(document, record):
    rid = record["identifier"]
    fit = record["reproduction"]["assessment"].get("fit", {})
    kind = fit.get("fit_kind")
    if kind == "simultaneous_cocompression_vinet":
        return dorfman_comparison(record)
    if kind == "acoustic_finite_strain_then_isothermal_conversion":
        return li_comparison(document, record)
    if rid == "kcl_b2_ma_2024_bm3_mgd":
        return ma_comparison(record)
    if kind == "joint_bm3_mgd_finite_strain_pressure_bulk_shear":
        return fu_comparison(record)
    if kind == "coupled_bm3_three_level_spin_crossover":
        return spin_comparison(record)
    if kind == "mie_gruneisen_hugoniot_bounded_reconstruction":
        return akins_comparison(record)
    if kind == "joint_static_pvt_shock_hugoniot_mgd":
        return mosenfelder_comparison(record)
    if kind == "sound_velocity_quasi_debye_bounded_partial":
        return luo_comparison(document, record)
    if kind == "approximate_coupled_vinet_mgd_global_reconstruction":
        return tange_comparison(record)
    if kind == "iterative_robust_bisquare_gruneisen_refit" and rid in {
        "gold_zhu_2025_pvt",
        "platinum_zhu_2025_pvt",
    }:
        return zhu_expansion_comparison(record)
    if record["eos"]["type"] == "LinearUsUpHugoniot":
        from scripts.validate_primary_eos_refits import _load_rows

        dataset = next(
            d
            for d in document["datasets"]
            if d["identifier"] == fit["dataset_identifiers"][0]
        )
        rows = _load_rows(dataset)
        columns = fit.get("columns", {})
        if rid == "iron_brown_2000_linear_hugoniot":
            xkey, ykey = "particle_velocity_km_s", "shock_velocity_km_s"
        else:
            xkey, ykey = columns["particle_velocity"], columns["shock_velocity"]
        if rid == "mgo_b1_duffy_ahrens_1995_hugoniot_5":
            rows = [r for r in rows if str(r["used_in_published_hugoniot_fit"]) == "1"]
        x, y = values(rows, xkey), values(rows, ykey)
        if len(x) != fit["observations"]:
            raise ValueError("Shock-velocity row selection differs from the audit")
        grid = np.linspace(x.min(), x.max(), 160)
        coefs = [coefficients(record, side) for side in ("published", "refit")]
        predictions = [c["c0"] + c["s"] * x for c in coefs]
        curves = [("Linear Us–up fit", grid, *[c["c0"] + c["s"] * grid for c in coefs])]
        return DiagnosticComparison(
            [
                panel(
                    "Shock-velocity fit",
                    x,
                    y,
                    *predictions,
                    "Particle velocity (km/s)",
                    "Shock velocity (km/s)",
                    "km/s",
                    curves=curves,
                    sigma=values(rows, columns["shock_velocity_sigma"])
                    if columns.get("shock_velocity_sigma")
                    else None,
                    x_sigma=values(rows, columns["particle_velocity_sigma"])
                    if columns.get("particle_velocity_sigma")
                    else None,
                )
            ],
            "The fitted quantity is shock velocity, not isothermal pressure. Residual = Us(model) − Us(input). All selected final shock states are retained; HEL precursor columns are excluded. RMS is at observed particle velocities and can differ from an errors-in-variables latent-coordinate statistic. No P–V branch extrapolation is used.",
        )
    return None


def dorfman_comparison(record):
    from scripts.reproduce_dorfman_2012_cocompression import (
        DATASET_FILENAME,
        PUBLISHED_FREE_PARAMETERS,
        REFERENCE_PARAMETERS,
        _pairs,
        load_dataset,
        vinet_pressure,
    )

    mode = "fixed" if "_k0_fixed_" in record["identifier"] else "free"
    report = json.loads(
        (ROOT / "docs/data/dorfman-2012-cocompression-refit.json").read_text(
            encoding="utf-8"
        )
    )
    fit = report["fits"][mode]
    source = REFERENCE_PARAMETERS if mode == "fixed" else PUBLISHED_FREE_PARAMETERS
    fitted = fit["parameters"]
    rows = load_dataset(ROOT / "peritheos/data/datasets" / DATASET_FILENAME)
    pairs = _pairs(rows)
    differences, normalized = [], []
    for model in (source, fitted):
        left = np.array(
            [
                vinet_pressure(vl, model[left_name])
                for left_name, vl, right_name, vr in pairs
            ]
        )
        right = np.array(
            [
                vinet_pressure(vr, model[right_name])
                for left_name, vl, right_name, vr in pairs
            ]
        )
        differences.append(left - right)
        normalized.append((left - right) / np.sqrt((left + right) / 2))
    for actual, expected in zip(
        normalized, (fit["published_objective"], fit["objective"])
    ):
        if not np.isclose(np.sum(actual**2), expected, rtol=1e-10):
            raise ValueError("Dorfman joint objective does not match stored audit")
    index = np.arange(1, len(pairs) + 1)
    # Pair equality is the measured co-compression constraint; there is no
    # independently observed pressure column shared by all phases.
    panels = [
        panel(
            "All phase-pair pressure differences",
            index,
            np.zeros(len(pairs)),
            *differences,
            "Source pair index",
            "P(left) − P(right) (GPa)",
            "GPa",
            observed_label="Equal-pressure constraint",
        ),
        panel(
            "Equation (3) residuals",
            index,
            np.zeros(len(pairs)),
            *normalized,
            "Source pair index",
            "(P(left) − P(right))/√Pmean (√GPa)",
            "√GPa",
            observed_label="Equal-pressure constraint",
        ),
    ]
    material = (
        "Au"
        if record["identifier"].startswith("gold")
        else "Mo"
        if record["identifier"].startswith("molybdenum")
        else "Pt"
    )
    volume = np.array(
        [
            r["values"][material]["volume_a3_conventional_cell"]
            for r in rows
            if material in r["values"]
        ]
    )
    grid = np.linspace(volume.min(), volume.max(), 160)
    # A curve-only panel is deliberately labelled as calculated, not as data.
    p0, p1 = [
        np.array([vinet_pressure(v, model[material]) for v in grid])
        for model in (source, fitted)
    ]
    panels.insert(
        0,
        panel(
            f"{material} curve difference (calculated, not data residual)",
            grid,
            p0,
            p0,
            p1,
            "Volume (Å³/conventional cell)",
            "Pressure (GPa)",
            "GPa",
            swap_axes=True,
            curves=[("Stored Vinet curves", grid, p0, p1)],
            observed_label="Calculated source baseline",
            residual_label="Curve − source curve",
            statistic_label="RMS curve difference",
        ),
    )
    return DiagnosticComparison(
        panels,
        f"Joint-system replay of all {len(rows)} source rows and {len(pairs)} phase pairs. The first panel compares calculated material curves, with the source curve as its baseline; its lower panel is a curve difference, not an experimental fit residual. The following two panels show the actual co-compression equality constraint and Equation (3) residuals; summing the squared normalized residuals reproduces the stored joint objective. No independent P–V observations are invented. The MgO anchor, source row order and fixed/free K0 mode are unchanged. Missing run AN012, source covariance and original unrounded inputs remain unresolved.",
    )


def li_comparison(document, record):
    from scripts.validate_primary_eos_refits import _load_rows

    row = ledger()[record["identifier"]]
    dataset = next(
        d
        for d in document["datasets"]
        if d["identifier"] == row["dataset_identifiers"][0]
    )
    rows = [
        r
        for r in _load_rows(dataset)
        if r["experimental_path"] in {"ambient", "decompression_after_annealing"}
    ]
    density = values(rows, "density_g_cm3")
    inputs = record["scientific_validation"]["isothermal_parameter_derivation"][
        "inputs"
    ]
    correction = 1 + inputs["alpha_per_k"] * inputs["gamma0"] * inputs["temperature_k"]

    def prediction(side, rho):
        c = {r["parameter"]: r[side] for r in row["acoustic_parameters"]}
        ks, kp, g, gp = [c[n] for n in ("K0S", "K0S_prime", "G0", "G0_prime")]
        strain = (1 - (rho / density[0]) ** (2 / 3)) / 2
        scale = (1 - 2 * strain) ** 2.5
        l1 = ks + 4 * g / 3
        l2 = 5 * l1 - 3 * (ks / correction) * (kp + 4 * gp / 3)
        m2 = 5 * g - 3 * (ks / correction) * gp
        return np.sqrt(scale * (l1 + l2 * strain) / rho), np.sqrt(
            scale * (g + m2 * strain) / rho
        )

    grid = np.linspace(density.min(), density.max(), 160)
    src, dst = prediction("published", density), prediction("refit", density)
    sc, dc = prediction("published", grid), prediction("refit", grid)
    panels = [
        panel(
            f"{wave}-wave acoustic fit",
            density,
            values(rows, name),
            src[i],
            dst[i],
            "Density (g/cm³)",
            "Velocity (km/s)",
            "km/s",
            curves=[("Acoustic finite-strain model", grid, sc[i], dc[i])],
        )
        for i, (wave, name) in enumerate(
            (("P", "p_wave_velocity_km_s"), ("S", "s_wave_velocity_km_s"))
        )
    ]
    return DiagnosticComparison(
        panels,
        "The 11 measured density/VP/VS rows are the regression inputs. Calculated pressures are excluded. Curves replay the audit's explicit finite-strain room-temperature approximation, including K0T in L2/M2. Residuals and displayed RMS cover all 11 velocity rows; the separately recorded audit RMS covers only the ten high-pressure rows, excluding the ambient anchor. The catalog's isothermal K0/K0′ are derived from this acoustic fit; their parameter comparison remains separate.",
    )


def ma_comparison(record):
    from scripts.reproduce_ma_2024_kcl import bulk_modulus, load, pressure

    acoustic, walker = load("acoustic"), load("walker")
    # Public coefficients use Å³; the source reader/equations use cm³/mol.
    vectors = []
    factor = 0.602214076  # KCl record is normalized to one formula unit.
    for side in ("published", "refit"):
        c = coefficients(record, side)
        vectors.append(np.array([c["V0"] * factor, c["K0"], c["K0_prime"]]))
    panels = []
    for title, rows, ykey, sigma, function, swap in (
        ("Acoustic KT–V rows", acoustic, "kt_gpa", "kt_error_gpa", bulk_modulus, False),
        (
            "Recalibrated Walker P–V rows",
            walker,
            "matsui_300k_pressure_gpa",
            "matsui_pressure_error_gpa",
            pressure,
            True,
        ),
    ):
        x, y = values(rows, "molar_volume_cm3_mol"), values(rows, ykey)
        grid = np.linspace(x.min(), x.max(), 160)
        panels.append(
            panel(
                title,
                x,
                y,
                *[function(x, c) for c in vectors],
                "Molar volume (cm³/mol)",
                "Pressure (GPa)" if swap else "Isothermal bulk modulus (GPa)",
                "GPa",
                swap_axes=swap,
                sigma=values(rows, sigma),
                curves=[
                    ("Stored BM3 model", grid, *[function(grid, c) for c in vectors])
                ],
            )
        )
    return DiagnosticComparison(
        panels,
        "Joint replay of 11 measured acoustic-volume KT rows and eight Matsui-recalibrated Walker pressure-volume rows. The two output residuals are shown separately in GPa. Their reported widths remain error bars; the audit's weighted objective and missing original covariance are not inferred from unweighted RMS. Model grids and modeled S1 pressures are excluded.",
    )


def fu_comparison(record):
    from scripts.audit_fu_2023_casio3_refit import (
        PARAMETERS,
        _arrays,
        _predictions,
        load_bundled_observations,
    )

    report = json.loads(
        (ROOT / "docs/data/fu-2023-casio3-refit-audit.json").read_text(encoding="utf-8")
    )
    fitted = next(
        f
        for f in report["sensitivity_fits"]
        if f["objective"] == "unweighted_absolute_gpa"
    )["parameters"]
    source = (
        fitted
        if record.get("record_kind") == "refit"
        else report["published_parameters"]
    )
    vectors = [np.array([c[n] for n in PARAMETERS]) for c in (source, fitted)]
    data = _arrays(*load_bundled_observations())
    panels = []
    for prefix, title in (("sun", "Sun P–V–T rows"), ("g", "Gréaux acoustic rows")):
        x, t = data[prefix + "_v"], data[prefix + "_t"]
        src, dst = [_predictions(c, x, t) for c in vectors]
        grid = np.linspace(x.min(), x.max(), 160)
        temperatures = np.unique(t)
        if len(temperatures) > 5:
            temperatures = temperatures[[0, len(temperatures) // 2, -1]]
        isotherms = [
            (
                temperature,
                [
                    _predictions(c, grid, np.full_like(grid, temperature))
                    for c in vectors
                ],
            )
            for temperature in temperatures
        ]
        outputs = (
            [(0, "Pressure", "p", "ps", True)]
            if prefix == "sun"
            else [
                (0, "Pressure", "p", "ps", True),
                (1, "Bulk modulus", "k", "ks", False),
                (2, "Shear modulus", "mu", "mus", False),
            ]
        )
        for index, name, key, err, swap in outputs:
            panels.append(
                panel(
                    title + ": " + name,
                    x,
                    data[prefix + "_" + key],
                    src[index],
                    dst[index],
                    "Volume (Å³/formula unit)",
                    name + " (GPa)",
                    "GPa",
                    swap_axes=swap,
                    sigma=data[prefix + "_" + err],
                    x_sigma=data[prefix + "_vs"],
                    curves=[
                        (f"{temperature:g} K isotherm", grid, a[index], b[index])
                        for temperature, (a, b) in isotherms
                    ],
                )
            )
    return DiagnosticComparison(
        panels,
        "Replay of the full unweighted seven-parameter joint audit: 140 Sun P–V–T rows and 34 Gréaux acoustic rows. Pressure, bulk modulus and shear modulus are separate outputs. Curves are labelled representative isotherms; residuals use every measured V,T. Mixed-temperature points are not connected as one isotherm. Elastic observations and propagated widths use the audit's density/VP/VS reductions. For the registered Peritheos record, the complete stored seven-parameter audit fit supplies the auxiliary shear parameters; the catalog record itself stores only the pressure-producing subset. Source weighting and covariance remain unknown; published non-parity and recommendation findings are unchanged.",
    )


def akins_comparison(record):
    from scripts.reproduce_akins_2004_mgsio3_liquid import (
        AVOGADRO_MOL,
        MOLAR_MASS_G_MOL,
        _load_fit_rows,
        candidate_hugoniot_pressure_gpa,
    )

    rows = _load_fit_rows()
    rho, rho_init, y = (
        rows["shock_density_g_cm3"],
        rows["initial_density_g_cm3"],
        rows["pressure_gpa"],
    )

    def prediction(side, density, initial):
        c = coefficients(record, side)
        rho0 = MOLAR_MASS_G_MOL * 1e24 / (c["V0"] * AVOGADRO_MOL)
        return candidate_hugoniot_pressure_gpa(initial, density, rho0, c["K0"])

    # Default K0S′/gamma0/q/transition energy are the audit's fixed constants.
    predictions = [prediction(s, rho, rho_init) for s in ("published", "refit")]
    curves = []
    for initial in np.unique(rho_init):
        grid = np.linspace(rho.min(), rho.max(), 160)
        curves.append(
            (
                f"Initial density {initial:g} g/cm³",
                grid,
                *[
                    prediction(s, grid, np.full_like(grid, initial))
                    for s in ("published", "refit")
                ],
            )
        )
    return DiagnosticComparison(
        [
            panel(
                "Energy-balanced liquid Hugoniot",
                rho,
                y,
                *predictions,
                "Shock density (g/cm³)",
                "Pressure (GPa)",
                "GPa",
                sigma=rows["pressure_uncertainty_gpa"],
                curves=curves,
            )
        ],
        "Only the three source-selected Akins/enstatite/melt shock states enter this diagnostic. Curves use each recorded initial density and the audit's fixed Mie–Grüneisen/transition-energy reduction. They are theoretical Hugoniots, not the catalog's reference isentrope fitted directly to shock states. Pressure residuals use the source shock pressures. The bounded, ill-conditioned reconstruction and its original qualifications remain applicable.",
    )


def mosenfelder_comparison(record):
    from scripts.reproduce_mosenfelder_2009_mgsio3_post_perovskite import (
        PARAMETER_NAMES,
        arrays,
        shock_pressure,
        static_pressure,
    )

    static, shock = arrays()
    vectors = [
        np.array([coefficients(record, side)[name] for name in PARAMETER_NAMES])
        for side in ("published", "refit")
    ]
    pv = [
        static_pressure(static["volume_a3"], static["temperature_k"], c)
        for c in vectors
    ]
    sh = [
        shock_pressure(
            shock["shock_density_mg_m3"],
            shock["initial_density_mg_m3"],
            shock["transition_energy_j_kg"],
            c,
        )
        for c in vectors
    ]
    grid = np.linspace(static["volume_a3"].min(), static["volume_a3"].max(), 160)
    temperatures = np.unique(static["temperature_k"])
    if len(temperatures) > 5:
        temperatures = temperatures[[0, len(temperatures) // 2, -1]]
    static_curves = [
        (f"{t:g} K isotherm", grid, *[static_pressure(grid, t, c) for c in vectors])
        for t in temperatures
    ]
    density_grid = np.linspace(
        shock["shock_density_mg_m3"].min(), shock["shock_density_mg_m3"].max(), 160
    )
    shock_curves = [
        (
            f"ρinitial={rho:g} g/cm³, ΔE={energy / 1e6:g} MJ/kg",
            density_grid,
            *[shock_pressure(density_grid, rho, energy, c) for c in vectors],
        )
        for rho, energy in np.unique(
            np.column_stack(
                (shock["initial_density_mg_m3"], shock["transition_energy_j_kg"])
            ),
            axis=0,
        )
    ]
    return DiagnosticComparison(
        [
            panel(
                "Static P–V–T observations",
                static["volume_a3"],
                static["pressure_gpa"],
                *pv,
                "Volume (Å³/conventional cell)",
                "Pressure (GPa)",
                "GPa",
                swap_axes=True,
                sigma=static["pressure_uncertainty_gpa"],
                curves=static_curves,
            ),
            panel(
                "Solid PPv shock observations",
                shock["shock_density_mg_m3"],
                shock["pressure_gpa"],
                *sh,
                "Shock density (g/cm³)",
                "Pressure (GPa)",
                "GPa",
                sigma=shock["pressure_uncertainty_gpa"],
                curves=shock_curves,
            ),
        ],
        "Joint audit replay of 48 Guignot P–V–T observations and six PPv shock states, with the audit's specific-heat Cvm convention and starting-material transition energies. Static curves are labelled representative isotherms; shock curves retain each initial density and transition energy. Residuals use each observation's actual state. The displayed vertical pressure RMS is at observed coordinates, distinct from the stored orthogonal objective with latent P,V,T/density adjustments. Shock temperatures are independent diagnostics and were not fitted to the solid PPv inversion.",
    )


def spin_comparison(record):
    from scipy.optimize import brentq

    from scripts.reproduce_solomatova_2016_ferropericlase import (
        CHEN,
        FEI,
        MAO,
        MARQUARDT,
        TABLE1,
        ZHURAVLEV,
        _initial_parameters,
        bm3_energy_ev,
        bm3_pressure,
        mixed_pressure,
        spin_populations,
    )
    from scripts.validate_primary_eos_refits import SOLOMATOVA_RECORD_SOURCES

    key = SOLOMATOVA_RECORD_SOURCES[record["identifier"]]
    specs = {
        "mg061_fei": (
            FEI,
            0.39,
            "sample_lattice_a_angstrom",
            "compression_path",
            "compression",
        ),
        "mg061_zhuravlev": (
            ZHURAVLEV,
            0.39,
            "volume_a3_conventional_cell",
            "used_in_solomatova_2016_fit",
            "1",
        ),
        "fp48": (TABLE1, 0.483, "volume_a3_conventional_cell", None, None),
        "mg075_mao": (MAO, 0.25, "cell_volume_a3", None, None),
        "mg090_marq": (
            MARQUARDT,
            0.10,
            "volume_a3_conventional_cell",
            "used_in_solomatova_2016_fit",
            "1",
        ),
        "mg065_chen": (
            CHEN,
            0.35,
            "volume_a3_conventional_cell",
            "used_in_solomatova_2016_fit",
            "1",
        ),
    }
    path, fraction, volume_key, selection, expected = specs[key]
    rows = [r for r in read_csv(path) if selection is None or r[selection] == expected]
    volume = values(rows, volume_key) ** (
        3 if volume_key == "sample_lattice_a_angstrom" else 1
    )
    pressure = values(rows, "pressure_gpa")
    source = _initial_parameters(key, fraction)
    fitted = source.copy()
    for rid, source_key in SOLOMATOVA_RECORD_SOURCES.items():
        if source_key != key:
            continue
        target = {p["parameter"]: p["refit"] for p in ledger()[rid]["parameters"]}
        if "_hs_" in rid:
            fitted[:3] = [target[n] for n in ("V0", "K0", "K0_prime")]
        else:
            fitted[3:5] = [target[n] for n in ("V0", "K0")]
    # Recover the omitted internal energy offset algebraically from its saved
    # 50% transition pressure. No data fitting or new physical assumption.
    midpoint = ledger()[record["identifier"]]["transition_pressure_50_percent_ls_gpa"]
    vmid = brentq(
        lambda v: (
            0.5
            * (bm3_pressure(v, *fitted[:3]) + bm3_pressure(v, fitted[3], fitted[4], 4))
            - midpoint
        ),
        35,
        83,
    )
    omega_half = brentq(
        lambda w: spin_populations(w) @ np.array([1, 0.5, 0]) - 0.5, -1, 1
    )
    fitted[5] = omega_half - (
        bm3_energy_ev(vmid, fitted[3], fitted[4], 4) - bm3_energy_ev(vmid, *fitted[:3])
    ) / (8 * fraction)
    grid = np.linspace(volume.min(), volume.max(), 160)
    p = panel(
        "Coupled three-level spin-crossover trajectory",
        volume,
        pressure,
        mixed_pressure(volume, source, fraction),
        mixed_pressure(volume, fitted, fraction),
        "Volume (Å³/conventional cell)",
        "Pressure (GPa)",
        "GPa",
        swap_axes=True,
        curves=[
            (
                "Coupled mixed-spin model",
                grid,
                mixed_pressure(grid, source, fraction),
                mixed_pressure(grid, fitted, fraction),
            )
        ],
    )
    return DiagnosticComparison(
        [p],
        "The plotted curve is the complete coupled mixed-spin trajectory; this record's coefficient table describes one HS/LS reference endmember. The source-constrained diagnostic uses both Table 7 endmembers and the Table 2 midpoint. The refit's internal omega0, omitted from the ledger, is recovered algebraically from both stored endmembers and its stored 50% transition pressure. No optimizer is run. Vertical pressure RMS is descriptive; the original audit minimized approximate volume residuals plus priors, and exact MINUTI weights/workflow remain unknown.",
        source_title="Source-constrained coupled model",
    )


def luo_comparison(document, record):
    from scipy.constants import Avogadro

    from peritheos.eos.rt import Vinet
    from peritheos.eos.thermal import SoundVelocityDebyeHelmholtz
    from scripts.validate_primary_eos_refits import _load_rows

    dataset = next(
        d
        for d in document["datasets"]
        if d["identifier"] == "mgo_luo_2023_table1_shock"
    )
    rows = _load_rows(dataset)
    rho, observed_p = values(rows, "density_g_cm3"), values(rows, "pressure_gpa")
    proxy = ledger()[record["identifier"]]["proxy_velocity_fits"]
    models = []
    scale = Avogadro * 1e-25 / float(document["formula_units_per_cell"])
    for side in ("published", "refit"):
        c = coefficients(record, side)
        models.append(
            SoundVelocityDebyeHelmholtz(
                Vinet(c["V0"] * scale, c["K0"], c["K0_prime"]),
                300,
                40.304,
                2,
                proxy["longitudinal_intercept_km_s"],
                proxy["longitudinal_slope_km_s_per_g_cm3"],
                proxy["shear_intercept_km_s"],
                proxy["shear_slope_km_s_per_g_cm3"],
            )
        )
    volume = 40.304 / (10 * rho)
    temperature = [
        np.array(
            [
                m.hugoniot_temperature_from_state(
                    v,
                    p,
                    initial_volume=40.304 / (10 * 3.590),
                    initial_temperature=300,
                    initial_pressure=0,
                )
                for v, p in zip(volume, observed_p)
            ],
            dtype=float,
        )
        for m in models
    ]
    predictions = [
        np.asarray(m.pressure(volume, t), dtype=float)
        for m, t in zip(models, temperature)
    ]
    panels = [
        panel(
            "Luo shock-pressure constraints",
            rho,
            observed_p,
            *predictions,
            "Shock density (g/cm³)",
            "Pressure (GPa)",
            "GPa",
            sigma=values(rows, "pressure_sigma_gpa"),
        )
    ]
    for title, field_name, sigma_name, unit, prediction in (
        (
            "Luo Hugoniot temperatures",
            "temperature_k",
            "temperature_2sigma_k",
            "K",
            temperature,
        ),
        (
            "Luo longitudinal velocities",
            "euler_sound_velocity_km_s",
            "euler_sound_velocity_sigma_km_s",
            "km/s",
            [np.asarray(m.longitudinal_velocity(volume), dtype=float) for m in models],
        ),
    ):
        y = np.array([float(r.get(field_name) or "nan") for r in rows])
        mask = np.isfinite(y)
        sigma = np.array([float(r.get(sigma_name) or "nan") for r in rows])[mask]
        if unit == "K":
            sigma = sigma / 2
        panels.append(
            panel(
                title,
                rho[mask],
                y[mask],
                *[a[mask] for a in prediction],
                "Shock density (g/cm³)",
                "Temperature (K)" if unit == "K" else "Velocity (km/s)",
                unit,
                sigma=sigma,
            )
        )
    return DiagnosticComparison(
        panels,
        "Replay of the existing bounded sensitivity diagnostic: five pressures, four temperatures and three velocities. The Li-derived proxy velocity-density regressions remain fixed; no regression is rerun. Temperature error bars are the audit's one-sigma conversion of the printed two-sigma widths. Model markers retain each measured shock state. This is a qualified proxy-model comparison, not recovery of Luo's unavailable global input set, Kono regressions or calibration covariance.",
        source_title="Published coefficients + audit proxy",
    )


def tange_comparison(record):
    from scripts.reproduce_tange_2009_mgo import (
        CELL_TO_MOLAR_J_PER_BAR,
        ELASTIC_DATASET,
        SHOCK_DATASET,
        THERMAL_DATASETS,
        V0_MOLAR,
        _rows,
        constrained_eos,
        predicted_hugoniot_pressure,
    )

    models = [
        constrained_eos(coefficients(record, side)) for side in ("published", "refit")
    ]
    thermal = [r for filename in THERMAL_DATASETS for r in _rows(filename)]
    volume = values(thermal, "conventional_cell_volume_a3") * CELL_TO_MOLAR_J_PER_BAR
    t = values(thermal, "temperature_k")
    panels = [
        panel(
            "Public thermal-expansion subset",
            t,
            np.zeros(len(t)),
            *[np.asarray(m.pressure(volume, t)) for m in models],
            "Temperature (K)",
            "Ambient pressure constraint (GPa)",
            "GPa",
            observed_label="Zero-pressure expansion constraint",
        )
    ]
    elastic = [r for r in _rows(ELASTIC_DATASET) if r["experimental_path"] != "ambient"]
    ratio = values(elastic, "volume_ratio")
    grid = np.linspace(ratio.min(), ratio.max(), 160)
    panels.append(
        panel(
            "Public Li elastic subset",
            ratio,
            values(elastic, "adiabatic_bulk_modulus_gpa"),
            *[
                np.asarray(m.adiabatic_bulk_modulus(ratio * V0_MOLAR, 300))
                for m in models
            ],
            "V/V0",
            "Adiabatic bulk modulus (GPa)",
            "GPa",
            curves=[
                (
                    "300 K adiabatic modulus",
                    grid,
                    *[
                        np.asarray(m.adiabatic_bulk_modulus(grid * V0_MOLAR, 300))
                        for m in models
                    ],
                )
            ],
        )
    )
    shock = _rows(SHOCK_DATASET)
    rho0, us, up = [
        values(shock, n)
        for n in (
            "initial_density_g_cm3",
            "shock_velocity_km_s",
            "particle_velocity_km_s",
        )
    ]
    ratio = 1 - up / us
    p = rho0 * us * up
    panels.append(
        panel(
            "Public LASL shock subset",
            ratio,
            p,
            *[
                np.array(
                    [predicted_hugoniot_pressure(m, v, pp) for v, pp in zip(ratio, p)]
                )
                for m in models
            ],
            "V/V0",
            "Hugoniot pressure (GPa)",
            "GPa",
        )
    )
    return DiagnosticComparison(
        panels,
        "PARTIAL REPLAY: these panels evaluate both stored coefficient sets only on the public bundled thermal-expansion, Li elastic and LASL shock subsets. They do not represent the complete 164-row approximate global reconstruction: restricted local legacy inputs are not redistributed, the Zha block is a surrogate, and original weights remain unresolved. The two dependent reference parameters are recomputed using the same documented thermodynamic constraints. RMS values apply only to the individually named public subsets.",
        partial_replay=True,
    )


def zhu_expansion_comparison(record):
    from scripts.reproduce_zhu_2025_pressure_standards import (
        THERMAL_RECORDS,
        _debye_energy,
        _gamma_theta,
        _load_thermal_csv,
        _source_excess_pressure,
        _vinet_pressure,
    )

    constants = THERMAL_RECORDS[record["identifier"]]
    prefix = "au" if record["identifier"].startswith("gold") else "pt"
    rows = _load_thermal_csv(f"zhu_2025_{prefix}_zero_pressure_thermal_expansion")
    t, volume = rows["temperature_k"], rows["volume_a3_conventional_cell"]
    ratio = volume / constants["V0"]
    predictions = []
    for side in ("published", "refit"):
        c = {**constants, **coefficients(record, side)}
        gamma, theta = _gamma_theta(ratio, c["gamma0"], c["a"], c["b"], c["theta0"])
        difference = _debye_energy(theta, t, c["n"]) - _debye_energy(theta, 300, c["n"])
        excess = _source_excess_pressure(ratio, t, c) - _source_excess_pressure(
            ratio, 300, c
        )
        predictions.append(
            _vinet_pressure(volume, c)
            + gamma * difference * 1e-9 / (volume * 6.022e-7 / 4)
            + excess
        )
    order = np.argsort(t)
    return DiagnosticComparison(
        [
            panel(
                "Released zero-pressure expansion subset",
                t,
                np.zeros(len(t)),
                *predictions,
                "Temperature (K)",
                "Pressure along input expansion path (GPa)",
                "GPa",
                observed_label="Zero-pressure expansion constraint",
                curves=[
                    (
                        "Stored coefficients on tabulated V(T) path",
                        t[order],
                        *[p[order] for p in predictions],
                    )
                ],
            )
        ],
        "PARTIAL REPLAY: the expansion panel uses every released expansion row and the audit's optimizer/property energy and excess-pressure conventions. It uses the measured V(T) path, not an isotherm. The inconsistent standalone calculators are excluded. Shock rows cannot yet be replayed as the identical fitted residual vector without the unsaved iterative shock-temperature state; they are omitted, not reconstructed as observations. The displayed RMS is for this expansion subset only, not the full robust-bisquare fit.",
        partial_replay=True,
    )
