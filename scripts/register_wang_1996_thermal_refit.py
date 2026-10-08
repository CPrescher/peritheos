"""Register the optional Wang-table thermal refit and its exact input subset.

Published records and source tables are preserved. Equal pressure weights are
an analyst choice; no reconstruction of the original source weights is claimed.
"""

import copy
import csv
import hashlib
import json

import numpy as np
import scipy

from scripts.reproduce_wang_1996_casio3 import DATASET_ID as SOURCE_ID
from scripts.reproduce_wang_1996_casio3 import ROOT, load_table, thermal_fit

RECORD_ID = "ca_perovskite_wang_1996_unweighted_bm3_linear_thermal_refit"
DATASET_ID = "ca_perovskite_wang_1996_selected64_thermal_refit"
FILENAME = "ca-perovskite-wang-1996-selected64-thermal-refit.csv"
NOTE = (
    "Independent Peritheos BM3 + linear thermal-pressure refit to 64 experimental "
    "Table 1 states. The two source-flagged sub-2-GPa amorphizing states are "
    "excluded. Equal pressure weights are an analyst choice, not recovered "
    "Wang weights. Conditional residual-scaled standard errors omit coordinate "
    "uncertainties, correlations, pressure-scale systematics and fixed-parameter "
    "uncertainty. This is not a reproduction of the published thermal fit or "
    "validation as a pressure standard."
)


def ledger_outcome(record):
    """Numerically reproduce this registered fit using its bundled 64 rows."""
    _, rows = load_table(DATASET_ID)
    result = thermal_fit(rows, selection="exclude_two_low_64")
    if result["observations"] != 64:
        raise AssertionError("Expected the registered 64-row selection")
    names = {
        "rt_eos.V0": "V0_a3",
        "rt_eos.K0": "K0_gpa",
        "alpha_KT": "dP_dT_v_gpa_per_k",
    }
    stored = {
        "rt_eos.V0": record["eos"]["parameters"]["V0"],
        "rt_eos.K0": record["eos"]["parameters"]["K0"],
        "alpha_KT": record["thermal"]["parameters"]["alpha_KT"],
    }
    comparisons = []
    for name, value in stored.items():
        fitted = result["parameters"][names[name]]
        comparisons.append(
            dict(
                parameter=name,
                published=value,
                published_error=None,
                refit=fitted,
                refit_error=result["standard_errors"][names[name]],
                difference=fitted - value,
                relative_difference=(fitted - value) / value,
                within_combined_2sigma=None,
                similar=bool(np.isclose(value, fitted, rtol=1e-6, atol=1e-10)),
            )
        )
    if not all(p["similar"] for p in comparisons):
        raise AssertionError("Registered Wang thermal refit is stale")
    return dict(
        status="parity",
        parity_basis="registered_peritheos_refit_reproduced",
        reproduction_status="refit_reproduced",
        original_publication_reproduction_status="not_reproduced",
        dataset_identifiers=[DATASET_ID],
        observations=64,
        fit_kind="selected_bm3_linear_thermal_pressure",
        objective="equal pressure weights",
        absolute_sigma=False,
        parameters=comparisons,
        rmse_gpa=result["pressure_rmse_gpa"],
        solver_success=result["solver_success"],
        qualification=NOTE,
    )


def register():
    source, rows = load_table(SOURCE_ID)
    fit = thermal_fit(rows, selection="exclude_two_low_64")
    selected = [r for r in rows if int(r["source_order"]) in fit["source_orders"]]
    assert len(selected) == 64
    assert {int(r["source_order"]) for r in rows} - set(fit["source_orders"]) == {
        33,
        34,
    }
    path = ROOT / "peritheos/data/datasets" / FILENAME
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    material_path = ROOT / "peritheos/data/materials/ca_perovskite.eosmat"
    doc = json.loads(material_path.read_text())
    published = next(
        r
        for r in doc["eos_records"]
        if r["identifier"] == "ca_perovskite_wang_1996_preferred_bm3"
    )
    p, e = fit["parameters"], fit["standard_errors"]
    pressure_range = [
        min(float(r["pressure_gpa"]) for r in selected),
        max(float(r["pressure_gpa"]) for r in selected),
    ]
    temperature_range = [
        min(float(r["temperature_k"]) for r in selected),
        max(float(r["temperature_k"]) for r in selected),
    ]
    validation = dict(
        status="primary_source_validated",
        audit_date="2026-10-08",
        note=NOTE,
        refit_status="independent_selected_data_refit",
        reproduction_status="refit_reproduced",
        original_publication_reproduction_status="not_reproduced",
        primary_source_check=dict(
            access_url=source["source_url"],
            doi="10.1029/95jb03254",
            locations=["Table 1", "Equation (1)", "Thermal Pressure"],
            finding="Published Table 1 observations and BM3 conventions checked; coefficients are an independent equal-weight Peritheos refit with constant thermal-pressure slope.",
        ),
        primary_data_check=dict(
            status="bundled",
            audit_date="2026-10-08",
            dataset_identifiers=[DATASET_ID],
            finding="Exact 64-row manifest retains source-row/run identities. Rows 33 and 34 are excluded without altering the 66-row source table. Runs 3/4/5 originate in Wang–Weidner (1994).",
        ),
        usage_recommendation="Optional comparison within sampled coverage; range bounds describe an observation envelope, not complete P–T coverage. Preserve the source pressure-scale and weighting qualifications. Do not extrapolate to ambient stability or present as Wang's published thermal coefficients.",
    )
    record = dict(
        identifier=RECORD_ID,
        label="Peritheos refit of Wang (1996), BM3 + linear thermal pressure (64 points)",
        record_kind="refit",
        default=False,
        determination_method="experimental",
        equation_kind="thermal",
        reference=copy.deepcopy(published["reference"]),
        derived_from_record=published["identifier"],
        source_lineage=[
            dict(
                role="experimental_data",
                citation="Wang et al. (1996), Table 1; runs 3/4/5 inherited from Wang–Weidner (1994)",
                doi="10.1029/95JB03254",
            )
        ],
        eos=dict(
            type="BM3",
            model="birch_murnaghan_3",
            parameters=dict(V0=p["V0_a3"], K0=p["K0_gpa"], K0_prime=4.8),
        ),
        parameter_errors=dict(V0=e["V0_a3"], K0=e["K0_gpa"], K0_prime=None),
        parameter_error_confidence=None,
        parameter_covariance=dict(
            parameter_order=["rt_eos.V0", "rt_eos.K0", "alpha_KT"],
            matrix=fit["covariance"],
        ),
        fixed_parameters=["K0_prime"],
        temperature_ref=300.0,
        thermal=dict(
            type="LinearThermalPressure",
            model="linear_thermal_pressure",
            parameters=dict(Tr=300.0, alpha_KT=p["dP_dT_v_gpa_per_k"]),
            parameter_errors=dict(Tr=None, alpha_KT=e["dP_dT_v_gpa_per_k"]),
            fixed_parameters=["Tr"],
        ),
        volume_basis=dict(
            kind="formula_units", formula_units=1.0, molar_mass_g_mol=116.162
        ),
        fit_datasets=[DATASET_ID],
        experimental_pressure_range_gpa=pressure_range,
        experimental_temperature_range_k=temperature_range,
        validity=dict(
            pressure_gpa=pressure_range,
            temperature_k=temperature_range,
            notes=[
                "Observed range envelope only; uneven P–T coverage. P=0 at fitted V0 and 300 K is an extrapolated reference, not an ambient stable-phase observation."
            ],
        ),
        pressure_calibration=copy.deepcopy(published["pressure_calibration"]),
        parameter_provenance={
            "determination_method": "Independent fit of experimental Table 1 P–V–T states; K0-prime=4.8 is adopted from Wang's experimental Mao-data reanalysis, and the thermal slope is fitted to these observations.",
            "eos.V0": "Joint unit-weight pressure fit; extrapolated 300 K formula-unit reference volume in Å³.",
            "eos.K0": "Joint unit-weight pressure fit; reference bulk modulus in GPa.",
            "eos.K0_prime": "Fixed at Wang's adopted 4.8, not re-estimated.",
            "thermal.alpha_KT": "Joint unit-weight pressure fit; constant (∂P/∂T)V in GPa/K.",
            "thermal.Tr": "Fixed at 300 K, matching the source thermal-pressure reference.",
        },
        fit_provenance=dict(
            software=dict(
                name="SciPy least_squares / Peritheos",
                version=scipy.__version__,
                version_date="2026-10-08",
            ),
            dataset=DATASET_ID,
            selection=dict(
                source_dataset=SOURCE_ID,
                count=64,
                source_orders=fit["source_orders"],
                excluded_source_orders=[33, 34],
                reason="Source-flagged sub-2-GPa amorphizing states; no values corrected.",
            ),
            objective="Minimize squared pressure residuals with unit weight per observation.",
            refined_parameters=["rt_eos.V0", "rt_eos.K0", "alpha_KT"],
            fixed_parameters=["K0_prime", "Tr"],
            model_assumptions="Constant volume-independent thermal-pressure slope; logarithmic volume and quadratic-temperature terms set to zero. K0-prime fixed at Wang's adopted 4.8; Tr fixed at 300 K.",
            statistics=dict(
                observations=64,
                free_parameters=3,
                degrees_of_freedom=61,
                pressure_rmse_gpa=fit["pressure_rmse_gpa"],
                residual_sum_squares_gpa2=fit["weighted_residual_sum_squares"],
                solver_success=True,
            ),
            covariance_scaling="Full local (J.T J)^-1 scaled by RSS/61; one-standard-error estimates conditional on fixed K0-prime/Tr, independent common-variance pressure residuals and the chosen model. Coordinate uncertainties, shared correlations and pressure-scale systematics omitted.",
            reproduction=dict(
                script="scripts/register_wang_1996_thermal_refit.py",
                diagnostic_script="scripts/reproduce_wang_1996_casio3.py",
                artifact="docs/data/wang-1996-casio3-refit.json",
                artifact_section="thermal_pressure.unweighted_exclude_low_64",
                initial_parameters=[45.58, 233.0, 0.0071],
                solver="scipy.optimize.least_squares",
                max_nfev=10000,
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            ),
        ),
        scientific_validation=validation,
        notes=NOTE,
    )
    dataset = copy.deepcopy(source)
    dataset.update(
        identifier=DATASET_ID,
        description="Exact 64 selected experimental states for the optional Peritheos BM3 + linear thermal-pressure refit; source uncertainties retained but not used as weights.",
        used_by_eos_records=[RECORD_ID],
        resource=dict(
            path="datasets/" + FILENAME,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            media_type="text/csv",
        ),
    )
    dataset["provenance"].update(
        type="selected_source_observations",
        source_dataset=SOURCE_ID,
        source_dataset_sha256=source["resource"]["sha256"],
        selection_script="scripts/register_wang_1996_thermal_refit.py",
        excluded_source_orders=[33, 34],
    )
    dataset["notes"] = (
        "Subset of the full 66-row Table 1 source, with rows 33/34 omitted. It overlaps the full source, the room-temperature subset and inherited Wang–Weidner observations; do not combine these as independent measurements. "
        + NOTE
    )
    doc["eos_records"] = [
        r for r in doc["eos_records"] if r["identifier"] != RECORD_ID
    ] + [record]
    doc["datasets"] = [d for d in doc["datasets"] if d["identifier"] != DATASET_ID] + [
        dataset
    ]
    material_path.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    register()
