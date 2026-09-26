"""Register the explicitly provisional refit and its exact 49 selected rows.

Run the combined-fit diagnostic first, then this module. Source observations
and the published catalog records are preserved.
"""

import copy
import csv
import hashlib
import json

import scipy

from scripts.fit_argon_errandonea_combined import PREFERRED, REPORT, observations
from scripts.reproduce_argon_errandonea_2006 import ROOT

RECORD_ID = "argon_fcc_ross_1986_errandonea_2006_bm3_refit"
DATASET_ID = "argon_fcc_ross_errandonea_selected49"
FILENAME = "argon-fcc-ross-errandonea-selected49.csv"


def register():
    report = json.loads(REPORT.read_text())
    fit = report["fits"][PREFERRED]["pressure"]
    rows = [r for r in observations() if r["included_in_preferred_fit"]]
    assert fit["count"] == len(rows) == 49
    assert fit["solver_success"] and fit["active_bound_mask"] == [0, 0, 0]
    path = ROOT / "peritheos/data/datasets" / FILENAME
    names = ["source", "source_row", "pressure_gpa", "volume_a3", "temperature_k"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=names, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    material_path = ROOT / "peritheos/data/materials/argon_fcc.eosmat"
    doc = json.loads(material_path.read_text())
    published = next(
        r
        for r in doc["eos_records"]
        if r["identifier"] == "argon_fcc_errandonea_2006_bm3"
    )
    lineage = [
        dict(
            role="experimental_data",
            citation="Ross et al. (1986), JCP 85, 1028–1033, Table I (298 K)",
            doi="10.1063/1.451346",
        ),
        dict(
            role="experimental_data",
            citation="Errandonea et al. (2006), PRB 73, 092106, Figure 5 (300 K)",
            doi="10.1103/PhysRevB.73.092106",
        ),
    ]
    record = dict(
        identifier=RECORD_ID,
        label="Provisional Peritheos BM3 refit: Ross (1986) + Errandonea (2006)",
        record_kind="refit",
        default=False,
        determination_method="experimental",
        reference=copy.deepcopy(published["reference"]),
        source_lineage=lineage,
        eos=dict(type="BM3", model="birch_murnaghan_3", parameters=fit["parameters"]),
        parameter_errors={},
        fixed_parameters=[],
        temperature_ref=300.0,
        experimental_pressure_range_gpa=[
            min(r["pressure_gpa"] for r in rows),
            max(r["pressure_gpa"] for r in rows),
        ],
        experimental_temperature_range_k=[298.0, 300.0],
        fit_datasets=[DATASET_ID],
        fit_provenance=dict(
            software=dict(
                name="SciPy least_squares / Peritheos diagnostic scripts",
                version=scipy.__version__,
                version_date="2026-09-26",
            ),
            dataset=DATASET_ID,
            selection=dict(
                count=49,
                ross_rows=41,
                errandonea_rows=8,
                excluded=[
                    dict(
                        source="Ross Table I",
                        source_row=17,
                        original_pressure_gpa=24.7,
                        reason="Suspected printed pressure error; user-requested exclusion, no corrected value asserted.",
                    ),
                    dict(
                        source="Anderson–Swenson Figure 5 diamond",
                        reason="Cryogenic reference, not room-temperature data.",
                    ),
                ],
                note="Errandonea points are a partial digitization; unresolved overlapping markers are absent. Source archive preserved in docs/data/argon-errandonea-2006-combined-inputs.json.",
            ),
            objective="Minimize sum of squared pressure residuals at the observed volumes, with equal weight per observation.",
            refined_parameters=["V0", "K0", "K0_prime"],
            fixed_parameters=[],
            bounds=report["fit_bounds"],
            statistics=dict(
                count=49,
                pressure_rmse_gpa=fit["pressure_rmse_gpa"],
                volume_rmse_a3=fit["volume_rmse_a3"],
                solver_success=True,
            ),
            reproduction=dict(
                script="scripts/fit_argon_errandonea_combined.py",
                report="docs/data/argon-errandonea-2006-combined-fit.json",
                subset=PREFERRED,
                objective="pressure",
            ),
            uncertainty_note="No calibrated parameter uncertainties or covariance assigned; equal weights are an analyst choice, and fits depend strongly on the residual objective.",
        ),
        parameter_provenance=dict(
            determination_method="Peritheos refit of the explicitly selected 49 experimental observations; these are not published source coefficients.",
            V0="Fitted zero-pressure extrapolated volume in A^3 per four-atom fcc cell.",
            K0="Fitted bulk modulus in GPa.",
            K0_prime="Fitted dimensionless pressure derivative.",
            reference_state="Nominal 300 K; 298 K Ross and 300 K Errandonea pooled without thermal correction. P=0 at extrapolated V0.",
        ),
        pressure_calibration=copy.deepcopy(published["pressure_calibration"]),
        scientific_validation=dict(
            status="primary_source_validated",
            audit_date="2026-09-26",
            refit_status="provisional",
            reproduction_status="refit_reproduced",
            note="Selected-data refit is reproducible; this is not reproduction of the original Errandonea fit or validation as a pressure standard. Source table and digitized points checked; weighting, temperature pooling, and partial recovery remain limitations.",
        ),
        notes="Provisional combined-data refit; not a default EOS or a reproduction of published coefficients. 41 Ross Table I observations at 298 K and 8 resolved Errandonea Fig.5 fcc markers at 300 K. Ross row17 (24.7 GPa) excluded as a suspected source error, retained in source archive. Cryogenic diamond excluded. Equal pressure weights; volume errors and covariance omitted. Strong objective sensitivity; parameter uncertainties not established. Fcc reflections retained through fcc/hcp coexistence. No independent hcp fit and no thermal EOS. Reference DOI identifies an input study, not authorship of these refitted coefficients.",
    )
    record["pressure_calibration"]["notes"] = (
        "Original reported pressures retained. Ross Table I uses ruby A=19.04 Mbar, B=7.665 (Eq.1); Errandonea cites Mao1986 ruby with W cross-check. No raw calibrant readings or pressure recalibration."
    )
    columns = [
        dict(name=n, quantity=q, unit=u, role="value")
        for n, q, u in [
            ("source", "source", "1"),
            ("source_row", "source_row", "1"),
            ("pressure_gpa", "pressure", "GPa"),
            ("volume_a3", "volume", "angstrom^3/conventional_unit_cell"),
            ("temperature_k", "temperature", "K"),
        ]
    ]
    dataset = dict(
        identifier=DATASET_ID,
        kind="pressure_volume_temperature",
        description="Exact 49 selected rows for the provisional combined BM3 refit; source and source-row identities retained.",
        reference=copy.deepcopy(published["reference"]),
        source_lineage=lineage,
        source_location="Ross (1986) Table I excluding row17; Errandonea (2006) Figure 5 eight resolved filled markers.",
        columns=columns,
        resource=dict(
            path="datasets/" + FILENAME,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            media_type="text/csv",
        ),
        used_by_eos_records=[RECORD_ID],
        provenance=dict(
            type="selected_multi_source_observations",
            source_archive="docs/data/argon-errandonea-2006-combined-inputs.json",
            selection_script="scripts/register_argon_combined_refit.py",
            notes="Original excluded values remain in the archive; no fitted-curve or theory samples included.",
        ),
    )
    doc["eos_records"] = [
        r for r in doc["eos_records"] if r["identifier"] != RECORD_ID
    ] + [record]
    doc["datasets"] = [d for d in doc["datasets"] if d["identifier"] != DATASET_ID] + [
        dataset
    ]
    material_path.write_text(json.dumps(doc, indent=1) + "\n")


if __name__ == "__main__":
    register()
