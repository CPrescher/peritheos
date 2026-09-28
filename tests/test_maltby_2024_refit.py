"""Research-trial identities, source discrepancy, and held-out provenance."""

import hashlib
import json
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import numpy as np
import pytest

from peritheos.eos.experimental.maltby_2024 import (
    Maltby2024Parameters,
    Maltby2024Published,
    Maltby2024Trial,
    fcc_shells,
)
from scripts.audit_maltby_2024_refit import (
    SI_COUNTS,
    load_observations,
    metrics,
    prediction,
)

ROOT = Path(__file__).resolve().parents[1]
REPORT = json.loads((ROOT / "docs/data/argon-maltby-2024-refit.json").read_text())


def test_full_si_table_comparison_records_real_geometric_discrepancy():
    actual = dict(zip(*fcc_shells(64)))
    assert [i for i, n in enumerate(SI_COUNTS, 1) if actual.get(i, 0) != n] == [14]
    assert SI_COUNTS[13] == 48
    assert REPORT["conventions"]["shell_table_discrepancies"] == [
        {"squared_distance": 14, "printed_count": 48, "geometric_count": 0}
    ]
    # A second geometric construction: conventional cubic cells with an fcc basis.
    counts = {}
    for i in range(-6, 7):
        for j in range(-6, 7):
            for k in range(-6, 7):
                for a, b, c in [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)]:
                    squared = 2 * ((i + a) ** 2 + (j + b) ** 2 + (k + c) ** 2)
                    if 0 < squared <= 64:
                        counts[squared] = counts.get(squared, 0) + 1
    assert counts == actual


def test_trial_does_not_change_published_parameters():
    original = Maltby2024Published(64)
    changed = replace(original.parameters, epsilon_k=130, rmin_a=3.81)
    trial = Maltby2024Trial(64, changed)
    assert trial.pressure(2.397, 70) != original.pressure(2.397, 70)
    assert original.parameters == Maltby2024Parameters()
    assert original.pressure(2.397, 70) == pytest.approx(
        -0.004225979767987841, abs=1e-14
    )
    with pytest.raises(FrozenInstanceError):
        trial.parameters.epsilon_k = 140
    with pytest.raises(AttributeError):
        trial.parameters = original.parameters


@pytest.mark.parametrize(
    "volume,temperature", [(1.2, 0), (1.8, 10), (2.397, 70), (1.8, 300)]
)
def test_changed_parameters_retain_helmholtz_identity_and_inverse(volume, temperature):
    # Exercise non-default potential, zero-point, vibrational and cold terms.
    params = replace(
        Maltby2024Parameters(),
        epsilon_k=130,
        alpha_r=14.5,
        rmin_a=3.81,
        lambda_k_a9=325000,
        z1_k=141,
        z3=0.7,
        z4_cm3=19.8,
        theta_d0_k=93,
        a_d_k=10.5,
        vref_cm3=22.7,
        gamma_d0=2.6,
        q_d=0.3,
        c2=-1.7,
        c4=-28,
        a0=0.027,
        theta2_k=46,
        gamma2=3.2,
        b1=-0.00045,
    )
    model = Maltby2024Trial(64, params)
    h = volume * 1e-5
    pressure = (
        -(
            model.molar_helmholtz_energy(volume + h, temperature)
            - model.molar_helmholtz_energy(volume - h, temperature)
        )
        / (2 * h)
        * 1e-4
    )
    assert model.pressure(volume, temperature) == pytest.approx(pressure, abs=1e-7)
    assert model.volume(pressure, temperature) == pytest.approx(volume, abs=1e-8)
    x = model.sigma_angstrom / params.rmin_a
    assert 6 * np.exp(params.alpha_r * (1 - x)) == pytest.approx(
        params.alpha_r / x**6, rel=1e-12
    )


@pytest.mark.parametrize(
    "changes",
    [
        {"epsilon_k": np.nan},
        {"epsilon_k": 0},
        {"alpha_r": 6},
        {"alpha_r": 7.5},
        {"theta_d0_k": 5},
        {"a_d_k": -1},
        {"b2": -1},
        {"a0": 1},
        {"a1": -1},
        {"lambda_k_a9": -1},
    ],
)
def test_invalid_trial_coefficients(changes):
    with pytest.raises(ValueError):
        replace(Maltby2024Parameters(), **changes)


def test_refit_report_predictions_are_at_original_observations_without_run_leakage():
    source = ROOT / "peritheos/data/datasets/argon-fcc-dewaele-2021-supplement.csv"
    assert (
        hashlib.sha256(source.read_bytes()).hexdigest() == REPORT["source_data_sha256"]
    )
    originals, _ = load_observations()
    by_id = {r["row_id"]: r for r in originals}
    for experiment in REPORT["experiments"].values():
        rows = experiment["predictions"]
        ids = {r["row_id"] for r in rows}
        test_union = []
        for row in rows:
            for key, value in by_id[row["row_id"]].items():
                if isinstance(value, float):
                    # Derived observation coordinates can differ by one or two
                    # ULPs across platforms; source bytes and row IDs stay exact.
                    assert row[key] == pytest.approx(value, rel=1e-14, abs=1e-14)
                else:
                    assert row[key] == value
        for fold in experiment["folds"]:
            train_ids, test_ids = set(fold["train_rows"]), set(fold["test_rows"])
            assert not train_ids & test_ids
            assert train_ids | test_ids == ids
            assert all(by_id[s]["run"] == fold["held_out_run"] for s in test_ids)
            assert all(by_id[s]["run"] != fold["held_out_run"] for s in train_ids)
            test_union.extend(test_ids)
            model = Maltby2024Trial(
                64, replace(Maltby2024Parameters(), **fold["fit"]["fitted_parameters"])
            )
            held_rows = [r for r in rows if r["row_id"] in test_ids]
            expected = prediction(model, held_rows)
            assert expected == pytest.approx(
                [r["out_of_fold_pressure_gpa"] for r in held_rows], abs=1e-12
            )
        assert len(test_union) == len(set(test_union)) == len(rows)
        expected = metrics(
            np.array([r["out_of_fold_pressure_gpa"] for r in rows]), rows
        )
        assert expected == pytest.approx(experiment["grouped_validation"], abs=1e-12)
    assert REPORT["experiments"]["all_130_equal_run"]["baseline"]["count"] == 130
    assert REPORT["experiments"]["source_informed_equal_run"]["baseline"]["count"] == 95
    assert REPORT["status"] == "diagnostic_refit_not_promoted"
