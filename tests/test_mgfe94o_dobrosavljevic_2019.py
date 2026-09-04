"""Primary-source checks for Dobrosavljevic et al. (2019) Mw94."""

import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document
from peritheos.eos.rt import BM3

ROOT = Path(__file__).resolve().parents[1]
DATASET_SHA256 = "4392a86c0f923eac8ccd6464422024e40dfed44f51617281b62e5d84a35847c1"
B1_ID = "mgfe94o_b1_dobrosavljevic_2019_bm3_1"
RHOMBOHEDRAL_ID = "mgfe94o_rhombohedral_dobrosavljevic_2019_bm3_1"


def _document_dataset_and_rows(identifier="mgfe94o_b1"):
    document = get_material_document(identifier)
    dataset = document["datasets"][0]
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, dataset, path, rows


def _minuti_refit(rows, selection_column, atom_scale, initial, prior=None):
    """Reproduce MINUTI's iteratively weighted volume-residual objective."""
    selected = [row for row in rows if row[selection_column] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in selected])
    pressure_sigma = np.array(
        [float(row["pressure_uncertainty_gpa"]) for row in selected]
    )
    volume = atom_scale * np.array(
        [float(row["atomic_volume_angstrom3"]) for row in selected]
    )
    volume_sigma = atom_scale * np.array(
        [float(row["atomic_volume_uncertainty_angstrom3"]) for row in selected]
    )
    parameters = np.array(initial, dtype=float)

    for _ in range(100):
        model = BM3(V0=parameters[0], K0=parameters[1], K0_prime=parameters[2])
        model_volume = np.asarray(model.volume(pressure), dtype=float)
        bulk_modulus = np.asarray(model.bulk_modulus(model_volume), dtype=float)
        sigma_effective = np.sqrt(
            volume_sigma**2 + (model_volume / bulk_modulus * pressure_sigma) ** 2
        )

        def residual(candidate):
            candidate_model = BM3(
                V0=candidate[0], K0=candidate[1], K0_prime=candidate[2]
            )
            candidate_volume = np.asarray(candidate_model.volume(pressure), dtype=float)
            values = (volume - candidate_volume) / sigma_effective
            if prior is not None:
                prior_value, prior_sigma = prior
                # MINUTI 2 weights each prior by N/n for N rows and n parameters.
                values = np.append(
                    values,
                    np.sqrt(len(selected) / 3.0)
                    * (candidate[2] - prior_value)
                    / prior_sigma,
                )
            return values

        optimization = least_squares(
            residual,
            parameters,
            bounds=([0.8 * initial[0], 50.0, 1.0], [1.2 * initial[0], 400.0, 8.0]),
            x_scale="jac",
            xtol=1.0e-14,
            ftol=1.0e-14,
            gtol=1.0e-14,
        )
        if np.allclose(optimization.x, parameters, rtol=1.0e-12, atol=1.0e-12):
            return optimization.x
        parameters = optimization.x

    raise AssertionError("MINUTI-style effective-variance refit did not converge")


def test_dobrosavljevic_table_s1_transcription_selection_and_checksum():
    _, dataset, path, rows = _document_dataset_and_rows()

    assert hashlib.sha256(path.read_bytes()).hexdigest() == DATASET_SHA256
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == 28
    assert [int(row["source_order"]) for row in rows] == list(range(1, 29))
    assert [row["measurement_id"] for row in rows[7:9]] == ["P8a", "P8b"]
    assert sum(row["phase_branch"] == "B1" for row in rows) == 16
    assert sum(row["phase_branch"] == "rhombohedral" for row in rows) == 12
    assert sum(int(row["used_in_b1_fit"]) for row in rows) == 14
    assert sum(int(row["used_in_rhombohedral_fit"]) for row in rows) == 12
    assert [row["measurement_id"] for row in rows if row["exclusion_reason"]] == [
        "P14",
        "P15",
    ]
    assert rows[0]["pressure_gpa"] == "1.5"
    assert rows[0]["atomic_volume_angstrom3"] == "9.751"
    assert rows[-1]["pressure_gpa"] == "88.7"
    assert rows[-1]["atomic_volume_angstrom3"] == "6.922"


def test_dobrosavljevic_phase_records_preserve_composition_and_volume_bases():
    b1, _, _, _ = _document_dataset_and_rows("mgfe94o_b1")
    rhombohedral, _, _, _ = _document_dataset_and_rows("mgfe94o_rhombohedral")

    assert b1["formula"] == rhombohedral["formula"] == "Mg0.058Fe0.942O"
    assert b1["space_group"] == "Fm-3m"
    assert b1["formula_units_per_cell"] == 4
    assert "space_group" not in rhombohedral
    assert rhombohedral["formula_units_per_cell"] == 1
    assert rhombohedral["atom_sites"] == []

    b1_parameters = b1["eos_records"][0]["eos"]["parameters"]
    rhombohedral_parameters = rhombohedral["eos_records"][0]["eos"]["parameters"]
    assert b1_parameters["V0"] == pytest.approx(8.0 * 9.860)
    assert rhombohedral_parameters["V0"] == pytest.approx(2.0 * 9.59)
    assert b1["eos_records"][0]["parameter_errors"]["V0"] == pytest.approx(8.0 * 0.007)
    assert rhombohedral["eos_records"][0]["parameter_errors"]["V0"] == (
        pytest.approx(2.0 * 0.10)
    )


def test_dobrosavljevic_published_phase_records_execute():
    b1 = Material.from_eosmat(get_material_document("mgfe94o_b1")).default_record()
    rhombohedral = Material.from_eosmat(
        get_material_document("mgfe94o_rhombohedral")
    ).default_record()

    assert b1.identifier == B1_ID
    assert rhombohedral.identifier == RHOMBOHEDRAL_ID
    assert b1.reference_temperature == rhombohedral.reference_temperature == 300.0

    # Independent non-reference-state checks against selected Table S1 rows.
    assert b1.pressure(8.0 * 9.303) == pytest.approx(9.7, abs=0.45)
    assert b1.pressure(8.0 * 8.757) == pytest.approx(23.3, abs=0.35)
    assert rhombohedral.pressure(2.0 * 7.963) == pytest.approx(46.4, abs=1.0)
    assert rhombohedral.pressure(2.0 * 6.922) == pytest.approx(88.7, abs=0.75)
    for record, volumes in (
        (b1, (8.0 * 9.751, 8.0 * 8.757)),
        (rhombohedral, (2.0 * 8.337, 2.0 * 6.922)),
    ):
        for volume in volumes:
            assert record.volume(record.pressure(volume)) == pytest.approx(volume)


def test_dobrosavljevic_minuti_refits_recover_both_published_branches():
    _, _, _, rows = _document_dataset_and_rows()
    b1_refit = _minuti_refit(
        rows,
        "used_in_b1_fit",
        8.0,
        [78.88, 155.3, 3.79],
        prior=(3.8, 0.3),
    )
    rhombohedral_refit = _minuti_refit(
        rows,
        "used_in_rhombohedral_fit",
        2.0,
        [19.18, 217.0, 2.06],
    )

    assert b1_refit == pytest.approx(
        [78.88845145, 155.15330821, 3.78834697], abs=2.0e-6
    )
    assert rhombohedral_refit == pytest.approx(
        [19.11914777, 222.59443231, 1.99412734], abs=2.0e-6
    )

    for identifier, refit in (
        ("mgfe94o_b1", b1_refit),
        ("mgfe94o_rhombohedral", rhombohedral_refit),
    ):
        published = get_material_document(identifier)["eos_records"][0]
        parameters = published["eos"]["parameters"]
        errors = published["parameter_errors"]
        for index, name in enumerate(("V0", "K0", "K0_prime")):
            assert abs(refit[index] - parameters[name]) <= errors[name]
        assert published["scientific_validation"]["independent_refit"]["result"] == (
            "parity"
        )


def test_dobrosavljevic_fit_protocol_calibration_and_spin_scope_are_explicit():
    for identifier in ("mgfe94o_b1", "mgfe94o_rhombohedral"):
        record = get_material_document(identifier)["eos_records"][0]
        assert record["fixed_parameters"] == []
        assert record["fit_provenance"]["refined_parameters"] == [
            "V0",
            "K0",
            "K0_prime",
        ]
        assert record["parameter_covariance"] is None
        assert record["pressure_calibration"]["status"] == "resolved"
        assert (
            record["pressure_calibration"]["methods"][0]["reference_calibration_record"]
            == "ruby_dewaele_2008"
        )
        exclusions = record["scientific_validation"]["excluded_branches"]
        assert any("spin" in item["branch"] for item in exclusions)

    b1 = get_material_document("mgfe94o_b1")["eos_records"][0]
    rhombohedral = get_material_document("mgfe94o_rhombohedral")["eos_records"][0]
    assert b1["fit_provenance"]["priors"] == {
        "K0_prime": {"value": 3.8, "uncertainty": 0.3}
    }
    assert rhombohedral["fit_provenance"]["priors"] == {}
    assert b1["experimental_pressure_range_gpa"] == [1.5, 23.3]
    assert b1["validity"]["pressure_gpa"] == [0.0, 23.3]
    assert rhombohedral["validity"]["pressure_gpa"] == [34.5, 88.7]
