import csv
import hashlib
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
DOI = "10.1073/pnas.0506811103"
FORMULA = "Mg0.90Fe0.09Al0.005Ca0.005SiO3"
PPV_MATERIAL = "mg090fe009al0005ca0005sio3_post_perovskite"
PV_MATERIAL = "mg090fe009al0005ca0005sio3_bridgmanite"
PPV_DATASET = "mg090fe009al0005ca0005sio3_shieh_2006_figure2_ppv_digitized"
PV_DATASET = "mg090fe009al0005ca0005sio3_shieh_2006_figure2_pv_digitized"
PPV_RESOURCE = "mg090fe009al0005ca0005sio3-shieh-2006-figure2-ppv-digitized.csv"
PV_RESOURCE = "mg090fe009al0005ca0005sio3-shieh-2006-figure2-pv-digitized.csv"


def _rows(resource):
    path = ROOT / "peritheos" / "data" / "datasets" / resource
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _bm3(volume, v0, k0, k0_prime):
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


@pytest.mark.parametrize(
    ("identifier", "space_group", "space_group_number", "volume", "record_count"),
    [
        (PPV_MATERIAL, "Cmcm", 63, 164.9, 2),
        (PV_MATERIAL, "Pbnm", 62, 163.3, 1),
    ],
)
def test_shieh_materials_preserve_identity_and_explicit_topology_proxy(
    identifier, space_group, space_group_number, volume, record_count
):
    document = get_material_document(identifier)
    assert document["formula"] == FORMULA
    assert document["space_group"] == space_group
    assert document["space_group_number"] == space_group_number
    assert document["formula_units_per_cell"] == 4
    assert len(document["eos_records"]) == record_count
    assert document["source"]["topology_proxy"] is True
    assert "no Fe/Al/Ca site occupancies" in document["source"]["topology_proxy_scope"]

    lattice = document["lattice"]
    assert lattice["a"] * lattice["b"] * lattice["c"] == pytest.approx(
        volume, abs=2.0e-8
    )
    assert all(
        "pure-MgSiO3 topology proxy" in site["source_label"]
        for site in document["atom_sites"]
    )

    topology = Counter()
    for site in document["atom_sites"]:
        topology[site["element"]] += site["site_multiplicity"] * site["occupancy"]
    assert topology == {"Mg": 4.0, "Si": 4.0, "O": 12.0}


def test_shieh_published_records_preserve_model_order_and_preference():
    ppv = get_material_document(PPV_MATERIAL)["eos_records"]
    preferred, sensitivity = ppv
    assert preferred["identifier"].endswith("shieh_2006_bm2_1")
    assert preferred["default_for"] == "equilibrium"
    assert preferred["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 164.9, "K0": 219.0},
    }
    assert preferred["implicit_parameters"] == {"K0_prime": 4.0}
    assert preferred["parameter_errors"] == {"V0": 0.6, "K0": 5.0}

    assert sensitivity["identifier"].endswith("shieh_2006_bm3_sensitivity_2")
    assert "default_for" not in sensitivity
    assert sensitivity["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 166.2, "K0": 198.0, "K0_prime": 4.4},
    }
    assert sensitivity["fixed_parameters"] == ["K0_prime"]

    pv = get_material_document(PV_MATERIAL)["eos_records"][0]
    assert pv["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 163.3, "K0": 255.0, "K0_prime": 3.7},
    }
    assert pv["fixed_parameters"] == ["V0", "K0_prime"]

    for record in [*ppv, pv]:
        assert record["reference"]["doi"].lower() == DOI
        assert record["parameter_error_confidence"] == pytest.approx(0.6826894921370859)
        assert record["scientific_validation"]["status"] == "primary_source_validated"
        assert record["pressure_calibration"]["status"] == "partially_resolved"
        assert len(record["pressure_calibration"]["methods"]) == 3


@pytest.mark.parametrize(
    ("material_identifier", "dataset_identifier", "resource", "checksum", "count"),
    [
        (
            PPV_MATERIAL,
            PPV_DATASET,
            PPV_RESOURCE,
            "617855e6a062074a7a7ac763c7786f239a4156ad1953ef7fd7f7961a3185a1e1",
            13,
        ),
        (
            PV_MATERIAL,
            PV_DATASET,
            PV_RESOURCE,
            "32e871601c21aec074bc134426f29791f14879b120cdb96b7b54866bbd8222ac",
            12,
        ),
    ],
)
def test_shieh_figure2_digitization_is_bundled_with_provenance(
    material_identifier, dataset_identifier, resource, checksum, count
):
    dataset = get_material_document(material_identifier)["datasets"][0]
    assert dataset["identifier"] == dataset_identifier
    assert dataset["resource"]["sha256"] == checksum
    path = ROOT / "peritheos" / "data" / "datasets" / resource
    assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum

    rows = _rows(resource)
    assert len(rows) == count
    assert all(row["fit_included"] == "1" for row in rows)
    assert all(
        float(row["pressure_digitization_uncertainty_gpa"]) <= 0.50 for row in rows
    )
    assert dataset["provenance"]["source_pdf_sha256"] == (
        "a52430c01097b8cb00eef8425803227e5f5f0336dd4c36090ccb8fcdc7d84849"
    )
    assert dataset["provenance"]["source_figure_crop_sha256"] == (
        "8eeb03d09c6592d4c54db25695456f1917ba458268c6dd61f00c845f68e93530"
    )
    assert dataset["provenance"]["type"] == "digitized_from_figure"


@pytest.mark.parametrize(
    ("material_identifier", "record_index", "resource", "expected_rmse"),
    [
        (PPV_MATERIAL, 0, PPV_RESOURCE, 1.2951501890),
        (PPV_MATERIAL, 1, PPV_RESOURCE, 1.3684713055),
        (PV_MATERIAL, 0, PV_RESOURCE, 1.2241948029),
    ],
)
def test_shieh_published_curves_reproduce_digitized_figure2(
    material_identifier, record_index, resource, expected_rmse
):
    document = get_material_document(material_identifier)
    stored = document["eos_records"][record_index]
    material = Material.from_eosmat(document, record_identifiers=[stored["identifier"]])
    record = material.get_eos_record(stored["identifier"])
    rows = _rows(resource)
    pressures = np.array([float(row["pressure_gpa"]) for row in rows])
    volumes = np.array([float(row["volume_a3_conventional_cell"]) for row in rows])
    residuals = np.asarray(record.pressure(volumes)) - pressures
    assert math.sqrt(float(np.mean(residuals**2))) == pytest.approx(
        expected_rmse, abs=5.0e-9
    )


def test_shieh_independent_refits_have_parameter_parity():
    ppv_rows = _rows(PPV_RESOURCE)
    ppv_pressure = np.array([float(row["pressure_gpa"]) for row in ppv_rows])
    ppv_volume = np.array(
        [float(row["volume_a3_conventional_cell"]) for row in ppv_rows]
    )
    ppv_document = get_material_document(PPV_MATERIAL)
    for index, expected in enumerate(
        [(164.78775093, 219.92599574), (165.64694843, 201.69342250)]
    ):
        stored = ppv_document["eos_records"][index]
        parameters = stored["eos"]["parameters"]
        k0_prime = stored.get("implicit_parameters", {}).get(
            "K0_prime", parameters.get("K0_prime")
        )
        result = least_squares(
            lambda free: _bm3(ppv_volume, free[0], free[1], k0_prime) - ppv_pressure,
            [parameters["V0"], parameters["K0"]],
        )
        assert result.x == pytest.approx(expected, abs=1.0e-6)
        errors = stored["parameter_errors"]
        assert abs(result.x[0] - parameters["V0"]) < errors["V0"]
        assert abs(result.x[1] - parameters["K0"]) < errors["K0"]

    pv_rows = _rows(PV_RESOURCE)
    pv_pressure = np.array([float(row["pressure_gpa"]) for row in pv_rows])
    pv_volume = np.array([float(row["volume_a3_conventional_cell"]) for row in pv_rows])
    pv = get_material_document(PV_MATERIAL)["eos_records"][0]
    parameters = pv["eos"]["parameters"]
    result = least_squares(
        lambda free: (
            _bm3(
                pv_volume,
                parameters["V0"],
                free[0],
                parameters["K0_prime"],
            )
            - pv_pressure
        ),
        [parameters["K0"]],
    )
    assert result.x[0] == pytest.approx(253.83836311, abs=2.0e-7)
    assert abs(result.x[0] - parameters["K0"]) < pv["parameter_errors"]["K0"]


def test_shieh_audit_disposes_all_nine_litcurate_candidates():
    audit = (
        ROOT
        / "docs"
        / "literature-reproductions"
        / "shieh-2006-natural-pyroxene-pv-ppv.md"
    ).read_text(encoding="utf-8")
    candidate_ids = {
        "litcurate_40925dc9160b7bdd",
        "litcurate_a84d5bbee1a84283",
        "litcurate_c2542f0b2490af33",
        "litcurate_b75920864fe802c4",
        "litcurate_9a4875b670cf3c0e",
        "litcurate_a439638bf2273241",
        "litcurate_d32fda0dff76cc71",
        "litcurate_e7fdbac2dda84b1e",
        "litcurate_766ca79517d1bb73",
    }
    assert all(audit.count(identifier) == 1 for identifier in candidate_ids)
