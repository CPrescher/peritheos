import csv
import hashlib
import math
import re
from collections import Counter
from pathlib import Path

import pytest

from peritheos import Material, load_eosmat, validate_eosmat_document
from peritheos.eosmat import validate_pressure_calibration_references
from scripts.reproduce_kubo_2000_aluminous_bridgmanite import (
    refit_model,
    reproduce_model,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIAL_PATH = (
    ROOT / "peritheos" / "data" / "materials" / "mg09al02si09o3_bridgmanite.eosmat"
)
DATA_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "mg09al02si09o3-bridgmanite-kubo-2000-table1-compression.csv"
)
DATASET_ID = "mg09al02si09o3_bridgmanite_kubo_2000_table1_compression"
FIXED_ID = "mg09al02si09o3_bridgmanite_kubo_2000_bm3_1"
FREE_ID = "mg09al02si09o3_bridgmanite_kubo_2000_bm3_2"
CHECKSUM = "1d18597504e34fe54cdbb839a29f63f7118324c046351c3174d7e2d923f7eb9f"


def _load():
    document = load_eosmat(MATERIAL_PATH)
    validate_eosmat_document(document)
    material = Material.from_eosmat(document)
    records = {record.identifier: record for record in material.eos_records}
    sources = {record["identifier"]: record for record in document["eos_records"]}
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, records, sources, rows


def test_kubo_material_identity_structure_and_proxy_are_explicit():
    document, _, sources, _ = _load()

    assert document["identifier"] == "mg09al02si09o3_bridgmanite"
    assert document["formula"] == "Mg0.9Al0.2Si0.9O3"
    assert document["space_group"] == "Pbnm"
    assert document["space_group_number"] == 62
    assert document["formula_units_per_cell"] == 4
    assert document["lattice"] == {
        "a": 4.774,
        "b": 4.941,
        "c": 6.935,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
    }
    assert math.prod(document["lattice"][axis] for axis in "abc") == pytest.approx(
        163.6, abs=0.02
    )
    assert "topology proxy" in document["cell_contents"]
    assert "do not encode the specimen composition" in document["cell_contents"]

    proxy_contents = Counter()
    for site in document["atom_sites"]:
        multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
        proxy_contents[site["element"]] += multiplicity * site["occupancy"]
    assert proxy_contents == {"Mg": 4.0, "Si": 4.0, "O": 12.0}
    assert set(sources) == {FIXED_ID, FREE_ID}


def test_kubo_fixed_fit_is_default_and_free_fit_is_secondary():
    _, records, sources, _ = _load()
    fixed = sources[FIXED_ID]
    free = sources[FREE_ID]

    assert fixed["default_for"] == "equilibrium"
    assert "default_for" not in free
    assert fixed["reference"]["doi"] == "10.2183/pjab.76.103"
    assert fixed["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 163.6, "K0": 225.5, "K0_prime": 4.0},
    }
    assert fixed["parameter_errors"] == {"V0": 0.4, "K0": 1.2}
    assert fixed["fixed_parameters"] == ["K0_prime"]
    assert fixed["parameter_covariance"] is None
    assert fixed["parameter_error_confidence"] is None
    assert fixed["author_fit_preference"].startswith("preferred;")

    assert free["eos"]["parameters"] == {
        "V0": 163.6,
        "K0": 215.4,
        "K0_prime": 7.2,
    }
    assert free["parameter_errors"] == {
        "V0": 0.4,
        "K0": 4.4,
        "K0_prime": 1.4,
    }
    assert free["fixed_parameters"] == []
    assert free["author_fit_preference"].startswith("nonpreferred;")

    for identifier in (FIXED_ID, FREE_ID):
        record = records[identifier]
        assert record.pressure(163.6, 300.0) == pytest.approx(0.0, abs=1.0e-13)
        pressure = record.pressure(159.9, 300.0, check_validity=True)
        assert record.volume(pressure, 300.0, check_validity=True) == pytest.approx(
            159.9, abs=1.0e-9
        )
        with pytest.raises(ValueError, match="outside the published calibration"):
            record.volume(9.0, 300.0, check_validity=True)


def test_kubo_table_i_is_complete_checksummed_and_lossless():
    document, _, _, rows = _load()
    dataset = document["datasets"][0]

    assert hashlib.sha256(DATA_PATH.read_bytes()).hexdigest() == CHECKSUM
    assert dataset["resource"]["sha256"] == CHECKSUM
    assert dataset["identifier"] == DATASET_ID
    assert dataset["used_by_eos_records"] == [FIXED_ID, FREE_ID]
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == 10
    assert [int(row["source_order"]) for row in rows] == list(range(1, 11))
    assert Counter(row["series"] for row in rows) == {"ap": 4, "ku": 6}
    assert all(row["used_in_published_fit"] == "1" for row in rows)
    assert max(float(row["pressure_ruby_gpa"]) for row in rows) == 8.71
    assert sum(bool(row["pressure_au_gpa"]) for row in rows) == 6

    assert rows[0] == {
        "source_order": "1",
        "run_id": "ap102",
        "series": "ap",
        "pressure_ruby_gpa": "5.31",
        "pressure_ruby_uncertainty_gpa": "0.01",
        "pressure_au_gpa": "",
        "pressure_au_uncertainty_gpa": "",
        "lattice_a_angstrom": "4.740",
        "lattice_a_uncertainty_angstrom": "0.005",
        "lattice_b_angstrom": "4.904",
        "lattice_b_uncertainty_angstrom": "0.004",
        "lattice_c_angstrom": "6.878",
        "lattice_c_uncertainty_angstrom": "0.006",
        "volume_a3_conventional_cell": "159.9",
        "volume_uncertainty_a3": "0.4",
        "b_over_a": "1.035",
        "c_over_a": "1.451",
        "run_reference_volume_a3": "163.6",
        "used_in_published_fit": "1",
    }
    assert rows[-1]["run_id"] == "ku117"
    assert rows[-1]["pressure_ruby_gpa"] == "0.00"
    assert rows[-1]["pressure_au_gpa"] == "0.00"
    assert rows[-1]["volume_a3_conventional_cell"] == "163.5"


@pytest.mark.parametrize(
    ("name", "rmse", "maximum"),
    (
        ("fixed", 0.0834714544, 0.1501695655),
        ("free", 0.0591814632, 0.1270385963),
    ),
)
def test_kubo_published_curves_reproduce_table(name, rmse, maximum):
    result = reproduce_model(name)

    assert result["observations"] == 10
    assert result["series_normalized_pressure_rmse_gpa"] == pytest.approx(
        rmse, abs=5.0e-10
    )
    assert result["series_normalized_max_abs_pressure_residual_gpa"] == (
        pytest.approx(maximum, abs=5.0e-10)
    )


def test_kubo_independent_refits_have_coefficient_parity():
    fixed = refit_model("fixed")
    fixed_axes = refit_model("fixed", use_axis_products=True)
    free = refit_model("free")

    assert fixed["observations"] == 8
    assert fixed["parameters"] == pytest.approx(
        {"K0": 226.86340115, "K0_prime": 4.0}, abs=5.0e-8
    )
    assert fixed_axes["parameters"]["K0"] == pytest.approx(224.83957212, abs=5.0e-8)
    assert free["parameters"] == pytest.approx(
        {"K0": 215.61948074, "K0_prime": 7.53074835}, abs=2.0e-6
    )

    assert abs(fixed["parameters"]["K0"] - 225.5) < 1.2 * 1.2
    assert abs(fixed_axes["parameters"]["K0"] - 225.5) < 1.2
    assert abs(free["parameters"]["K0"] - 215.4) < 4.4
    assert abs(free["parameters"]["K0_prime"] - 7.2) < 1.4


def test_kubo_pressure_calibrations_resolve_and_limits_are_explicit():
    _, _, sources, _ = _load()

    for source in sources.values():
        calibration = source["pressure_calibration"]
        assert calibration["status"] == "partially_resolved"
        assert calibration["methods"][0]["reference_calibration_record"] == (
            "ruby_mao_1986"
        )
        assert calibration["methods"][1]["reference_eos_record"] == (
            "gold_anderson_1989_bm3_1"
        )
        assert calibration["recalculation"]["status"] == (
            "missing_calibrant_observations"
        )
        assert source["fit_datasets"] == [DATASET_ID]
        assert source["scientific_validation"]["status"] == ("primary_source_validated")

    validate_pressure_calibration_references()
