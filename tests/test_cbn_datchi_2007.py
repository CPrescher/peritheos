import csv
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.materials import Material
from scripts.reproduce_datchi_2007_cbn import load_data, refit_q

RECORD_ID = "boron_nitride_datchi_2007_vinet_mgd_2"
DATASET_ID = "cubic_boron_nitride_datchi_2007_table4_pvt"


def _source_record():
    document = get_material_document("boron_nitride")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    return document, record


def test_datchi_thermal_record_preserves_published_cold_curve_convention():
    document, record = _source_record()

    assert document["phase"] == "cubic zinc-blende"
    assert document["space_group"] == "F-43m"
    assert document["space_group_number"] == 216
    assert document["formula_units_per_cell"] == 4
    assert document["lattice"]["a"] == pytest.approx(3.6152)
    assert {(site["element"], site["wyckoff"]) for site in document["atom_sites"]} == {
        ("B", "4a"),
        ("N", "4c"),
    }
    assert record["eos"]["parameters"] == {
        "V0": 47.2208,
        "K0": 397.0,
        "K0_prime": 3.62,
    }
    assert record["fixed_parameters"] == ["V0", "K0", "K0_prime"]
    assert record["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert record["thermal"]["thermal_pressure_reference"] == "absolute_zero"
    assert record["thermal"]["parameters"] == {
        "Tr": 295.0,
        "theta0": 1700.0,
        "gamma0": 1.04,
        "q": 4.0,
        "n": 2.0,
    }
    assert record["thermal"]["fixed_parameters"] == [
        "Tr",
        "theta0",
        "gamma0",
        "n",
    ]


def test_datchi_absolute_mgd_reproduces_table_vi_300_k_zero_pressure_volume():
    document, _ = _source_record()
    loaded = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]

    # Table VI independently reports V(P=0, 300 K) = 5.9055 A^3/atom.
    assert loaded.volume(0.0, 300.0, check_validity=False) / 8.0 == pytest.approx(
        5.9055, abs=1.0e-5
    )
    assert loaded.eos.thermal_pressure_reference == "absolute_zero"
    assert loaded.eos.thermal_pressure(loaded.eos.rt_eos.V0, 295.0) > 0.0
    assert loaded.eos.thermal_pressure_increment(loaded.eos.rt_eos.V0, 295.0) == (
        pytest.approx(0.0, abs=1.0e-14)
    )


def test_reference_subtracted_interpretation_does_not_reproduce_table_vi_volume():
    document, _ = _source_record()
    loaded = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]
    wrong = MieGruneisenDebye(
        Vinet(loaded.eos.rt_eos.V0, 397.0, 3.62),
        Tr=295.0,
        theta0=1700.0,
        gamma0=1.04,
        q=4.0,
        n=2.0,
        debye_temperature_law="integrated_gruneisen",
        thermal_pressure_reference="reference_temperature",
    )

    wrong_atomic_volume = wrong.volume(0.0, 300.0) / loaded.volume_scale / 8.0
    assert abs(wrong_atomic_volume - 5.9055) > 0.002


def test_datchi_table_iv_is_complete_and_published_curve_has_reported_rmse():
    document, record = _source_record()
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 66
    assert sum(float(row["temperature_k"]) == 295.0 for row in rows) == 38
    assert dataset["used_by_eos_records"] == [
        "boron_nitride_datchi_2007_vinet_1",
        RECORD_ID,
    ]
    assert min(float(row["pressure_gpa"]) for row in rows) == 0.0
    assert max(float(row["pressure_gpa"]) for row in rows) == 162.5
    assert max(float(row["temperature_k"]) for row in rows) == 948.0

    loaded = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]
    cell_volume = np.array([float(row["lattice_a_angstrom"]) ** 3 for row in rows])
    temperature = np.array([float(row["temperature_k"]) for row in rows])
    observed = np.array([float(row["pressure_gpa"]) for row in rows])
    predicted = loaded.pressure(cell_volume, temperature, check_validity=False)
    rmse = float(np.sqrt(np.mean((predicted - observed) ** 2)))

    # Figure 4(b) and Section V report an rms pressure deviation of 0.6 GPa.
    assert rmse == pytest.approx(0.6, abs=0.02)
    assert record["fit_datasets"] == [DATASET_ID]


def test_datchi_q_only_refit_is_reproducible_but_not_a_replacement_record():
    result = refit_q(load_data())

    assert result.success
    assert result.parameters["q"] == pytest.approx(7.34734, abs=1.0e-5)
    assert result.standard_errors["q"] == pytest.approx(1.88832, abs=1.0e-5)
    assert np.sqrt(np.mean(result.residuals**2)) == pytest.approx(0.565419, abs=1.0e-6)
