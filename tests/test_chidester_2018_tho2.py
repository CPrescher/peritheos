"""Primary ThO2 data transcription, molar normalization, and source reproduction."""

import csv
import hashlib
import json
from collections import Counter

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, eosmat_schema, get_material, get_material_document
from scripts.reproduce_chidester_2018_tho2 import (
    CELL_PER_MOLAR,
    ERRORS,
    PUBLISHED,
    ROOT,
    ledger_outcome,
    reproduce,
    source_pressure,
)


@pytest.mark.parametrize(
    "material,index,count,selected",
    [("thorianite", 1, 41, 29), ("tho2_cotunnite", 2, 59, 59)],
)
def test_complete_official_rows_and_errors(material, index, count, selected):
    folder = ROOT / "peritheos/data/datasets/chidester_2018_sources"
    manifest = json.loads((folder / "manifest.json").read_text())
    for entry in manifest["files"]:
        assert (
            hashlib.sha256((folder / entry["path"]).read_bytes()).hexdigest()
            == entry["sha256"]
        )
    raw = list(
        csv.reader(
            (folder / f"6212_supp{index}.csv")
            .read_text(encoding="utf-8-sig")
            .splitlines(),
            delimiter="," if index == 1 else "\t",
        )
    )[4:]
    path = ROOT / f"peritheos/data/datasets/chidester_2018_{material}_pvt.csv"
    with path.open() as stream:
        reader = csv.reader(stream)
        next(reader)
        rows = list(reader)
    assert len(rows) == count
    assert [r[:-1] for r in rows] == raw  # Every original field, including errors.
    assert sum(int(r[-1]) for r in rows) == selected
    ds = get_material(material).get_dataset(f"{material}_chidester_2018_table_s{index}")
    assert len(ds) == count
    assert ds.checksum_verified
    assert ds["pressure_standard"][-1] == ("KBr" if index == 1 else "KCl")
    # Independent duplicate volume columns establish cm3/mol of ThO2 and Z=4.
    assert (
        np.max(np.abs(ds["molar_volume_cm3_mol"] * CELL_PER_MOLAR - ds["volume_a3"]))
        < 0.10
    )
    assert (
        ds.columns[2].role == "uncertainty"
    )  # Unspecified confidence, not invented sigma.
    if index == 1:
        assert ds["pressure_gpa"].max() == 54.6
        assert ds["included_in_fit"][ds["sample"] == "B12_30"] == 0
    else:
        assert min(ds["temperature_k"]) == 1018
        assert max(ds["temperature_k"]) == 2553
        assert max(ds["pressure_gpa"]) == 62


@pytest.mark.parametrize("material", PUBLISHED)
def test_identity_source_parameters_and_interchange(material):
    doc = get_material_document(material)
    Draft202012Validator(eosmat_schema()).validate(doc)
    contents = Counter()
    for site in doc["atom_sites"]:
        contents[site["element"]] += int(site["wyckoff"][:-1]) * site["occupancy"]
    assert contents == {"Th": 4, "O": 8}
    assert doc["formula_units_per_cell"] == 4
    lattice = doc["lattice"]
    cell = lattice["a"] * lattice["b"] * lattice["c"]
    assert cell == pytest.approx(
        175.22 if material == "thorianite" else 145.7, abs=0.05
    )
    record = doc["eos_records"][0]
    v0, k0, c = PUBLISHED[material]
    assert record["eos"]["parameters"] == pytest.approx(
        dict(V0=v0 * CELL_PER_MOLAR, K0=k0, K0_prime=4)
    )
    assert record["parameter_errors"]["V0"] == pytest.approx(
        ERRORS[material][0] * CELL_PER_MOLAR
    )
    assert record["thermal"]["parameters"] == dict(Tr=300, alpha_KT=c)
    assert record["parameter_errors"]["K0_prime"] is None
    assert record["fixed_parameters"] == (
        ["V0", "K0_prime"] if material == "thorianite" else ["K0_prime"]
    )
    assert (
        record["pressure_calibration"]["recalculation"]["status"]
        == "missing_calibrant_observations"
    )
    obj = Material.from_eosmat(doc)
    exported = obj.to_eosmat()
    for key in ("datasets", "source", "atom_sites", "lattice"):
        assert exported[key] == doc[key]
    for key in (
        "scientific_validation",
        "parameter_provenance",
        "pressure_calibration",
        "fit_datasets",
    ):
        assert exported["eos_records"][0][key] == record[key]
    restored = Material.from_eosmat(obj.to_eosmat())
    assert restored.eos_records[0].pressure(150, 1500) == pytest.approx(
        obj.eos_records[0].pressure(150, 1500)
    )


@pytest.mark.parametrize("material", PUBLISHED)
def test_independent_thermal_equation_roundtrips_and_validity(material):
    record = get_material(material).eos_records[0]
    pub = PUBLISHED[material]
    volumes = np.array([pub[0], pub[0] * 0.95, pub[0] * 0.85])[:, None]
    temperatures = np.array([300, 1500, 2500])[None, :]
    expected = source_pressure(volumes, temperatures, pub)
    assert record.pressure(volumes * CELL_PER_MOLAR, temperatures) == pytest.approx(
        expected, abs=1e-10
    )
    assert record.eos.thermal_pressure(volumes * CELL_PER_MOLAR, 300) == pytest.approx(
        np.zeros((3, 1)), abs=1e-12
    )
    assert record.pressure(pub[0] * CELL_PER_MOLAR, 300) == pytest.approx(0, abs=1e-12)
    assert record.eos.bulk_modulus(pub[0] * CELL_PER_MOLAR, 300) == pytest.approx(
        pub[1]
    )
    pressures = [18, 27] if material == "thorianite" else [20, 60]
    temps = [1200, 1900]
    v = record.volume(pressures, temps, check_validity=True)
    assert record.pressure(v, temps, check_validity=True) == pytest.approx(
        pressures, abs=1e-8
    )
    with pytest.raises(ValueError):
        record.volume(80, 1500, check_validity=True)
    if material == "tho2_cotunnite":
        with pytest.raises(ValueError):
            record.volume(20, 300, check_validity=True)


def test_source_benchmarks_and_all_weighting_diagnostics():
    results = reproduce()
    for material, result in results.items():
        assert result["native_max_difference_gpa"] < 1e-10
        assert all(b["within_tolerance"] for b in result["benchmarks"])
        assert all(f["solver_success"] for f in result["refits"].values())
        assert result["published_rmse_gpa"] < 1.0
        record = get_material_document(material)["eos_records"][0]
        outcome = ledger_outcome(record)
        assert outcome["status"] == "similar"
        assert all(p["within_combined_2sigma"] is None for p in outcome["parameters"])
    assert results["thorianite"]["refits"]["unweighted"]["parameters_molar"][
        1
    ] == pytest.approx(200.09516, abs=1e-4)
    assert results["tho2_cotunnite"]["refits"]["effective_errors"][
        "parameters_molar"
    ] == pytest.approx([24.763289, 191.039762, 0.003524437], rel=1e-6)
