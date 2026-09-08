import csv
import hashlib
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
IDENTIFIER = "ca_perovskite_tetragonal_chen_2018_vinet"
DATASET_IDENTIFIER = "ca_perovskite_tetragonal_chen_2018_table1_compression"


def test_chen_2018_record_preserves_fit_and_executes():
    document = get_material_document("ca_perovskite_tetragonal")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == IDENTIFIER
    )
    assert record["reference"]["doi"] == "10.2138/am-2018-6087"
    assert record["eos"]["type"] == "Vinet"
    assert record["eos"]["model"] == "vinet"
    assert record["eos"]["parameters"] == {"V0": 185.2, "K0": 223.0, "K0_prime": 4.0}
    assert record["parameter_errors"] == {"V0": 0.4, "K0": 6.0, "K0_prime": None}
    assert record["parameter_error_confidence"] == pytest.approx(0.9544997361036416)
    assert record["fixed_parameters"] == ["K0_prime"]
    assert record["fit_datasets"] == [DATASET_IDENTIFIER]
    assert record["experimental_pressure_range_gpa"] == [28.824, 62.477]
    assert record["pressure_calibration"]["methods"][0]["reference"]["doi"] == (
        "10.1002/2016JB013811"
    )
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    assert eos.pressure(185.2) == pytest.approx(0.0, abs=1e-12)
    assert eos.pressure(0.9 * 185.2) == pytest.approx(28.9286907881)
    assert eos.volume(eos.pressure(0.85 * 185.2)) == pytest.approx(
        0.85 * 185.2, rel=1e-10
    )


def test_chen_2018_table1_transcription_is_complete_and_checksum_tracked():
    document = get_material_document("ca_perovskite_tetragonal")
    dataset = next(
        row for row in document["datasets"] if row["identifier"] == DATASET_IDENTIFIER
    )
    path = ROOT / "peritheos/data" / dataset["resource"]["path"]
    assert (
        hashlib.sha256(path.read_bytes()).hexdigest() == dataset["resource"]["sha256"]
    )

    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert [row["run_id"] for row in rows] == [
        "71013",
        "81021",
        "81030",
        "71050",
        "81066",
        "81074",
        "71088",
    ]
    assert [float(row["pressure_pt_ye_gpa"]) for row in rows] == [
        28.824,
        29.138,
        32.307,
        41.972,
        47.878,
        53.884,
        62.477,
    ]
    assert {row["used_in_published_fit"] for row in rows} == {"1"}
    for row in rows:
        i4_volume = float(row["i4mcm_lattice_a_angstrom"]) ** 2 * float(
            row["i4mcm_lattice_c_angstrom"]
        )
        p4_volume = float(row["p4mmm_lattice_a_angstrom"]) ** 2 * float(
            row["p4mmm_lattice_c_angstrom"]
        )
        assert float(row["i4mcm_volume_a3_conventional_cell"]) == pytest.approx(
            i4_volume, abs=5e-13
        )
        assert float(row["p4mmm_volume_a3_conventional_cell"]) == pytest.approx(
            p4_volume, abs=5e-13
        )


def test_chen_2018_reproduction_script():
    path = ROOT / "scripts/reproduce_chen_2018_casio3.py"
    spec = importlib.util.spec_from_file_location("reproduce_chen_2018", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    protocol = result["source_protocol"]
    assert protocol["equation"] == "Vinet"
    assert protocol["observations"] == 7
    assert protocol["excluded_runs"] == []
    assert protocol["fixed_parameters"] == {"K0_prime": 4.0}
    assert protocol["reported_weights"] is None
    assert protocol["reported_pressure_uncertainties"] is None

    published = result["published_curve"]
    assert published["V0_a3_conventional_cell"] == pytest.approx(185.2)
    assert published["K0_gpa"] == pytest.approx(223.0)
    assert published["pressure_rmse_gpa"] == pytest.approx(0.688528394433)
    assert published["pressure_max_abs_residual_gpa"] == pytest.approx(1.62701596404)

    refits = result["source_scope_refits"]
    pressure_fit = refits["unweighted_pressure_residual"]
    assert [
        pressure_fit["V0_a3_conventional_cell"],
        pressure_fit["K0_gpa"],
    ] == pytest.approx([185.978178086, 215.290968798])
    weighted = refits["propagated_volume_1sigma_weighted"]
    assert [weighted["V0_a3_conventional_cell"], weighted["K0_gpa"]] == pytest.approx(
        [185.624949627, 218.268891250]
    )
    assert np.max(np.abs(published["pressure_residuals_gpa"])) == pytest.approx(
        published["pressure_max_abs_residual_gpa"]
    )
