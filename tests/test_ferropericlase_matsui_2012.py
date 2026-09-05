"""Primary-source regressions for Matsui et al. (2012) ferropericlase."""

import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import Avogadro, R
from scipy.integrate import quad
from scipy.optimize import least_squares

from peritheos.eosmat import validate_eosmat_document
from peritheos.materials import Material

ROOT = Path(__file__).resolve().parents[1]
DOI = "10.2138/am.2012.3937"

CASES = (
    {
        "material": "ferropericlase_mg83fe17",
        "record": "ferropericlase_mg83fe17_matsui_2012_bm3_mgd_1",
        "dataset": "ferropericlase_mg83fe17_matsui_2012_table1_pvt",
        "csv": "ferropericlase-mg83fe17-matsui-2012-table1-pvt.csv",
        "sha256": "b933ef0489ea28a42917006da0de84ef92142c721418ce1fa218c31e5d8849dd",
        "rows": 34,
        "fit_rows": 23,
        "composition": {"Mg": 3.32, "Fe": 0.68, "O": 4.0},
        "lattice": 4.2331,
        "parameters": (75.849, 4.08, 1.53, 0.7),
        "errors": (0.011, 0.02, 0.04, 0.2),
        "refit": (75.80704803, 4.13453651, 1.52863286, 0.68229851),
        "max_printed_delta": 0.014,
        "rmse": 0.123302748733,
        "last_row": (1100.0, 61.542, 55.47, 56.63),
    },
    {
        "material": "ferropericlase_mg75fe25",
        "record": "ferropericlase_mg75fe25_matsui_2012_bm3_mgd_1",
        "dataset": "ferropericlase_mg75fe25_matsui_2012_table2_pvt",
        "csv": "ferropericlase-mg75fe25-matsui-2012-table2-pvt.csv",
        "sha256": "22c7a9930b69528f635593e376c8124daf1e6e004a8727acc9cc10c768384cf8",
        "rows": 39,
        "fit_rows": 30,
        "composition": {"Mg": 3.0, "Fe": 1.0, "O": 4.0},
        "lattice": 4.2427,
        "parameters": (76.372, 4.22, 1.64, 0.7),
        "errors": (0.016, 0.03, 0.04, 0.2),
        "refit": (76.37279132, 4.22184341, 1.64549639, 0.74501033),
        "max_printed_delta": 0.006,
        "rmse": 0.142609785218,
        "last_row": (1100.0, 62.110, 55.52, 56.92),
    },
)


def _load(case):
    path = ROOT / "peritheos" / "data" / "materials" / f"{case['material']}.eosmat"
    document = json.loads(path.read_text(encoding="utf-8"))
    csv_path = ROOT / "peritheos" / "data" / "datasets" / case["csv"]
    rows = list(csv.DictReader(io.StringIO(csv_path.read_text(encoding="utf-8"))))
    return document, csv_path, rows


def _bm3_pressure(volume, v0, k0_prime):
    ratio = v0 / volume
    return (
        1.5
        * 160.0
        * (ratio ** (7.0 / 3.0) - ratio ** (5.0 / 3.0))
        * (1.0 + 0.75 * (k0_prime - 4.0) * (ratio ** (2.0 / 3.0) - 1.0))
    )


def _debye_function_3(argument):
    integral = quad(
        lambda value: value**3 / np.expm1(value),
        0.0,
        argument,
        epsabs=1.0e-12,
        epsrel=1.0e-12,
    )[0]
    return 3.0 * integral / argument**3


def _source_pressure(volume, temperature, v0, k0_prime, gamma0, q):
    gamma = gamma0 * (volume / v0) ** q
    theta = 500.0 * np.exp((gamma0 - gamma) / q)

    def energy(at_temperature):
        return (
            3.0 * 2.0 * R * at_temperature * _debye_function_3(theta / at_temperature)
        )

    molar_volume_j_per_bar = volume * Avogadro * 1.0e-25 / 4.0
    thermal_pressure = (
        gamma * (energy(temperature) - energy(300.0)) / molar_volume_j_per_bar / 1.0e4
    )
    return _bm3_pressure(volume, v0, k0_prime) + thermal_pressure


def _staged_unweighted_refit(rows, published):
    fit_rows = [row for row in rows if row["used_in_published_fit"] == "1"]
    room_temperature = [row for row in fit_rows if row["temperature_k"] == "300"]
    heated = [row for row in fit_rows if row["temperature_k"] != "300"]

    cold = least_squares(
        lambda parameters: np.array(
            [
                _bm3_pressure(float(row["volume_a3"]), parameters[0], parameters[1])
                - float(row["pressure_gpa"])
                for row in room_temperature
            ]
        ),
        published[:2],
    ).x
    thermal = least_squares(
        lambda parameters: np.array(
            [
                _source_pressure(
                    float(row["volume_a3"]),
                    float(row["temperature_k"]),
                    published[0],
                    published[1],
                    parameters[0],
                    parameters[1],
                )
                - float(row["pressure_gpa"])
                for row in heated
            ]
        ),
        published[2:],
    ).x
    return np.concatenate((cold, thermal))


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["material"])
def test_matsui_2012_records_are_executable_and_source_reproduced(case):
    document, _, rows = _load(case)
    validate_eosmat_document(document)
    assert document["formula_units_per_cell"] == 4
    assert document["space_group"] == "Fm-3m"
    assert document["space_group_number"] == 225
    assert document["lattice"]["a"] == case["lattice"]

    atom_counts = {"Mg": 0.0, "Fe": 0.0, "O": 0.0}
    for site in document["atom_sites"]:
        multiplicity = int(site["wyckoff"][:-1])
        atom_counts[site["element"]] += multiplicity * site["occupancy"]
    assert atom_counts == pytest.approx(case["composition"])

    source_record = document["eos_records"][0]
    assert source_record["identifier"] == case["record"]
    assert source_record["reference"]["doi"] == DOI
    assert source_record["eos"]["parameters"] == {
        "V0": case["parameters"][0],
        "K0": 160.0,
        "K0_prime": case["parameters"][1],
    }
    assert source_record["fixed_parameters"] == ["K0"]
    assert source_record["thermal"]["debye_temperature_law"] == ("integrated_gruneisen")
    assert source_record["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 500.0,
        "gamma0": case["parameters"][2],
        "q": case["parameters"][3],
        "n": 2.0,
    }
    assert source_record["thermal"]["fixed_parameters"] == ["Tr", "theta0", "n"]

    record = Material.from_eosmat(document).get_eos_record(case["record"])
    assert record.pressure(record.reference_volume, 300.0) == pytest.approx(
        0.0, abs=1.0e-12
    )
    assert record.thermal_pressure_increment(0.85 * record.reference_volume, 300.0) == (
        pytest.approx(0.0, abs=1.0e-12)
    )

    fit_rows = [row for row in rows if row["used_in_published_fit"] == "1"]
    calculated = np.array(
        [
            record.pressure(
                float(row["volume_a3"]),
                float(row["temperature_k"]),
                check_validity=False,
            )
            for row in fit_rows
        ]
    )
    printed = np.array([float(row["calculated_pressure_gpa"]) for row in fit_rows])
    observed = np.array([float(row["pressure_gpa"]) for row in fit_rows])
    assert len(fit_rows) == case["fit_rows"]
    assert np.max(np.abs(calculated - printed)) < case["max_printed_delta"]
    assert np.sqrt(np.mean((calculated - observed) ** 2)) == pytest.approx(
        case["rmse"], abs=1.0e-10
    )

    high_temperature_row = next(
        row for row in reversed(fit_rows) if row["temperature_k"] == "1100"
    )
    volume = float(high_temperature_row["volume_a3"])
    pressure = record.pressure(volume, 1100.0, check_validity=False)
    assert record.volume(pressure, 1100.0, check_validity=False) == pytest.approx(
        volume
    )
    assert record.within_validity(volume, 1100.0)
    assert not record.within_validity(volume, 1200.0)


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["material"])
def test_matsui_2012_primary_tables_are_transcribed_verbatim(case):
    document, csv_path, rows = _load(case)
    dataset = document["datasets"][0]
    assert dataset["identifier"] == case["dataset"]
    assert dataset["used_by_eos_records"] == [case["record"]]
    assert hashlib.sha256(csv_path.read_bytes()).hexdigest() == case["sha256"]
    assert len(rows) == case["rows"]
    assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
    assert rows[0]["temperature_k"] == "300"
    assert rows[0]["gold_volume_ratio"] == "1.0"
    assert rows[0]["pressure_gpa"] == "0.0"
    assert rows[0]["pressure_uncertainty_gpa"] == ""
    last_temperature, last_volume, last_pressure, last_calculated = case["last_row"]
    assert float(rows[-1]["temperature_k"]) == last_temperature
    assert float(rows[-1]["volume_a3"]) == last_volume
    assert float(rows[-1]["pressure_gpa"]) == last_pressure
    assert float(rows[-1]["calculated_pressure_gpa"]) == last_calculated
    assert sum(row["used_in_published_fit"] == "1" for row in rows) == case["fit_rows"]


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["material"])
def test_matsui_2012_published_staged_fit_is_independently_reproduced(case):
    _, _, rows = _load(case)
    refit = _staged_unweighted_refit(rows, case["parameters"])
    assert refit == pytest.approx(case["refit"], abs=5.0e-8)

    published = np.array(case["parameters"])
    errors = np.array(case["errors"])
    if case["material"] == "ferropericlase_mg83fe17":
        assert np.all(np.abs(refit[2:] - published[2:]) <= errors[2:])
        assert np.any(np.abs(refit[:2] - published[:2]) > errors[:2])
    else:
        assert np.all(np.abs(refit - published) <= errors)
