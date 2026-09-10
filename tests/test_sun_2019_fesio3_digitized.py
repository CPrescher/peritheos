import csv
import hashlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos/data/datasets"


def _module():
    path = ROOT / "scripts/reproduce_sun_2019_fesio3_digitized.py"
    spec = importlib.util.spec_from_file_location("sun_2019_digitized", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_digitized_s3_s4_dataset_is_complete():
    path = DATA / "fesio3-liquid-sun-2019-figures-s3-s4-digitized.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 9
    assert [float(row["volume_ratio_to_vx"]) for row in rows] == pytest.approx(
        [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    )
    assert min(float(row["cv_j_mol_k"]) for row in rows) == pytest.approx(170.078159)
    assert max(float(row["gamma"]) for row in rows) == pytest.approx(1.112172)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == (
        "a59993ad87fa47f4ed2a9d5ddb75ce951e30e9b46805ec9de78e56663243cdc6"
    )


def test_digitized_thermal_refit_and_bm4_diagnostic():
    result = _module().reproduce()
    comparison = result["thermal_parameter_comparison"]
    assert comparison["Cv"]["Cv_prime"]["digitized_refit"] == pytest.approx(
        205.3991, abs=1.0e-3
    )
    assert comparison["Cv"]["Vc"]["digitized_refit"] == pytest.approx(
        29.0843, abs=1.0e-3
    )
    assert comparison["gamma"]["gamma_Vx"]["digitized_refit"] == pytest.approx(
        0.347132, abs=1.0e-5
    )
    free = result["bm4_free_sensitivity"]
    assert free["k0_double_prime_lower_audit_bound_hit"] is True
    assert free["jacobian_condition_number"] > 1.0e7
    stable = result["bm4_stabilized_check"]
    assert stable["parameters"]["K0"] == pytest.approx(2.24587, abs=1.0e-4)
    assert stable["parameters"]["K0_prime"] == pytest.approx(23.0422, abs=1.0e-3)
