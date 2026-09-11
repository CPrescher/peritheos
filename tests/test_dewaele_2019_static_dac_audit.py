from __future__ import annotations

import csv
import importlib.util
import json
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_dewaele_2019_static_dac.py"
FE_DATA = (
    ROOT / "peritheos" / "data" / "datasets" / "iron-dewaele-2006-epaps-compression.csv"
)
REFIT_DATA = ROOT / "docs" / "data" / "dewaele-2019-static-dac-refit.json"


def _load_script():
    spec = importlib.util.spec_from_file_location("dewaele_2019_static_dac", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_dewaele_2019_audit_preserves_partial_reproduction_status(assert_audit_close):
    result = _load_script().reproduce()

    stored = json.loads(REFIT_DATA.read_text(encoding="utf-8"))
    # Subtraction amplifies relative noise in near-zero parameter differences.
    # Check those against each report's coefficients, then compare fitted
    # parameters and all remaining metrics at the usual solver precision.
    for audit in (stored, result):
        for group in ("row_level_refits", "source_publication_checks"):
            for record in audit[group].values():
                for fit in record.values():
                    if (
                        not isinstance(fit, dict)
                        or "absolute_parameter_difference" not in fit
                    ):
                        continue
                    differences = fit.pop("absolute_parameter_difference")
                    assert differences == pytest.approx(
                        [
                            abs(fitted - published)
                            for fitted, published in zip(
                                fit["parameters"], record["published_atomic_parameters"]
                            )
                        ],
                        rel=1e-12,
                        abs=1e-12,
                    )
    assert_audit_close(result, stored)

    assert result["catalog_record_count"] == 30
    assert len(result["row_level_refits"]) == 28
    assert result["complete_family_reproducible"] is False
    assert set(result["source_gaps"]) == {"beryllium_hcp", "lead_hcp"}

    for refit in result["row_level_refits"].values():
        assert refit["unweighted_pressure_residual_fit"]["within_reported_95pct_errors"]

    gold = result["row_level_refits"]["gold_dewaele_2019_mao_vinet"]
    assert gold["rows"] == 86
    assert gold["unweighted_pressure_residual_fit"]["parameters"] == pytest.approx(
        [16.98316389, 166.36665274, 5.46613027], rel=2e-8
    )

    lead = result["row_level_refits"]["lead_hcp_dewaele_2019_dor_vinet"]
    assert lead["rows"] == 17
    assert lead["scope"].startswith("current_run_only")

    iron_mao = result["row_level_refits"]["iron_dewaele_2019_mao_vinet"]
    assert iron_mao["rows"] == 53
    assert iron_mao["unweighted_pressure_residual_fit"]["parameters"] == pytest.approx(
        [11.20800948, 164.63787093, 4.96211218], rel=2e-8
    )

    iron_dor = result["row_level_refits"]["iron_dewaele_2019_dor_vinet"]
    assert iron_dor["rows"] == 53
    assert iron_dor["unweighted_pressure_residual_fit"]["parameters"] == pytest.approx(
        [11.17599418, 168.64954660, 5.32808152], rel=2e-8
    )

    iron_2006 = result["source_publication_checks"]["iron_dewaele_2006_helium_vinet"]
    assert iron_2006["rows"] == 37
    assert iron_2006["unweighted_pressure_residual_fit"]["within_reported_95pct_errors"]


def test_complete_fe_epaps_table_preserves_runs_phases_and_raw_gauges():
    with FE_DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 63
    assert Counter(row["phase"] for row in rows) == {"hcp": 53, "bcc": 10}
    assert Counter(int(row["run"]) for row in rows) == {
        1: 18,
        2: 6,
        3: 11,
        4: 12,
        5: 16,
    }
    assert all(row["ruby_wavelength_angstrom"] for row in rows if int(row["run"]) <= 3)
    assert all(
        row["tungsten_lattice_parameter_angstrom"]
        for row in rows
        if int(row["run"]) >= 4
    )
