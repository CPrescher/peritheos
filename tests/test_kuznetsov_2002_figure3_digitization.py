from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_kuznetsov_2002_figure3_digitization.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("kuznetsov_figure3", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_kuznetsov_figure3_hcp_digitization_is_reproducibly_calibrated():
    result = _load_script().check()

    assert result["rows"] == 21
    assert result["open_symbols"] == 14
    assert result["solid_symbols"] == 7
    assert result["pressure_range_gpa"] == pytest.approx([10.527, 35.121])
    assert result["maximum_stored_rounding_difference"] < 5e-4
