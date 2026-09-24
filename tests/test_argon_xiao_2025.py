"""Primary workbook and source-table checks for the bundled Xiao record."""

import numpy as np
import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_argon_xiao_2025 import reproduce


def test_xiao_workbook_and_neutron_reproduction():
    report = reproduce()
    assert (
        max(
            abs(x)
            for x in report["workbook_checkpoint"]["relative_differences"].values()
        )
        < 1e-12
    )
    assert report["neutron_table5"]["rows"] == 22
    assert report["neutron_table5"]["relative_volume_rms_percent"] == pytest.approx(
        0.09517622713
    )
    assert report["global_refit"]["status"] == "not_reproduced"


def test_xiao_bundled_cell_molar_conversion():
    material = Material.from_eosmat(get_material_document("argon_fcc"))
    record = material.get_eos_record("argon_fcc_xiao_2025_helmholtz")
    assert record.pressure(23 * 4 / 0.602214076, 70) == pytest.approx(
        0.0783517375843629, abs=1e-10
    )
    for t, p, v in [
        (70, 0.01, 23.827471251113984),
        (83.806, 0.000068891, 24.604805616824237),
    ]:
        assert record.volume(p, t) == pytest.approx(v * 4 / 0.602214076, rel=1e-10)
    assert np.isfinite(record.pressure(15 * 4 / 0.602214076, 760))
