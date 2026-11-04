"""Scientific checks of the staged Brown calibration reconstruction."""

import numpy as np
import pytest

from scripts.ferrous_ringwoodite_brown_replay import (
    TABLE1,
    TABLE2,
    Brown1999,
    replay,
    source_pdf_check,
    specific_volume,
    vibrational_energy,
)
from scripts.ferrous_ringwoodite_recover import read_csv
from scripts.ferrous_ringwoodite_reproduce import STAGE, R


@pytest.fixture(scope="module")
def replayed():
    return replay()


def test_brown_primary_tables_cell_contents_and_thermal_normalization():
    assert len(read_csv(TABLE1)) == 21
    assert len(read_csv(TABLE2)) == 63
    # The preliminary secondary code transcription incorrectly had 10.07.
    assert read_csv(TABLE1)[4]["static_pressure_gpa"] == "10.70"
    assert specific_volume(5.2405) == pytest.approx(0.37076512762511543)
    # NaCl has two atoms per formula unit: high-T energy tends to 6RT.
    assert vibrational_energy(1e8, 375) / (6 * R * 1e8) == pytest.approx(1, abs=2e-6)
    model = Brown1999()
    assert model.table1(0.3709, 300) == pytest.approx(9.069532127311282)
    assert model.table2(0.3709, 300) == pytest.approx(9.19)
    assert model.table2(0.3709, 1200) == pytest.approx(11.60)
    # A real published-table disagreement must remain visible, not be calibrated away.
    assert model.table1(0.3709, 300) - model.table2(0.3709, 300) == pytest.approx(
        -0.120467872688718
    )


def test_primary_table_thermal_increment_check_and_explicit_domain(replayed):
    result = replayed["primary_table_consistency"]
    assert result["table2_pressure_cells_checked"] == 571
    assert 0.15 < result["maximum_abs_table1_minus_table2_gpa"] < 0.17
    assert result["maximum_abs_300k_anchored_thermal_increment_difference_gpa"] < 0.016
    model = Brown1999()
    with pytest.raises(ValueError, match="explicit extension"):
        model.table2(0.371, 1407.5)
    with pytest.raises(ValueError, match="Specific volume"):
        model.table1(0.30, 300)
    with pytest.raises(ValueError, match="Specific volume"):
        model.table2(0.50, 300, extend=True)
    # Extension must be continuous and equal the independently integrated Debye increment.
    assert model.table2(0.371, 1407.5, extend=True) - model.table2(
        0.371, 1200
    ) == pytest.approx(model.thermal(0.371, 1407.5) - model.thermal(0.371, 1200))


def test_every_marker_reading_preserved_and_hot_residuals_remain(replayed):
    rows = replayed["rows"]
    assert len(rows) == 18
    assert sum(r["included_by_table_s1"] for r in rows) == 14
    assert replayed["table2_temperature_extension_measurements"] == [
        "Fay025heat",
        "Fay028heat",
    ]
    by_name = {r["measurement"]: r for r in rows}
    hot = by_name["Fay025heat"]
    assert hot["source_pressure_gpa"] == 11.48
    assert hot["sample_temperature_k"] == 2515
    assert hot["midpoint_marker_temperature_k"] == 1407.5
    assert hot["brown_table1_pressure_gpa"] == pytest.approx(11.955303153332636)
    assert hot["brown_table2_pressure_gpa"] == pytest.approx(12.064494040048901)
    assert hot["table1_minus_ringwoodite_gpa"] < -2
    assert hot["table2_minus_ringwoodite_gpa"] < -1.9
    assert by_name["Fay028heat"]["table2_minus_ringwoodite_gpa"] < -1
    assert (
        by_name["Fay008heat"]["marker_temperature_matching_ringwoodite_table2_k"]
        is None
    )
    assert replayed["calibration_explains_hot_residuals"] is False
    assert replayed["author_reduction_reproduced"] is False
    for r in rows:
        assert r["table2_lattice_only_pressure_sd_gpa"] > 0
        assert 0.0024 < r["table2_dp_dt_gpa_k"] < 0.0029


def test_temperature_diagnostics_are_inversions_not_observations(replayed):
    model = Brown1999()
    for row in replayed["rows"]:
        volume = row["nacl_specific_volume_cm3_g"]
        np.testing.assert_allclose(
            model.table1(volume, row["marker_temperature_matching_source_table1_k"]),
            row["source_pressure_gpa"],
            atol=1e-10,
        )
        if row["included_by_table_s1"]:
            assert (
                row["half_sample_temperature_diagnostic_k"]
                == row["midpoint_marker_temperature_k"] - 150
            )
    stats = replayed["accepted_14_statistics_gpa"]
    assert stats["half_sample_temperature_table1_minus_source_gpa"]["rms"] < 0.11
    assert stats["table1_minus_source_gpa"]["rms"] > 0.45
    # The temperature convention is unresolved despite the smaller diagnostic RMS.
    assert replayed["author_reduction_reproduced"] is False


def test_optional_complete_brown_pdf_transcription():
    if not (STAGE / "sources/brown-1999.pdf").is_file():
        pytest.skip(
            "User-supplied copyrighted PDF is not a redistributable test fixture"
        )
    result = source_pdf_check(STAGE / "sources")
    assert result["all_source_cells_and_blanks_match"]
    assert result["table1_numeric_cells"] == 126
    assert result["table2_populated_pressure_cells"] == 571
