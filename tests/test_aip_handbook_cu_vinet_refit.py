import pytest

from scripts.reproduce_aip_handbook_cu_vinet_refit import (
    load_observations,
    refit,
)


def test_checked_cu_transcription_includes_corrected_printing_error():
    pressure, relative_volume = load_observations()
    assert len(pressure) == 25
    assert pressure[0] == pytest.approx(1.5)
    assert relative_volume[0] == pytest.approx(0.9900)
    assert pressure[-1] == pytest.approx(34.0)
    assert relative_volume[-1] == pytest.approx(0.849)


def test_cu_vinet_l8_refit_is_reproducible_and_improves_the_objective():
    result = refit()
    fitted = result["refit"]
    published = result["published"]

    assert result["observations"] == 25
    assert fitted["k0_gpa"] == pytest.approx(139.23431, abs=2.0e-5)
    assert fitted["k0_prime"] == pytest.approx(5.001127, abs=2.0e-6)
    assert fitted["rmse_gpa"] == pytest.approx(0.07015421, abs=1.0e-8)
    assert fitted["objective8_gpa8"] < published["objective8_gpa8"]
