"""Published zero errors and state-only Monte Carlo uncertainty."""

import numpy as np
import pytest

from peritheos import get_eos_record
from peritheos.eos.rt import BM3
from peritheos.uncertainty import EOSUncertainty


@pytest.mark.parametrize("material", ["mg075fe025o", "mg0490fe0483ti0027o"])
def test_published_zero_errors_propagate_with_precision_qualification(material):
    record = get_eos_record(f"{material}_solomatova_2016_hs_bm3_reference_1")
    assert record.parameter_errors["V0"] == 0
    result = record.pressure_with_uncertainty(70)
    assert np.isfinite(result.standard_error) and result.standard_error > 0
    assert any("does not establish exactness" in s for s in result.assumptions)


def test_state_only_monte_carlo_matches_independent_local_volume_derivative():
    eos = BM3(10, 100, 4)
    uncertainty = EOSUncertainty.state_only(eos)
    result = uncertainty.pressure(
        9, volume_sigma=0.001, method="monte_carlo", sample_count=4000, random_state=1
    )
    assert result.standard_error == pytest.approx(
        eos.bulk_modulus(9) / 9 * 0.001, rel=0.04
    )
    assert not any("parameter uncertainty sampled" in s for s in result.assumptions)
    assert any("state-variable errors sampled" in s for s in result.assumptions)
    deterministic = uncertainty.pressure(
        9, method="monte_carlo", sample_count=10, random_state=1
    )
    assert deterministic.standard_error == 0
