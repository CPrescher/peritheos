"""The separately named refit must be usable without Brown's calculated tables."""

import json

import numpy as np
import pytest

from scripts import brown_nacl_1999_independent_eos as candidate
from scripts import brown_nacl_1999_reconstruct as base


@pytest.fixture(scope="module")
def artifact(tmp_path_factory):
    original = candidate.DEST
    candidate.DEST = tmp_path_factory.mktemp("nacl-refit") / "eos.json"
    original_read = base.read

    def guarded(name):
        assert not name.startswith("brown-table")
        return original_read(name)

    base.read = guarded
    try:
        record = candidate.build()
        return record, candidate.DEST
    finally:
        base.read = original_read
        candidate.DEST = original


def test_saved_eos_needs_no_source_data_or_refit(artifact, monkeypatch):
    record, path = artifact

    def forbidden(*args, **kwargs):
        raise AssertionError("Loading/evaluation must not refit or read source CSVs")

    monkeypatch.setattr(base, "read", forbidden)
    monkeypatch.setattr(base.Model, "__init__", forbidden)
    saved = candidate.load_eos(path)
    assert record["record_kind"] == "refit"
    assert saved.pressure(base.VREF, base.TREF) == pytest.approx(base.PREF, abs=1e-12)
    for T in [300, 600, 1000]:
        for P in [base.PREF, 3, 10, 20]:
            v = saved.volume(P, T)
            assert float(saved.pressure(v, T)) == pytest.approx(P, abs=1e-8)


def test_saved_spline_matches_refit_and_has_positive_compressibility(artifact):
    record, path = artifact
    saved = candidate.load_eos(path)
    refit = base.Model(base.Settings(**record["settings"]))
    v = np.linspace(saved.v.min() + 1e-5, saved.v.max() - 1e-5, 501)
    for T in [293, 600, 1023.15]:
        np.testing.assert_allclose(
            saved.pressure(v, T), refit.pressure(v, T), atol=1e-12
        )
        assert np.all(saved.kt(v, T) > 0)
    assert record["validation"]["cold_kprime_range"][0] > 0
    assert record["fit_provenance"]["covariance"] is None


def test_saved_domain_rejects_extrapolation_and_mismatched_constants(
    artifact, tmp_path
):
    record, path = artifact
    saved = candidate.load_eos(path)
    with pytest.raises(ValueError, match="Volume"):
        saved.pressure(saved.v.min() - 0.001, 300)
    with pytest.raises(ValueError, match="Temperature"):
        saved.pressure(0.4, 1200)
    record = json.loads(json.dumps(record))
    record["constants"]["N"] = 1
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="constants"):
        candidate.load_eos(bad)
