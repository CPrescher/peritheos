import pytest


def test_audit_comparison_allows_only_small_float_drift(assert_audit_close):
    assert_audit_close({"fit": [158.00001], "rows": 49}, {"fit": [158.0], "rows": 49})
    with pytest.raises(AssertionError, match=r"report.fit\[0\]"):
        assert_audit_close({"fit": [159.0]}, {"fit": [158.0]})


@pytest.mark.parametrize(
    ("actual", "expected"),
    [
        ({"rows": 50}, {"rows": 49}),
        ({"status": "parity"}, {"status": "similar"}),
        ({"rows": 49.0}, {"rows": 49}),
        ({"fit": [1.0, 2.0]}, {"fit": [1.0]}),
        ({"extra": None}, {}),
    ],
)
def test_audit_comparison_keeps_metadata_and_structure_exact(
    assert_audit_close, actual, expected
):
    with pytest.raises(AssertionError):
        assert_audit_close(actual, expected)
