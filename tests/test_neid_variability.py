import numpy as np
import pytest

from research.neid_variability import nightly_variability, structure_function


def _fixture():
    return {
        "bjd": np.array([1.1, 1.2, 2.1, 2.2]),
        "rvRelativeMetersPerSecond": np.array([-1.0, 1.0, 2.0, 4.0]),
        "uncertaintyMetersPerSecond": np.ones(4),
        "snr": np.full(4, 100.0),
        "sunAltitudeDegrees": np.full(4, 50.0),
    }


def test_nightly_weighted_statistics() -> None:
    rows = nightly_variability(_fixture())
    assert [row["weighted_mean_rv_mps"] for row in rows] == pytest.approx([0.0, 3.0])
    assert all(row["formal_mean_error_mps"] == pytest.approx(1 / np.sqrt(2)) for row in rows)


def test_structure_function_uses_unique_pairs() -> None:
    rows = structure_function(_fixture(), np.array([0.05, 0.15, 2.0]))
    assert sum(row["pair_count"] for row in rows) == 6
    assert all(row["median_abs_difference_over_sqrt2_mps"] >= 0 for row in rows)


def test_structure_function_rejects_unsorted_edges() -> None:
    with pytest.raises(ValueError, match="positive and increasing"):
        structure_function(_fixture(), np.array([0.1, 0.05]))

