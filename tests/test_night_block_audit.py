import numpy as np
import pytest

from research.night_block_audit import (
    quality_cut_decomposition, variance_components, leave_one_night_out,
    night_cluster_bootstrap,
)


def fixture():
    return {
        "bjd": np.repeat([10.1, 11.1, 12.1, 13.1], 8) + np.tile(np.arange(8) * 0.002, 4),
        "rvRelativeMetersPerSecond": np.repeat([0.0, 1.0, -1.0, 8.0], 8)
        + np.tile(np.array([-1, 1, -.5, .5, -1, 1, -.5, .5]), 4),
        "uncertaintyMetersPerSecond": np.ones(32),
        "snr": np.repeat([100, 250, 350, 450], 8).astype(float),
    }


def test_weighted_variance_identity_and_zero_between():
    rv = np.array([0., 2., 1., 3.])
    e = np.array([1., 2., 1., 2.])
    d = variance_components(rv, e, np.array([0, 0, 1, 1]))
    assert d["total_rms_mps"] ** 2 == pytest.approx(
        d["within_night_rms_mps"] ** 2 + d["between_night_rms_mps"] ** 2)
    assert d["between_night_rms_mps"] == pytest.approx(0.5)


def test_quality_cut_tracks_composition_and_invalid_cut():
    rows = quality_cut_decomposition(fixture(), cuts=(10, 200, 400, 999))
    assert rows[0]["n"] == 32 and rows[0]["valid"]
    assert rows[1]["n"] == 24
    assert not rows[2]["valid"] and not rows[3]["valid"]
    assert rows[0]["total_rms_mps"] > rows[0]["within_night_rms_mps"]


def test_leave_one_out_detects_outlying_night():
    rows = leave_one_night_out(fixture())
    assert len(rows) == 4
    assert min(rows, key=lambda x: x["delta_total_rms_mps"])["omitted_bjd_day"] == 13


def test_night_bootstrap_seed_and_finite_interval():
    a = night_cluster_bootstrap(fixture(), n_resamples=200, seed=7)
    b = night_cluster_bootstrap(fixture(), n_resamples=200, seed=7)
    assert a == b
    assert a["bootstrap_p025_mps"] <= a["within_night_rms_mps"] <= a["bootstrap_p975_mps"]


def test_duplicate_nights_do_not_change_within_scatter():
    a = fixture()
    b = {k: np.concatenate([v, v]) for k, v in a.items()}
    assert quality_cut_decomposition(a)[0]["within_night_rms_mps"] == pytest.approx(
        quality_cut_decomposition(b)[0]["within_night_rms_mps"])


def test_invalid_inputs_and_insufficient_nights():
    a = fixture()
    a["uncertaintyMetersPerSecond"][0] = 0
    with pytest.raises(ValueError, match="uncertainties"):
        quality_cut_decomposition(a)
    with pytest.raises(ValueError, match="three nights"):
        leave_one_night_out(fixture(), snr_cut=300)
    with pytest.raises(ValueError, match="resamples"):
        night_cluster_bootstrap(fixture(), n_resamples=3)
