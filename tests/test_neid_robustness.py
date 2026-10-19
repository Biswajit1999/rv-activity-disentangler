"""Deterministic numerical checks, not observational claims."""
import numpy as np
import pytest
from research.neid_robustness import profile_intrinsic_jitter, night_cluster_bootstrap, leave_one_night_out


def fixture(nights=5, samples=20, jitter=0.0, seed=3):
    rng = np.random.default_rng(seed)
    day = np.repeat(np.arange(nights) + 2459700, samples)
    t = day + np.tile(np.linspace(0.1, 0.4, samples), nights)
    offsets = np.repeat(np.linspace(-2, 2, nights), samples)
    sigma = np.full(len(day), 0.2)
    rv = offsets + rng.normal(0, np.hypot(sigma, jitter))
    return {"bjd": t, "rvRelativeMetersPerSecond": rv, "uncertaintyMetersPerSecond": sigma}


def test_jitter_zero_boundary_for_exact_constant_nights():
    d = fixture()
    d["rvRelativeMetersPerSecond"] = np.floor(d["bjd"]) - 2459700.0
    result = profile_intrinsic_jitter(d)
    assert result["at_zero_boundary"]
    assert result["within_night_reduced_chi2"] == pytest.approx(0)


def test_jitter_recovers_injected_noise_scale():
    result = profile_intrinsic_jitter(fixture(nights=8, samples=100, jitter=1.2))
    assert 1.0 < result["intrinsic_white_jitter_mps"] < 1.4
    assert result["profile_delta_log_likelihood_vs_zero"] > 100


def test_jitter_handles_heteroscedastic_errors():
    d = fixture()
    d["uncertaintyMetersPerSecond"][::2] *= 3
    assert np.isfinite(profile_intrinsic_jitter(d)["intrinsic_white_jitter_mps"])


def test_bootstrap_seed_reproducible():
    d = fixture()
    assert night_cluster_bootstrap(d, seed=77) == night_cluster_bootstrap(d, seed=77)


def test_bootstrap_exact_nightly_constant():
    d = fixture(nights=3)
    d["rvRelativeMetersPerSecond"] = np.floor(d["bjd"]) - 2459700.0
    result = night_cluster_bootstrap(d)
    assert result["night_equal_weight_mean_mps"] == pytest.approx(1)
    assert result["bootstrap_standard_error_mps"] > result["naive_exposure_weighted_formal_error_mps"]


def test_leave_one_night_out_identifies_outlier_night():
    d = fixture(nights=4, samples=10)
    d["rvRelativeMetersPerSecond"] = np.zeros(40)
    d["rvRelativeMetersPerSecond"][-10:] = 10
    rows = leave_one_night_out(d)
    assert len(rows) == 4
    assert rows[-1]["between_night_scatter_mps"] == pytest.approx(0)
    assert rows[-1]["delta_scatter_from_full_mps"] < 0


def test_bad_input_rejected():
    d = fixture()
    d["uncertaintyMetersPerSecond"][0] = 0
    with pytest.raises(ValueError, match="strictly positive"):
        profile_intrinsic_jitter(d)


def test_leave_one_night_out_requires_three_nights():
    with pytest.raises(ValueError, match="at least three nights"):
        leave_one_night_out(fixture(nights=2))


def test_bootstrap_minimum_replicates():
    with pytest.raises(ValueError, match="replicates"):
        night_cluster_bootstrap(fixture(), replicates=20)
