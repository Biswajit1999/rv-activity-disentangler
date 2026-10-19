"""Night-clustered uncertainty diagnostics for irregularly sampled solar RV data.

All times are BJDs in days, velocities and uncertainties are in m/s. Results
are descriptive: no planet, solar rotation, or instrumental cause is inferred.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from scipy.optimize import minimize_scalar


def _validated(data: Mapping[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return time, RV, sigma and floor(BJD) night identifiers."""
    try:
        t = np.asarray(data["bjd"], dtype=float)
        v = np.asarray(data["rvRelativeMetersPerSecond"], dtype=float)
        sigma = np.asarray(data["uncertaintyMetersPerSecond"], dtype=float)
    except KeyError as exc:
        raise ValueError(f"missing required field: {exc.args[0]}") from exc
    if any(a.ndim != 1 for a in (t, v, sigma)) or len(t) != len(v) or len(t) != len(sigma):
        raise ValueError("time, RV and uncertainty must be matching one-dimensional arrays")
    if len(t) < 4 or not all(np.all(np.isfinite(a)) for a in (t, v, sigma)):
        raise ValueError("at least four finite observations are required")
    if np.any(sigma <= 0):
        raise ValueError("formal RV uncertainties must be strictly positive")
    nights = np.floor(t).astype(np.int64)
    if np.unique(nights).size < 2:
        raise ValueError("at least two observed BJD-day groups are required")
    return t, v, sigma, nights


def _night_means(v: np.ndarray, sigma2: np.ndarray, nights: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    groups = np.unique(nights)
    means = np.empty(groups.size)
    variances = np.empty(groups.size)
    for i, day in enumerate(groups):
        use = nights == day
        w = 1.0 / sigma2[use]
        means[i] = np.dot(w, v[use]) / w.sum()
        variances[i] = 1.0 / w.sum()
    return groups, means, variances


def profile_intrinsic_jitter(data: Mapping[str, np.ndarray]) -> dict[str, float | int | bool]:
    """Profile a shared additive white-noise variance after fitting night offsets.

    Gaussian maximum likelihood includes log-variance normalization and refits
    each nightly offset at every candidate jitter. It is *not* an estimator of
    uniquely stellar jitter: correlated activity and instrumental terms remain.
    """
    _, v, sigma, nights = _validated(data)
    unique = np.unique(nights)
    if len(v) <= len(unique):
        raise ValueError("nightly offsets leave no within-night degrees of freedom")

    def nll(jitter2: float) -> float:
        variance = sigma**2 + jitter2
        _, means, _ = _night_means(v, variance, nights)
        model = means[np.searchsorted(unique, nights)]
        return 0.5 * float(np.sum(np.log(2 * np.pi * variance) + (v - model)**2 / variance))

    # Broad deterministic bracket, with explicit comparison to the zero boundary.
    span = float(np.ptp(v))
    upper = max(span**2 * 4, float(np.max(sigma**2)) * 100, 1e-6)
    result = minimize_scalar(nll, bounds=(0.0, upper), method="bounded", options={"xatol": upper * 1e-12})
    if not result.success:
        raise RuntimeError("jitter profile optimization did not converge")
    best2 = float(result.x) if result.fun < nll(0.0) else 0.0
    variance = sigma**2 + best2
    _, means, _ = _night_means(v, variance, nights)
    model = means[np.searchsorted(unique, nights)]
    chi2 = float(np.sum((v - model)**2 / variance))
    return {
        "observations": len(v),
        "nights": len(unique),
        "intrinsic_white_jitter_mps": float(np.sqrt(best2)),
        "at_zero_boundary": bool(best2 == 0.0),
        "profile_delta_log_likelihood_vs_zero": float(nll(0.0) - nll(best2)),
        "within_night_reduced_chi2": chi2 / (len(v) - len(unique)),
    }


def night_cluster_bootstrap(
    data: Mapping[str, np.ndarray], *, replicates: int = 2000, seed: int = 20261019
) -> dict[str, float | int | list[float]]:
    """Bootstrap the equally night-weighted mean by resampling whole nights.

    Night-level means, not individual exposures, are resampled. The interval
    describes sensitivity to the observed nights, not the full solar activity
    cycle, and is not an instrumental zero-point calibration.
    """
    if not isinstance(replicates, int) or replicates < 100:
        raise ValueError("replicates must be an integer >= 100")
    _, v, sigma, nights = _validated(data)
    groups, means, _ = _night_means(v, sigma**2, nights)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(groups), size=(replicates, len(groups)))
    draws = means[indices].mean(axis=1)
    interval = np.percentile(draws, [2.5, 97.5])
    return {
        "observed_nights": len(groups),
        "replicates": replicates,
        "night_equal_weight_mean_mps": float(means.mean()),
        "bootstrap_standard_error_mps": float(np.std(draws, ddof=1)),
        "bootstrap_95pct_percentile_interval_mps": [float(x) for x in interval],
        "naive_exposure_weighted_formal_error_mps": float(1.0 / np.sqrt(np.sum(1 / sigma**2))),
    }


def leave_one_night_out(data: Mapping[str, np.ndarray]) -> list[dict[str, float | int]]:
    """Assess sensitivity of between-night scatter to each observed BJD day."""
    _, v, sigma, nights = _validated(data)
    groups, means, _ = _night_means(v, sigma**2, nights)
    if len(groups) < 3:
        raise ValueError("at least three nights are required for leave-one-night-out")
    baseline = float(np.std(means, ddof=1))
    rows: list[dict[str, float | int]] = []
    for i, day in enumerate(groups):
        retained = np.delete(means, i)
        scatter = float(np.std(retained, ddof=1))
        rows.append({
            "omitted_bjd_day": int(day),
            "remaining_nights": len(retained),
            "between_night_scatter_mps": scatter,
            "delta_scatter_from_full_mps": scatter - baseline,
            "remaining_night_equal_mean_mps": float(retained.mean()),
        })
    return rows
