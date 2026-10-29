"""Independent night-block sensitivity checks for the retained NEID solar RV sample.

The unit of independent resampling is an observed night, NOT an exposure.
All variance terms are descriptive; no activity model or correlated-noise
likelihood is fitted. Weighted variances use population denominators.
"""
from __future__ import annotations

import numpy as np


def _validated(data: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    keys = ("bjd", "rvRelativeMetersPerSecond", "uncertaintyMetersPerSecond", "snr")
    arrays = [np.asarray(data[k], dtype=float) for k in keys]
    if any(a.ndim != 1 for a in arrays) or len({a.size for a in arrays}) != 1:
        raise ValueError("all input columns must be equally sized one-dimensional arrays")
    if arrays[0].size < 3 or any(not np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("at least three finite observations are required")
    if np.any(arrays[2] <= 0) or np.any(arrays[3] < 0):
        raise ValueError("uncertainties must be positive and S/N nonnegative")
    return np.floor(arrays[0]).astype(np.int64), arrays[1], arrays[2], arrays[3]


def variance_components(rv: np.ndarray, errors: np.ndarray, nights: np.ndarray) -> dict[str, float]:
    """Exactly partition weighted total RV variance into within/between-night terms.

    The per-night intercept is estimated from the same observations; no
    degrees-of-freedom correction is claimed for these descriptive variances.
    """
    if rv.size < 3 or np.unique(nights).size < 2:
        raise ValueError("need at least three exposures and two observed nights")
    w = 1.0 / errors**2
    overall = np.average(rv, weights=w)
    within = 0.0
    between = 0.0
    for night in np.unique(nights):
        sel = nights == night
        nw = np.sum(w[sel])
        mean = np.average(rv[sel], weights=w[sel])
        within += np.sum(w[sel] * (rv[sel] - mean) ** 2)
        between += nw * (mean - overall) ** 2
    denom = np.sum(w)
    total = float(np.sum(w * (rv - overall) ** 2) / denom)
    within /= denom
    between /= denom
    if not np.isclose(total, within + between, rtol=1e-11, atol=1e-12):
        raise ArithmeticError("weighted variance decomposition failed")
    return {
        "total_rms_mps": float(np.sqrt(total)),
        "within_night_rms_mps": float(np.sqrt(within)),
        "between_night_rms_mps": float(np.sqrt(between)),
        "within_variance_fraction": float(within / total) if total > 0 else 0.0,
    }


def quality_cut_decomposition(data: dict[str, np.ndarray], cuts=(10, 200, 300, 400)) -> list[dict]:
    """Compare S/N cuts after removing night means and reporting night composition."""
    nights, rv, errors, snr = _validated(data)
    rows = []
    for cut in cuts:
        if not np.isfinite(cut) or cut < 0:
            raise ValueError("S/N cuts must be finite and nonnegative")
        mask = snr >= cut
        if mask.sum() < 3 or np.unique(nights[mask]).size < 2:
            rows.append({"snr_cut": float(cut), "valid": False, "n": int(mask.sum())})
            continue
        result = variance_components(rv[mask], errors[mask], nights[mask])
        rows.append({
            "snr_cut": float(cut), "valid": True, "n": int(mask.sum()),
            "night_counts": {str(int(n)): int(np.sum(nights[mask] == n)) for n in np.unique(nights[mask])},
            **result,
        })
    return rows


def leave_one_night_out(data: dict[str, np.ndarray], snr_cut: float = 10) -> list[dict]:
    """Quantify leverage of each observed night on weighted RV scatter."""
    nights, rv, errors, snr = _validated(data)
    mask = snr >= snr_cut
    unique = np.unique(nights[mask])
    if unique.size < 3:
        raise ValueError("leave-one-night-out requires at least three nights")
    baseline = variance_components(rv[mask], errors[mask], nights[mask])
    result = []
    for night in unique:
        keep = mask & (nights != night)
        reduced = variance_components(rv[keep], errors[keep], nights[keep])
        result.append({
            "omitted_bjd_day": int(night), "n_retained": int(keep.sum()),
            "delta_total_rms_mps": reduced["total_rms_mps"] - baseline["total_rms_mps"],
            "delta_within_rms_mps": reduced["within_night_rms_mps"] - baseline["within_night_rms_mps"],
        })
    return result


def night_cluster_bootstrap(
    data: dict[str, np.ndarray], snr_cut: float = 10, n_resamples: int = 2000,
    seed: int = 20261029,
) -> dict:
    """Resample entire observed nights, preserving exposure-level dependence.

    Repeated selected nights are treated as separate blocks when estimating
    within-night scatter. Only 4 real nights exist; bootstrap percentiles are
    consequently exploratory and not calibrated confidence intervals.
    """
    nights, rv, errors, snr = _validated(data)
    selected = snr >= snr_cut
    unique = np.unique(nights[selected])
    if unique.size < 3 or n_resamples < 20:
        raise ValueError("need >=3 nights and >=20 resamples")
    blocks = [np.flatnonzero(selected & (nights == n)) for n in unique]
    rng = np.random.default_rng(seed)
    values = np.empty(n_resamples)
    for j in range(n_resamples):
        picks = rng.integers(0, len(blocks), size=len(blocks))
        idx = np.concatenate([blocks[i] for i in picks])
        # Repeated nights get new block IDs, avoiding false pooling.
        block_id = np.concatenate([np.full(blocks[i].size, k) for k, i in enumerate(picks)])
        w = 1.0 / errors[idx] ** 2
        resid = rv[idx].copy()
        for k in range(len(picks)):
            m = block_id == k
            resid[m] -= np.average(rv[idx][m], weights=w[m])
        values[j] = np.sqrt(np.sum(w * resid**2) / np.sum(w))
    point = variance_components(rv[selected], errors[selected], nights[selected])
    return {
        "n_nights": int(unique.size), "n_resamples": int(n_resamples),
        "within_night_rms_mps": point["within_night_rms_mps"],
        "bootstrap_p025_mps": float(np.quantile(values, 0.025)),
        "bootstrap_p975_mps": float(np.quantile(values, 0.975)),
        "resampling_unit": "observed night",
    }
