"""Robust RV location and grouped-night sensitivity; no new NEID measurement."""
from __future__ import annotations
import math

def weighted_mean(values, errors):
    if len(values) != len(errors) or not values:
        raise ValueError("matching nonempty arrays required")
    if any(not math.isfinite(x) for x in values) or any(not math.isfinite(e) or e <= 0 for e in errors):
        raise ValueError("finite values and positive errors required")
    weights = [1 / e**2 for e in errors]
    return sum(x*w for x,w in zip(values,weights)) / sum(weights)

def huber_location(values, errors, cutoff=1.5, iterations=200, tolerance=1e-11):
    """IRLS Huber location with standardized residuals (x-mu)/sigma.
    
    Errors are treated as fixed known independent Gaussian scales; robust
    weights limit influence of large standardized residuals. No covariance
    or uncertainty interval is estimated.
    """
    mu = weighted_mean(values, errors)
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("positive finite cutoff required")
    for _ in range(iterations):
        w = [(1/e**2)*min(1., cutoff*e/max(abs(x-mu), 1e-300))
             for x,e in zip(values,errors)]
        new = sum(x*wi for x,wi in zip(values,w))/sum(w)
        if abs(new-mu) <= tolerance*max(1.,abs(mu)):
            return new
        mu = new
    raise RuntimeError("Huber IRLS did not converge")

def leave_one_group_out(values, errors, groups, cutoff=1.5):
    """Return group-deletion Huber estimates; require >=2 groups and >=1 retained point."""
    if len(groups) != len(values):
        raise ValueError("group length mismatch")
    unique = list(dict.fromkeys(groups))
    if len(unique) < 2:
        raise ValueError("at least two groups required")
    return {str(g): huber_location([v for v,k in zip(values,groups) if k != g],
                                   [e for e,k in zip(errors,groups) if k != g],cutoff)
            for g in unique}
