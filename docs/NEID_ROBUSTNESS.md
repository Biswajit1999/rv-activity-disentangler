# NEID solar RV robustness: night-clustered uncertainty

Run `python scripts/analyze_neid_robustness.py` after installing `requirements-analysis.txt` and `scipy>=1.10`. This extends the existing 13-day NEID Level-2 observational-metadata audit without replacing it. It reads only rows already passing `analysisPass` and writes `research/neid-robustness-summary.json` (the JSON output is not claimed here until the runner is executed).

Three independent diagnostics:

1. **Profiled white jitter.** Maximize a Gaussian likelihood with variance `sigma_i²+s²` and one refitted weighted RV offset per `floor(BJD)` group. The objective includes `log(variance)`; choosing `s` only by reduced chi-square would be statistically inconsistent. The zero-jitter boundary is tested explicitly. Because night offsets are fitted, the estimate is conditional and potentially biased in short groups; it does not isolate stellar activity from instrumental effects.
2. **Whole-night bootstrap.** Compute formal-error-weighted means within each BJD-day group, then resample the group means with replacement. Report an equal-night-weight mean, percentile interval and bootstrap standard error, alongside the naive inverse-variance formal exposure error (not directly interchangeable estimands). Exposure-level resampling would underestimate correlation within nights.
3. **Leave-one-night-out sensitivity.** Recompute the sample standard deviation of the night means and the equal-night mean after omitting each observed day. This identifies undue dependence on a particular night, not a causal explanation.

All time units are BJD days and all velocities are m/s. The data span is far shorter than a solar rotation, and bootstrap coverage is conditional on observed nights being exchangeable. Results should not be interpreted as detections, calibrated RV precision or independently measured stellar jitter.
