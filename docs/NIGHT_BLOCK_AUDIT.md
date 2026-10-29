# NEID night-block robustness audit (29 October 2026)

Run `python scripts/run_night_block_audit.py` after installing `requirements-analysis.txt`.
The input is the repository's **existing public NEID solar metadata**, and the
script writes `research/night-block-audit-results.json` (not used by the browser).

Three independent questions:

1. **S/N selection and nightly composition:** weighted total RV variance is
   decomposed exactly into within-night and between-night terms, with
   \(w_i=1/\sigma_i^2\). Report retained counts by night at each cut; a
   reduction in pooled scatter can reflect changing night composition.
2. **Leave-one-night-out influence:** recompute weighted RMS after removing
   each complete night, not one highly dependent exposure at a time.
3. **Night-cluster bootstrap:** resample whole observed nights and recompute
   within-night weighted RMS, assigning distinct block IDs to repeated nights.

The four-night sample is **too small** for reliable frequentist uncertainty
coverage; bootstrap percentile endpoints are exploratory. Per-exposure formal
errors omit correlated calibration and solar effects. These diagnostics are
not a planet search, stellar activity correction, or period measurement.
No observational results are claimed until the runner executes against the
pinned input in a reproducible environment.
