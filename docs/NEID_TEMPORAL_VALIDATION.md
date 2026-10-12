# NEID solar RV temporal validation

This analysis operates on the 635 accepted rows in the committed, provenance-
preserving NEID Level-2 solar metadata bundle. These are real public archive
measurements; no new observations are claimed.

Three complementary diagnostics are generated:

1. Nightly inverse-variance means and reduced chi-square. For measurements
   `v_i` with formal errors `sigma_i`, weights are `w_i = 1/sigma_i^2` and the
   mean uncertainty is `1/sqrt(sum(w_i))`. Reduced chi-square much larger than
   one means the formal errors do not describe the total within-night scatter;
   it does not identify whether the excess is solar or instrumental.
2. An irregular-cadence first-order structure function, reported as the median
   `|v_j-v_i|/sqrt(2)` in logarithmic time-lag bins. Every unordered exposure
   pair is used exactly once.
3. A predeclared S/N-threshold sweep that recomputes sample size, RV RMS, and
   nuisance correlations. This reveals selection sensitivity instead of
   presenting the original quality cut as uniquely correct.

Run `python scripts/analyze_neid_variability.py`. The 13-day window is far too
short for a solar-rotation measurement, formal archive errors do not include
stellar variability or shared systematics, and Pearson correlations with S/N
or solar altitude must not be interpreted causally.

