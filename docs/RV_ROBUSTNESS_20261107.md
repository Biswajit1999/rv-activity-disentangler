# Robust RV sensitivity (2026-11-07)

Two synthetic-only numerical diagnostics: Huber IRLS location under known per-exposure formal scales, and leave-one-night-out influence. The robust location is not a corrected solar RV or a precision estimate. Standardized residuals assume independent scales; shared instrumental or solar systematics are omitted. Reproduce: `python -m unittest discover -s tests -p 'test_rv_robustness_20261107.py'`. No NEID observational results were regenerated. AI-assisted implementation; review before scientific use.
