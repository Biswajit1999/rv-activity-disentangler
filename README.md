# RV Activity Disentangler

An interactive laboratory for testing when a low-amplitude Keplerian radial-velocity (RV) signal can survive a simple activity-indicator decorrelation.

[![CI](https://github.com/Biswajit1999/rv-activity-disentangler/actions/workflows/ci.yml/badge.svg)](https://github.com/Biswajit1999/rv-activity-disentangler/actions/workflows/ci.yml)
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**[Launch the interactive laboratory →](https://biswajit1999.github.io/rv-activity-disentangler/)**

## Motivation

Stellar photospheres generate RV structure through spots, faculae, convective-blueshift suppression, rotation, and magnetic cycles. A planet candidate is credible only if its Doppler signal can be distinguished from those processes and from the observing window. This repository provides a controlled synthetic experiment: the planet and activity truth are known, so the failure mode of linear decorrelation is visible rather than hidden.

## Research question

How do orbital period, stellar rotation period, activity amplitude, and proxy fidelity affect recovery of a planet's RV semi-amplitude?

## Implemented experiment

- circular single-planet injection;
- rotational activity fundamental plus first harmonic;
- an imperfect activity proxy with configurable fidelity;
- deterministic sub-metre-per-second measurement structure;
- ordinary-least-squares activity slope;
- raw and corrected RV trajectories;
- recovered semi-amplitude, residual RMS, separation gain, and period ratio;
- accessible data table, responsive UI, reduced-motion support;
- tests for finite output, useful-proxy improvement, and the zero-activity limit.

## Mathematical outline

```text
v_planet(t) = K sin(2πt/P_orb + φ)

v_activity(t) = A [sin(2πt/P_rot) + 0.32 sin(4πt/P_rot + 0.8)]

v_raw = v_planet + v_activity + ε

v_corrected = v_raw - β(I - mean(I))
```

`β` is fit by ordinary least squares against the synthetic indicator `I`. The recovered `K` is the quadrature amplitude at the known injected orbital period.

## Run locally

```bash
git clone https://github.com/Biswajit1999/rv-activity-disentangler.git
cd rv-activity-disentangler
npm install
npm run dev
```

## Verify

```bash
npm run check
```

This runs OXlint, Vitest, TypeScript, and the production Vite build. CI runs the same stages on GitHub.

## Repository structure

```text
src/science.ts       injection, OLS separation, recovery diagnostics
src/science.test.ts  numerical and scientific invariants
src/project.ts       experiment metadata and declared assumptions
src/App.tsx          interactive workbench with Motion transitions
src/Chart.tsx        SVG series and accessible table
docs/METHODS.md      derivation and research caveats
design-system/       persisted UI design contract
```

## Scientifically useful sweeps

- move `P_rot/P_orb` toward one and inspect leakage;
- lower indicator fidelity until correction becomes misleading;
- increase activity amplitude at fixed planetary `K`;
- compare an activity-free baseline before interpreting recovery;
- identify cases where residual RMS improves but recovered `K` becomes biased.

## What this demonstrates

The app demonstrates controlled signal injection, a transparent nuisance regression, and recovery diagnostics. It is useful for intuition, teaching, test-driven algorithm development, and defining harder experiments.

## What this does not demonstrate

It is not a planet detection pipeline and does not establish a false-alarm probability. It has no irregular cadence, seasonal gaps, heteroscedastic uncertainties, tellurics, instrumental offsets, eccentric or multiple planets, evolving active regions, or correlated stochastic process.

Linear decorrelation can also absorb a real planet when the proxy and planetary basis are correlated. That is a central result to investigate, not a nuisance to hide.

## Validation and next steps

1. Add irregular timestamp upload and instrument offsets.
2. Implement a generalized Lomb–Scargle window-function view.
3. Compare OLS with joint Keplerian + quasi-periodic Gaussian-process inference.
4. Run injection–recovery over phase, cadence, activity amplitude, and hyperparameters.
5. Report completeness and reliability rather than a single best fit.
6. Test on public HARPS/HARPS-N targets with provenance and published benchmarks.

## References

- Dumusque, X. et al. (2017), *Radial-velocity fitting challenge*, [A&A 598, A133](https://doi.org/10.1051/0004-6361/201628671).
- Haywood, R. D. et al. (2014), *Planets and stellar activity: hide and seek in the CoRoT-7 system*, [MNRAS 443, 2517](https://doi.org/10.1093/mnras/stu1320).
- Rajpaul, V., Aigrain, S. & Roberts, S. (2015), *A Gaussian process framework for modelling stellar activity signals in radial velocity data*, [MNRAS 452, 2269](https://doi.org/10.1093/mnras/stv1428).

## Citation and license

Citation metadata is in [`CITATION.cff`](CITATION.cff). Released under the [MIT License](LICENSE).
