# Methods

## Sampling and truth

The model evaluates 481 evenly spaced epochs from 0 to 120 days at a six-hour cadence. Uniform sampling intentionally removes the spectral window so the first experiment isolates activity correlation.

## Activity proxy

The selected fidelity multiplies the true activity term. The remaining fraction is filled by a deterministic, differently phased oscillation. This makes every UI state reproducible while allowing the proxy to become incomplete.

## OLS coefficient

The slope is the centred covariance of indicator and RV divided by indicator variance. An epsilon guard covers the zero-activity limit. The correction subtracts only the centred proxy term, preserving the RV mean.

## Semi-amplitude recovery

Sine and cosine projections at the *known injected period* provide a compact recovery diagnostic. This is not a blind period search. A real analysis must fit period, phase, eccentricity, offsets, jitter, and activity jointly and propagate posterior uncertainty.

## Validation target

The recommended extension is a seeded injection–recovery matrix with held-out phases and irregular public observation timestamps. Recovery criteria should be declared before inspecting results.
