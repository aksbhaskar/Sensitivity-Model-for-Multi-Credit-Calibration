# Research plan

## Question

In Baker & Capponi's multi-credit model, every firm receives the same common compound Poisson
jump. We give firm *i* a loading ω<sub>i</sub> > 0, so its driver becomes X̃<sup>i</sup> + ω<sub>i</sub>Y.
We then ask whether loadings estimated under the historical measure P, and held fixed, improve
out-of-sample tranche pricing without the convergence and over-parameterisation costs of
calibrating the loadings under Q.

## What is already settled

1. **Single-name tractability.** Firm *i*'s jump law is the two-phase mixture
   λ<sub>J</sub>/(λ<sub>J</sub>+λ<sub>c</sub>)·Exp(η) + λ<sub>c</sub>/(λ<sub>J</sub>+λ<sub>c</sub>)·Exp(1/(ω<sub>i</sub>γ)).
   This is phase-type with m = 2, so Theorem 2.3 applies unchanged (`hetload.marginal_jump_law`).
2. **Simulation.** The WHMC scheme only changes at common arrivals: one shared draw E, and firm *i*
   moves by ω<sub>i</sub>E.
3. **Identification.** Only ω<sub>i</sub>γ is identified, so we normalise mean(ω) = 1.

Still to check:

- Talbot-contour safety. The roots β<sub>j</sub>(q) now depend on ω<sub>i</sub>γ. Check that
  fitted ranges stay clear of root coalescence and of β<sub>j</sub>(q) ∈ {0, 1} (Remark 2.5).
- The subordinator driver. Its closed-form marginal must use the loaded mixture.

## Estimating ω under P

Candidate estimators, from cheapest to richest. All are normalised to mean 1 and floored at a
small positive value.

| Estimator | Input | Notes |
|---|---|---|
| Sector/rating lookup | historical default-correlation studies by sector and rating | coarse, no market data needed |
| Equity factor beta | daily equity returns vs. a market or credit-sensitive index | firm-level; needs public tickers |
| Equity tail beta | co-exceedance of large negative returns with the index | targets jump-type co-movement specifically |
| First principal component | panel of equity returns or realised-volatility changes | data-driven; check sign and stability |

Estimation windows end **before** each calibration fold, so loadings are never fitted on
test-period data.

## Experiments

For each board (NAHY Series 36 cohort, NAHY full basket, NAIG full basket) and each driver
(LBM+Exp, Cramér–Lundberg, subordinator):

1. **Homogeneous**: ω ≡ 1. This must reproduce the paper's Table 12 within Monte Carlo noise.
2. **Free buckets (Q)**: ω by sector and/or rating, calibrated in the outer loop.
3. **Fixed prior (P)**: ω from each estimator above, frozen.

Metrics, per fold and as the mean over folds:

- summed absolute tranche error at the *frozen*, *frequency*, and *full* marking levels;
- median single-name hazard RMSE (it must not get worse);
- outer-loop iterations, wall time, and the condition number / smallest Hessian eigenvalue
  of the outer objective at the optimum (to measure flat directions);
- stability of fitted (λ<sub>c</sub>, γ) across folds.

## Success criteria

The fixed-prior specification
(a) has a lower out-of-sample tranche error than homogeneous on at least two of three boards,
(b) matches or beats the free-bucket version out of sample, and
(c) converges in about the same number of outer iterations as homogeneous.

If (a) fails, the result is still informative: under Q, the shared-jump assumption is not the
binding constraint.
