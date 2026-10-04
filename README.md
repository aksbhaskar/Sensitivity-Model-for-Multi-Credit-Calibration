# Sensitivity Model for Multi-Credit Calibration

[![arXiv](https://img.shields.io/badge/arXiv-2608.10321-b31b1b.svg)](https://arxiv.org/abs/2608.10321)
[![CI](https://github.com/aksbhaskar/sensitivity-model-for-multi-credit-calibration/actions/workflows/ci.yml/badge.svg)](https://github.com/aksbhaskar/sensitivity-model-for-multi-credit-calibration/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**Firm-specific sensitivities to the common jump factor in the elastically stopped Lévy multi-credit model.**

This project extends Baker & Capponi,
[*Multi-Credit Calibration via Elastically Stopped Lévy Processes*](https://arxiv.org/abs/2608.10321) (2026).
That paper models each firm's default as the time when the running supremum of a latent,
spectrally positive distress process crosses an independent exponential barrier. Adding **one common
compound Poisson jump factor** to every firm's driver gives a parsimonious joint model with
simultaneous defaults. The model is priced exactly by a Wiener–Hopf Monte Carlo scheme and
calibrated to CDX North American High-Yield and Investment-Grade tranches.

In the baseline, **every firm receives the same common jump**. Real portfolios are not uniformly
exposed to systemic shocks. This repository gives each firm *i* its own sensitivity
**ω<sub>i</sub>** to the common factor. The loadings are **estimated under the historical measure P
and held fixed**, so the risk-neutral calibration gets no new free parameters.

---

## The model

**Baseline (Baker & Capponi, eq. 10).** Firm *i* has an idiosyncratic spectrally positive driver
X̃<sup>i</sup> (LBM+Exp, Cramér–Lundberg, or subordinator), a distance-to-distress δ<sub>i</sub> ≥ 0,
and an independent barrier ξ<sub>i</sub> ~ Exp(1). Y is a compound Poisson process with rate
λ<sub>c</sub> and Exp(1/γ) jumps (mean γ), shared by all firms:

```
Λ^i_t = sup_{s≤t} ( X̃^i_s + Y_s ),        τ_i = inf{ t ≥ 0 : Λ^i_t ≥ δ_i + ξ_i }
```

**This project.** Each firm takes its own loading on the common factor:

```
Λ^i_t = sup_{s≤t} ( X̃^i_s + ω_i Y_s ),    ω_i > 0
```

At each common arrival, firm *i*'s distress jumps by ω<sub>i</sub>·E, where E ~ Exp(1/γ) is a single
draw shared by all firms. Simultaneous defaults are preserved, but the systemic shock now hits
different firms with different force.

### Tractability survives

Seen from a single firm, the jumps of X̃<sup>i</sup> + ω<sub>i</sub>Y are still a two-component
hyperexponential mixture:

```
λ_J/(λ_J+λ_c) · Exp(η)  +  λ_c/(λ_J+λ_c) · Exp(1/(ω_i γ))
```

This is phase-type with *m* = 2, because ω<sub>i</sub>·Exp(1/γ) = Exp(1/(ω<sub>i</sub>γ)). So:

- **Single-name pricing is unchanged.** Theorem 2.3's finite partial-fraction formula applies to each
  firm as written, with d = 4 for LBM+Exp and d = 3 for Cramér–Lundberg. Only the common component's
  mean changes, from γ to ω<sub>i</sub>γ, so the root system becomes firm-specific. It is already
  firm-specific through each firm's own static parameters, so the computational cost stays the same.
- **Wiener–Hopf Monte Carlo is unchanged.** Common arrivals still form an exponential clock. At an
  arrival, draw one E ~ Exp(1/γ) and add ω<sub>i</sub>E to firm *i*, instead of E to every firm.
- **Identification needs a normalisation.** Only the products ω<sub>i</sub>γ enter the joint law, so
  the scale of ω is confounded with γ. We fix mean(ω) = 1, which nests the baseline at ω ≡ 1.

So the open question is empirical rather than mathematical. Free per-firm (or per-sector/rating)
loadings calibrated under Q add parameters and flat directions to the bi-level calibration, and they
risk overfitting. **Can P-measure loadings, fixed in advance, capture real heterogeneity and beat the
homogeneous model out of sample without those costs?**

## Specifications compared

| | Homogeneous (baseline) | Free loadings (Q) | **Fixed-prior loadings (P)** |
|---|---|---|---|
| Common-jump size for firm *i* | Exp(1/γ) | Exp(1/(ω<sub>i</sub>γ)), ω calibrated to tranches | Exp(1/(ω<sub>i</sub>γ)), **ω estimated from history, frozen** |
| Dependence parameters (Q) | (λ<sub>c</sub>, γ) | (λ<sub>c</sub>, γ) + sector/rating or per-firm ω | (λ<sub>c</sub>, γ) |
| Expected calibration behaviour | reference | slower convergence, flat directions | same cost as baseline |

## Evaluation protocol

We follow the paper's design so the results are directly comparable:

- **Boards:** CDX-NAHY Series 36 cohort, CDX-NAHY full basket, CDX-NAIG full basket.
- **Single names:** forward hazards bootstrapped ISDA-style from Markit CDS composites; per-firm statics on
  a trailing 3-month window, re-fitted quarterly; daily state δ<sub>i</sub>.
- **Tranches:** 5Y on-the-run upfronts (0/15/25/35/100% HY, 0/3/7/15/100% IG); target = summed
  absolute tranche error; 10 series-aligned rolling folds per board (3m fit / 3m test).
- **Marking levels:** *frozen*, *frequency* (re-mark λ<sub>c</sub>), and *full* (re-mark λ<sub>c</sub>, γ).
- **Simulation:** Wiener–Hopf Monte Carlo at q = 16 with common random numbers.
- **Report:** out-of-sample tranche error, single-name RMSE (to check the marginals are not
  corrupted), convergence iterations and wall time, and the conditioning of the outer objective.

## Roadmap

- [ ] **M0 — Baseline replication.** Phase-type root system and Theorem 2.3 transform, Talbot inversion,
      single-name panel calibration, WHMC, and bi-level tranche calibration with homogeneous jumps.
- [ ] **M1 — Loadings in the pricer.** Firm-specific ω<sub>i</sub> in the marginal mixture and the WHMC
      common-arrival step; unit tests showing ω ≡ 1 reproduces the baseline exactly.
- [ ] **M2 — P-measure loadings.** Estimate ω from public historical data: equity-return and default
      co-movement, and sector/rating default-correlation studies. Normalise to mean 1.
- [ ] **M3 — Calibration.** Recalibrate (λ<sub>c</sub>, γ) and the statics with ω frozen; also run the
      free-loading benchmark (sector/rating buckets).
- [ ] **M4 — Evaluation and write-up.** Compare all three specifications on all three boards.

See [`docs/research-plan.md`](docs/research-plan.md) for details.

## Repository layout

```
.
├── src/hetload/          # Python package
│   ├── phase_type.py     # phase-type laws, scaling, the per-firm marginal jump mixture
│   ├── loadings.py       # P-measure estimation and normalisation of loadings ω
│   ├── simulation.py     # Wiener–Hopf Monte Carlo with heterogeneous common jumps
│   └── calibration.py    # bi-level calibration to CDX index tranches
├── tests/                # pytest suite
├── notebooks/            # exploratory analysis and figures
├── data/                 # data instructions (no market data is committed)
└── docs/                 # research plan and derivations
```

## Getting started

Requires Python ≥ 3.10.

```bash
git clone https://github.com/aksbhaskar/sensitivity-model-for-multi-credit-calibration.git
cd sensitivity-model-for-multi-credit-calibration
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
python -m pytest
```

## Data

The original study uses S&P Global / Markit CDS and tranche data via WRDS, and FRED Treasury
par yields. None of this is redistributed here. See [`data/README.md`](data/README.md).

## Contributing

Work on feature branches and open pull requests into `main`. See [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Team

- **Akshat Bhaskar** ([@aksbhaskar](https://github.com/aksbhaskar))
- **Ishaan Harish**
- **Saksham Arora** ([@saksham10arora-dotcom](https://github.com/saksham10arora-dotcom))

## Reference

```bibtex
@misc{baker2026multicredit,
  title         = {Multi-Credit Calibration via Elastically Stopped L\'evy Processes},
  author        = {Baker, Graeme and Capponi, Agostino},
  year          = {2026},
  eprint        = {2608.10321},
  archivePrefix = {arXiv},
  primaryClass  = {q-fin.MF},
  url           = {https://arxiv.org/abs/2608.10321}
}
```

Citation metadata for this repository is in [`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE).
