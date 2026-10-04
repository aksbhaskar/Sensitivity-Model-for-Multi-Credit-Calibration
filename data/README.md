# Data

No market data is committed to this repository (see `.gitignore`). Put local copies under
`data/raw/` and derived files under `data/processed/`.

| Dataset | Source | Used for |
|---|---|---|
| Single-name CDS composites (CDX-NAHY, CDX-NAIG constituents) | S&P Global / Markit via WRDS | forward-hazard bootstrap, single-name calibration |
| Index tranche upfronts (5Y, on-the-run) | S&P Global via WRDS | joint calibration target |
| Constant-maturity Treasury par yields | FRED | discount curve |
| Equity prices / returns | CRSP via WRDS, or public sources | P-measure loadings |
| Default-correlation studies by sector and rating | rating-agency annual default studies | P-measure loadings (bucket prior) |

Reference sample in the original paper: 2021-05-27 to 2025-02-13. HY recovery is 0.30 and
the HY coupon is 500 bp; IG recovery is 0.40 and the IG coupon is 100 bp.

WRDS and S&P Global data are licensed. Do not commit or share raw extracts.
