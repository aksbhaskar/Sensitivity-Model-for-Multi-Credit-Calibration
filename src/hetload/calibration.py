"""Bi-level calibration to CDX index tranches.

Outer loop: dependence parameters (lambda_c, gamma) against the summed absolute
tranche error. Inner loop: per-firm statics and daily states delta_i against
bootstrapped forward hazards. Loadings omega are fixed inputs.
"""

from __future__ import annotations


def calibrate(market_data, omega=None, **options):
    """Calibrate with loadings ``omega`` held fixed (``None`` means homogeneous)."""
    raise NotImplementedError("M0/M3: bi-level calibration")
