"""Wiener-Hopf Monte Carlo for the common-jump multi-credit model.

Segments of length Exp(q + lambda_c) end either in a pricing step or, with
probability lambda_c / (q + lambda_c), in a common arrival. At a common arrival
a single E ~ Exp(1/gamma) is drawn and firm i's position moves by omega_i * E
(omega = 1 for every firm recovers the homogeneous scheme).
"""

from __future__ import annotations


def simulate_default_times(
    drivers, deltas, lambda_c, gamma, omega, horizon, n_paths, q=16.0, rng=None
):
    """Simulate joint default times for all firms in the portfolio."""
    raise NotImplementedError("M0/M1: Wiener-Hopf Monte Carlo")
