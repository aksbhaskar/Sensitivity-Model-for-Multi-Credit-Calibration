"""Phase-type jump laws and the per-firm marginal jump mixture.

A phase-type law PH(alpha, T) is the absorption time of a continuous-time Markov
chain with initial distribution ``alpha`` over transient states and sub-generator
``T``; the exit-rate vector is ``t = -T 1`` and the Laplace transform is
``alpha (theta I - T)^{-1} t``.

Two facts make firm-specific loadings tractable:

* scaling: if Z ~ PH(alpha, T) and c > 0 then c Z ~ PH(alpha, T / c);
* mixing: superposing the idiosyncratic jump stream (rate lambda_J, Exp(eta))
  with the loaded common stream (rate lambda_c, omega_i * Exp(1/gamma)) gives a
  compound Poisson process whose jump law is a two-phase hyperexponential.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PhaseType:
    """Phase-type law PH(alpha, T) with no atom at zero."""

    alpha: np.ndarray
    T: np.ndarray

    def __post_init__(self) -> None:
        alpha = np.asarray(self.alpha, dtype=float)
        T = np.asarray(self.T, dtype=float)
        if T.ndim != 2 or T.shape[0] != T.shape[1]:
            raise ValueError("T must be a square matrix")
        if alpha.shape != (T.shape[0],):
            raise ValueError("alpha must have one entry per phase")
        if not np.isclose(alpha.sum(), 1.0):
            raise ValueError("alpha must sum to one (no atom at zero)")
        object.__setattr__(self, "alpha", alpha)
        object.__setattr__(self, "T", T)

    @classmethod
    def exponential(cls, rate: float) -> "PhaseType":
        """Exp(rate), the one-phase case."""
        return cls(np.array([1.0]), np.array([[-rate]]))

    @classmethod
    def hyperexponential(cls, weights, rates) -> "PhaseType":
        """Mixture sum_k weights[k] Exp(rates[k])."""
        return cls(np.asarray(weights, dtype=float), -np.diag(np.asarray(rates, dtype=float)))

    @property
    def exit_rates(self) -> np.ndarray:
        """Exit-rate vector t = -T 1."""
        return -self.T.sum(axis=1)

    def mean(self) -> float:
        """E[Z] = -alpha T^{-1} 1."""
        ones = np.ones(self.T.shape[0])
        return float(-self.alpha @ np.linalg.solve(self.T, ones))

    def laplace(self, theta: complex) -> complex:
        """Laplace transform E[exp(-theta Z)] = alpha (theta I - T)^{-1} t."""
        n = self.T.shape[0]
        return complex(self.alpha @ np.linalg.solve(theta * np.eye(n) - self.T, self.exit_rates))

    def scale(self, c: float) -> "PhaseType":
        """Law of c Z, which is PH(alpha, T / c)."""
        if c <= 0:
            raise ValueError("scale factor must be positive")
        return PhaseType(self.alpha, self.T / c)

    def sample(self, size: int, rng: np.random.Generator | None = None) -> np.ndarray:
        """Draw samples by simulating the underlying Markov chain."""
        rng = np.random.default_rng() if rng is None else rng
        n = self.T.shape[0]
        rates = -np.diag(self.T)
        # Row i: probabilities of moving to each transient state, then absorption.
        jump = np.zeros((n, n + 1))
        jump[:, :n] = self.T / rates[:, None]
        np.fill_diagonal(jump[:, :n], 0.0)
        jump[:, n] = self.exit_rates / rates

        out = np.empty(size)
        for k in range(size):
            state = rng.choice(n, p=self.alpha)
            total = 0.0
            while state < n:
                total += rng.exponential(1.0 / rates[state])
                state = rng.choice(n + 1, p=jump[state])
            out[k] = total
        return out


def marginal_jump_law(
    lambda_j: float, eta: float, lambda_c: float, gamma: float, omega: float = 1.0
) -> tuple[float, PhaseType]:
    """Total jump intensity and jump law of firm i's driver X^i + omega_i Y.

    Idiosyncratic jumps arrive at rate ``lambda_j`` with Exp(eta) sizes; common
    jumps arrive at rate ``lambda_c`` with mean ``gamma`` and are scaled by the
    firm's loading ``omega``. With ``omega = 1`` this is eq. (11) of Baker and
    Capponi. The result plugs straight into the single-name pricing theorem.
    """
    if omega <= 0:
        raise ValueError("loading omega must be positive")
    total = lambda_j + lambda_c
    law = PhaseType.hyperexponential(
        weights=[lambda_j / total, lambda_c / total],
        rates=[eta, 1.0 / (omega * gamma)],
    )
    return total, law
