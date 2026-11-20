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
    def exponential(cls, rate: float) -> PhaseType:
        """Exp(rate), the one-phase case."""
        return cls(np.array([1.0]), np.array([[-rate]]))

    @classmethod
    def hyperexponential(cls, weights, rates) -> PhaseType:
        """Mixture sum_k weights[k] Exp(rates[k])."""
        return cls(np.asarray(weights, dtype=float), -np.diag(np.asarray(rates, dtype=float)))

    @property
    def exit_rates(self) -> np.ndarray:
        """Exit-rate vector t = -T 1."""
        return -self.T.sum(axis=1)

    def mean(self) -> float:
        """E[Z] = -alpha T^{-1} 1."""
        rates = -np.diag(self.T)
        if np.array_equal(self.T, -np.diag(rates)):
            return float(np.sum(self.alpha / rates))
        ones = np.ones(self.T.shape[0])
        return float(-self.alpha @ np.linalg.solve(self.T, ones))

    def laplace(self, theta: complex | np.ndarray) -> complex | np.ndarray:
        """Laplace transform E[exp(-theta Z)] = alpha (theta I - T)^{-1} t."""
        n = self.T.shape[0]
        theta = np.asarray(theta, dtype=complex)
        rates = -np.diag(self.T)
        if np.array_equal(self.T, -np.diag(rates)):
            result = np.sum(self.alpha * rates / (theta[..., None] + rates), axis=-1)
        else:
            matrices = theta[..., None, None] * np.eye(n) - self.T
            rhs = np.broadcast_to(self.exit_rates, theta.shape + (n,))[..., None]
            result = np.linalg.solve(matrices, rhs)[..., 0] @ self.alpha
        return complex(result) if theta.ndim == 0 else result

    def scale(self, c: float) -> PhaseType:
        """Law of c Z, which is PH(alpha, T / c)."""
        if c <= 0:
            raise ValueError("scale factor must be positive")
        return PhaseType(self.alpha, self.T / c)

    def sample(self, size: int, rng: np.random.Generator | None = None) -> np.ndarray:
        """Draw samples by simulating the underlying Markov chain."""
        rng = np.random.default_rng() if rng is None else rng
        n = self.T.shape[0]
        rates = -np.diag(self.T)
        if n == 1:
            return rng.exponential(1.0 / rates[0], size=size)
        if np.array_equal(self.T, -np.diag(rates)):
            state = rng.choice(n, size=size, p=self.alpha)
            return rng.exponential(1.0 / rates[state])
        # Row i: probabilities of moving to each transient state, then absorption.
        jump = np.zeros((n, n + 1))
        jump[:, :n] = self.T / rates[:, None]
        np.fill_diagonal(jump[:, :n], 0.0)
        jump[:, n] = self.exit_rates / rates

        cumulative = np.cumsum(jump, axis=1)
        cumulative[:, -1] = 1.0
        out = np.zeros(size)
        state = rng.choice(n, size=size, p=self.alpha)
        active = np.arange(size)
        while active.size:
            current = state[active]
            out[active] += rng.exponential(1.0 / rates[current])
            uniforms = rng.random(active.size)
            state[active] = np.sum(uniforms[:, None] >= cumulative[current], axis=1)
            active = active[state[active] < n]
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
