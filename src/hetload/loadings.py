"""Firm-level loadings omega_i on the common jump factor, estimated under P.

Only the products omega_i * gamma enter the joint law, so the scale of the
loadings is not identified separately from the common jump mean gamma. Loadings
are therefore normalised to mean one, which nests the homogeneous model at
omega = 1 and keeps gamma comparable across specifications.
"""

from __future__ import annotations

import numpy as np


def normalise_loadings(raw) -> np.ndarray:
    """Rescale positive raw loadings to have mean one."""
    raw = np.asarray(raw, dtype=float)
    if np.any(raw <= 0):
        raise ValueError("loadings must be positive")
    return raw / raw.mean()


def estimate_loadings(panel, method: str = "first_pc") -> np.ndarray:
    """Estimate loadings from a (date x firm) panel of historical co-movement data."""
    raise NotImplementedError("M2: P-measure loading estimation")
