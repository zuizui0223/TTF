from __future__ import annotations

import numpy as np


def _validated_probability_vector(value: np.ndarray, *, atol: float = 1e-8) -> np.ndarray:
    x = np.asarray(value, dtype=float)
    if x.shape != (4,):
        raise ValueError("colour vector must have exactly four components")
    if not np.isfinite(x).all() or np.any(x < 0.0):
        raise ValueError("colour vector must be finite and nonnegative")
    if not np.isclose(float(x.sum()), 1.0, atol=atol, rtol=0.0):
        raise ValueError("colour vector must sum to one within the frozen tolerance")
    return x


def four_component_colour_jsd(a: np.ndarray, b: np.ndarray) -> float:
    """Jensen-Shannon divergence for the frozen four-component soft colour vector.

    Natural logarithms are used. Zero-probability terms contribute exactly zero.
    Inputs are validated but never renormalized; invalid measured vectors must remain
    terminal missingness under the empirical support firewall rather than being
    silently repaired here.
    """

    p = _validated_probability_vector(a)
    q = _validated_probability_vector(b)
    m = 0.5 * (p + q)

    def kl_term(x: np.ndarray) -> float:
        mask = x > 0.0
        return float(np.sum(x[mask] * np.log(x[mask] / m[mask])))

    value = 0.5 * kl_term(p) + 0.5 * kl_term(q)
    # Roundoff can produce a tiny negative value around exact equality.
    if value < 0.0 and value >= -1e-15:
        value = 0.0
    if not np.isfinite(value) or value < 0.0:
        raise RuntimeError("Jensen-Shannon divergence produced an invalid value")
    return float(value)
