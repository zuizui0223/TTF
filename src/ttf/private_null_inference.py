from __future__ import annotations

from collections.abc import Mapping

import numpy as np


def upper_monte_carlo_pvalue(
    observed: float,
    reference: np.ndarray,
) -> float:
    """Exact-style one-sided Monte Carlo rank p-value against one reference null.

    The reference draws must be generated independently of the observed statistic
    under the null configuration. The +1 correction makes the finite-reference
    test conservative rather than allowing a zero p-value.
    """
    value = float(observed)
    ref = np.asarray(reference, dtype=float)
    if ref.ndim != 1 or len(ref) < 99 or not np.isfinite(ref).all():
        raise ValueError("reference must be a finite 1D array with at least 99 draws")
    if not np.isfinite(value):
        raise ValueError("observed statistic must be finite")
    return float((1 + np.count_nonzero(ref >= value)) / (len(ref) + 1))


def envelope_upper_pvalue(
    observed: float,
    references: Mapping[str, np.ndarray],
) -> tuple[float, str, dict[str, float]]:
    """Worst-case upper-tail p-value across a frozen private-null family.

    Returns ``(p_envelope, least_favourable_label, component_pvalues)``. Taking
    the maximum p-value tests against the union of the supplied private-null
    configurations without selecting a favourable nuisance setting after seeing
    the statistic.
    """
    if not references:
        raise ValueError("at least one private-null reference is required")
    component = {
        str(label): upper_monte_carlo_pvalue(observed, values)
        for label, values in references.items()
    }
    least = max(sorted(component), key=lambda label: component[label])
    return float(component[least]), least, component


def envelope_upper_pvalues(
    observed: np.ndarray,
    references: Mapping[str, np.ndarray],
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Vectorized convenience wrapper for a sequence of observed statistics."""
    values = np.asarray(observed, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("observed must be a finite 1D array")
    p = np.empty(len(values), dtype=float)
    labels: list[str] = []
    for i, value in enumerate(values):
        p[i], label, _ = envelope_upper_pvalue(float(value), references)
        labels.append(label)
    return p, tuple(labels)
