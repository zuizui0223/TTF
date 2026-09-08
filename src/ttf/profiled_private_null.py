from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .private_null_inference import upper_monte_carlo_pvalue


def robust_profile_distance(
    value: float,
    reference: np.ndarray,
    *,
    scale_floor: float = 1e-6,
) -> tuple[float, float, float]:
    """Robust distance from a training-only nuisance summary to one null model."""
    x = np.asarray(reference, dtype=float)
    if x.ndim != 1 or len(x) < 20 or not np.isfinite(x).all():
        raise ValueError("profile reference must be a finite 1D array with >=20 draws")
    if scale_floor <= 0:
        raise ValueError("scale_floor must be positive")
    median = float(np.median(x))
    mad = float(np.median(np.abs(x - median)))
    scale = max(1.4826 * mad, float(scale_floor))
    return float(abs(float(value) - median) / scale), median, scale


def profiled_private_pvalue(
    observed_statistic: float,
    observed_strength: float,
    references: Mapping[str, tuple[np.ndarray, np.ndarray]],
    *,
    profile_draws: int = 999,
    selected_configs: int = 2,
    scale_floor: float = 1e-6,
) -> tuple[float, tuple[str, ...], dict[str, float], dict[str, float]]:
    """Profile nuisance strength using training-only coherence, then test raw T.

    Each reference entry is ``(strength, statistic)`` from the same ordered set
    of private-null worlds.  The first ``profile_draws`` strengths are used only
    to rank nuisance configurations.  Statistics from the remaining worlds are
    used only for Monte Carlo p-values.  Held-out transfer ``T`` therefore never
    chooses the nuisance configuration.

    The returned p-value is the maximum component p-value among the fixed number
    of most strength-compatible configurations.
    """
    if not references:
        raise ValueError("at least one profiled private reference is required")
    if selected_configs < 1 or selected_configs > len(references):
        raise ValueError("selected_configs is out of range")

    distances: dict[str, float] = {}
    calibration: dict[str, np.ndarray] = {}
    for label, pair in references.items():
        strength = np.asarray(pair[0], dtype=float)
        statistic = np.asarray(pair[1], dtype=float)
        if strength.shape != statistic.shape or strength.ndim != 1:
            raise ValueError(f"reference shape mismatch for {label}")
        if not profile_draws < len(strength):
            raise ValueError("profile_draws must leave independent calibration draws")
        distance, _, _ = robust_profile_distance(
            observed_strength,
            strength[:profile_draws],
            scale_floor=scale_floor,
        )
        distances[str(label)] = distance
        calibration[str(label)] = statistic[profile_draws:]

    selected = tuple(
        sorted(distances, key=lambda label: (distances[label], label))[:selected_configs]
    )
    component = {
        label: upper_monte_carlo_pvalue(observed_statistic, calibration[label])
        for label in selected
    }
    p = float(max(component.values()))
    return p, selected, component, distances
