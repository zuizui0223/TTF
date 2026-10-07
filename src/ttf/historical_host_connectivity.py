from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .lgm_climate_envelope import CommonClimateSpace, gaussian_kde_density, whiten_environment


RESISTANCE_EPSILON = 1e-6


def host_kernel_affinity(
    training_current_environment: np.ndarray,
    query_environment: np.ndarray,
    space: CommonClimateSpace,
    *,
    chunk_size: int = 4096,
) -> np.ndarray:
    training = whiten_environment(np.asarray(training_current_environment, float), space)
    query = whiten_environment(np.asarray(query_environment, float), space)
    out = gaussian_kde_density(training, query, chunk_size=chunk_size)
    if np.any(out < 0.0) or np.any(out > 1.0 + 1e-12):
        raise RuntimeError("host kernel affinity left frozen [0,1] range")
    return np.clip(out, 0.0, 1.0)


def host_union(suitabilities: np.ndarray) -> np.ndarray:
    x = np.asarray(suitabilities, dtype=float)
    if x.ndim != 2 or x.shape[0] < 1 or x.shape[1] < 1:
        raise ValueError("host suitability must be host x point")
    if not np.isfinite(x).all() or np.any(x < 0.0) or np.any(x > 1.0):
        raise ValueError("host suitability must be finite in [0,1]")
    return 1.0 - np.prod(1.0 - x, axis=0)


def resistance_from_union(
    union_suitability: np.ndarray,
    *,
    epsilon: float = RESISTANCE_EPSILON,
) -> np.ndarray:
    x = np.asarray(union_suitability, dtype=float)
    if epsilon <= 0.0 or epsilon >= 1.0:
        raise ValueError("epsilon must lie in (0,1)")
    if not np.isfinite(x).all() or np.any(x < 0.0) or np.any(x > 1.0):
        raise ValueError("union suitability must be finite in [0,1]")
    return -np.log(epsilon + x)


def edge_resistance(point_union_suitability: np.ndarray, *, points_per_edge: int = 5) -> np.ndarray:
    x = np.asarray(point_union_suitability, dtype=float)
    if x.ndim != 2 or x.shape[1] != int(points_per_edge):
        raise ValueError("point suitability must be edge x frozen segment-point count")
    return resistance_from_union(x).mean(axis=1)


def zscore_population(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("predictor must be a finite vector with >=2 edges")
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0.0:
        raise ValueError("predictor has zero/non-finite population SD")
    return (x - float(np.mean(x))) / sd


@dataclass(frozen=True)
class SpeciesHistoricalInformation:
    species: str
    edges: int
    unique_fraction: float
    signal_share: float = float("nan")


def historical_unique_fraction(lgm: np.ndarray, present: np.ndarray) -> float:
    z_lgm = zscore_population(lgm)
    z_now = zscore_population(present)
    design = np.column_stack([np.ones(len(z_now)), z_now])
    beta, *_ = np.linalg.lstsq(design, z_lgm, rcond=None)
    residual = z_lgm - design @ beta
    q = float(np.mean(residual * residual))
    if q < -1e-12 or q > 1.0 + 1e-10:
        raise RuntimeError("historical unique fraction left expected [0,1] range")
    return float(np.clip(q, 0.0, 1.0))


def information_gate(
    by_species: dict[str, tuple[np.ndarray, np.ndarray]],
    *,
    minimum_total_unique_fraction: float = 0.15,
    maximum_single_species_signal_share: float = 0.15,
    minimum_species: int = 90,
) -> dict[str, object]:
    if len(by_species) < int(minimum_species):
        return {
            "status": "NOT_EVALUABLE_HISTORICAL_HOST_CONNECTIVITY_PREDICTOR_REDUNDANT",
            "species": len(by_species),
            "reason": "too_few_predictor_admissible_species",
        }
    items = []
    for species in sorted(by_species):
        lgm, present = by_species[species]
        q = historical_unique_fraction(lgm, present)
        items.append((species, len(np.asarray(lgm)), q))
    total_signal = float(sum(q for _, _, q in items))
    shares = [q / total_signal if total_signal > 0 else float("inf") for _, _, q in items]
    total_unique = float(np.mean([q for _, _, q in items]))
    max_share = float(max(shares))
    passed = (
        total_unique >= float(minimum_total_unique_fraction)
        and max_share <= float(maximum_single_species_signal_share)
    )
    return {
        "status": (
            "PASS_TO_HISTORICAL_HOST_SYNTHETIC_QUALIFICATION"
            if passed
            else "NOT_EVALUABLE_HISTORICAL_HOST_CONNECTIVITY_PREDICTOR_REDUNDANT"
        ),
        "species": len(items),
        "total_unique_fraction": total_unique,
        "max_single_species_signal_share": max_share,
        "species_information": [
            {
                "species": species,
                "edges": edges,
                "unique_fraction": q,
                "signal_share": share,
            }
            for (species, edges, q), share in zip(items, shares)
        ],
    }


__all__ = [
    "RESISTANCE_EPSILON",
    "edge_resistance",
    "historical_unique_fraction",
    "host_kernel_affinity",
    "host_union",
    "information_gate",
    "resistance_from_union",
    "zscore_population",
]
