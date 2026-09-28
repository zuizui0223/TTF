from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .relational_dyadic import prepare_two_way_absorber


@dataclass(frozen=True)
class KnownTruthDyadicDesign:
    source: np.ndarray
    target: np.ndarray
    primary: np.ndarray
    baseline_control: np.ndarray
    geographic_control: np.ndarray
    unique_component: np.ndarray
    endpoint_component: np.ndarray
    truth: dict[str, float | str | int]


def complete_bipartite(n_source: int, n_target: int) -> tuple[np.ndarray, np.ndarray]:
    if n_source < 3 or n_target < 3:
        raise ValueError("known-truth benchmark requires at least three source and target nodes")
    source = np.repeat(np.arange(int(n_source), dtype=np.int64), int(n_target))
    target = np.tile(np.arange(int(n_target), dtype=np.int64), int(n_source))
    return source, target


def _unit_ss(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x - float(x.mean())
    ss = float(np.dot(x, x))
    if ss <= np.finfo(float).eps:
        raise ValueError("cannot normalize zero-variance component")
    return x * np.sqrt(len(x) / ss)


def _orthogonalize(values: np.ndarray, basis: list[np.ndarray]) -> np.ndarray:
    x = np.asarray(values, dtype=float).copy()
    for column in basis:
        den = float(np.dot(column, column))
        if den <= np.finfo(float).eps:
            raise ValueError("degenerate benchmark basis")
        x -= column * (float(np.dot(column, x)) / den)
    return x


def _fe_component(
    absorber,
    rng: np.random.Generator,
    *,
    basis: list[np.ndarray],
    scale_by_source: np.ndarray | None = None,
    scale_by_target: np.ndarray | None = None,
) -> np.ndarray:
    raw = rng.normal(size=len(absorber.source))
    if scale_by_source is not None:
        raw *= np.asarray(scale_by_source, dtype=float)[absorber.source]
    if scale_by_target is not None:
        raw *= np.asarray(scale_by_target, dtype=float)[absorber.target]
    residual = absorber.residualize(raw)
    residual = _orthogonalize(residual, basis)
    return _unit_ss(residual)


def known_truth_design(
    *,
    n_source: int = 20,
    n_target: int = 20,
    endpoint_retained_fraction: float = 0.50,
    baseline_survival_after_fe: float = 0.60,
    geography_survival_after_baseline: float = 0.80,
    concentration: str = "broad",
    seed: int = 1,
) -> KnownTruthDyadicDesign:
    """Construct a dyadic predictor with exact nested information survival.

    The focal relation is assembled from mutually orthogonal components in the
    exact source/target-FE residual subspace plus one endpoint-only component.

    Let f be the fraction surviving endpoint identity, b the fraction of that
    FE-residual information surviving the baseline nuisance control, and g the
    fraction of C1 information surviving addition of geography. Then the true
    cumulative retained fractions are exactly

        C1 = f * b
        C2 = f * b * g

    up to floating-point error.
    """
    f = float(endpoint_retained_fraction)
    b = float(baseline_survival_after_fe)
    g = float(geography_survival_after_baseline)
    if not (0 < f <= 1 and 0 < b <= 1 and 0 < g <= 1):
        raise ValueError("all survival fractions must lie in (0, 1]")
    if concentration not in {"broad", "source", "target"}:
        raise ValueError("concentration must be broad, source, or target")

    source, target = complete_bipartite(n_source, n_target)
    absorber = prepare_two_way_absorber(source, target)
    rng = np.random.default_rng(int(seed))

    source_effect = rng.normal(size=n_source)
    target_effect = rng.normal(size=n_target)
    endpoint = _unit_ss(source_effect[source] + target_effect[target])
    endpoint_residual = absorber.residualize(endpoint)
    if float(np.dot(endpoint_residual, endpoint_residual)) > 1e-18 * len(source):
        raise RuntimeError("endpoint-only benchmark component was not absorbed")

    baseline = _fe_component(absorber, rng, basis=[])
    geography = _fe_component(absorber, rng, basis=[baseline])

    source_scale = None
    target_scale = None
    if concentration == "source":
        source_scale = np.full(n_source, 0.05, dtype=float)
        source_scale[: max(1, n_source // 5)] = 20.0
    elif concentration == "target":
        target_scale = np.full(n_target, 0.05, dtype=float)
        target_scale[: max(1, n_target // 5)] = 20.0

    unique = _fe_component(
        absorber,
        rng,
        basis=[baseline, geography],
        scale_by_source=source_scale,
        scale_by_target=target_scale,
    )

    weight_baseline = 1.0 - b
    weight_unique = b * g
    weight_geography = b * (1.0 - g)
    fe_relation = (
        np.sqrt(weight_baseline) * baseline
        + np.sqrt(weight_geography) * geography
        + np.sqrt(weight_unique) * unique
    )
    primary = np.sqrt(f) * fe_relation + np.sqrt(1.0 - f) * endpoint

    truth = {
        "n_source": int(n_source),
        "n_target": int(n_target),
        "dyads": int(len(source)),
        "endpoint_retained_fraction": f,
        "baseline_survival_after_fe": b,
        "geography_survival_after_baseline": g,
        "C1_total_unique_variance_fraction": f * b,
        "C2_total_unique_variance_fraction": f * b * g,
        "baseline_attributable_fraction_after_fe": 1.0 - b,
        "geography_attributable_fraction_of_C1": 1.0 - g,
        "concentration": concentration,
        "seed": int(seed),
    }
    return KnownTruthDyadicDesign(
        source=source,
        target=target,
        primary=primary,
        baseline_control=baseline,
        geographic_control=geography,
        unique_component=unique,
        endpoint_component=endpoint,
        truth=truth,
    )
