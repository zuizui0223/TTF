from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import knn_edges
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test
from .mismatch import PairedSpeciesSample, pointwise_mismatch, pointwise_relative_state
from .mismatch_inference import PairedHeldoutInferenceResult, paired_heldout_species_bootstrap_test


@dataclass(frozen=True)
class OrientedPairedSpeciesSample:
    """Paired state sample plus a predeclared signed orientation coordinate.

    ``orientation`` must be fixed independently of observed mismatch outcomes.
    It may represent signed distance, increasing exposure, isolation, or another
    scientifically defensible prospectively declared direction axis.
    """

    sample: PairedSpeciesSample
    orientation: np.ndarray

    def __post_init__(self) -> None:
        g = np.asarray(self.orientation, dtype=float)
        if g.shape != (len(self.sample.coordinates),):
            raise ValueError("orientation must have one value per observation")
        if not np.isfinite(g).all():
            raise ValueError("orientation must be finite")
        if np.ptp(g) <= np.finfo(float).eps:
            raise ValueError("orientation must contain a signed range")
        object.__setattr__(self, "orientation", g)

    @property
    def species(self) -> str:
        return self.sample.species


@dataclass(frozen=True)
class DirectionalCouplingResult:
    """Direction-free TTF-C result plus held-out oriented mismatch inference."""

    coupling: PairedHeldoutInferenceResult
    front_center: float
    species_delta: Mapping[str, float]
    breakdown_bootstrap: MeanBootstrapResult
    recoupling_bootstrap: MeanBootstrapResult
    positive_fraction: float
    negative_fraction: float
    relation_detected: bool
    breakdown_label: bool
    recoupling_label: bool


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    if x.ndim != 1 or w.shape != x.shape or len(x) == 0:
        raise ValueError("weighted median inputs must be equal non-empty vectors")
    if not np.isfinite(x).all() or not np.isfinite(w).all() or np.any(w < 0):
        raise ValueError("weighted median inputs must be finite with non-negative weights")
    if float(w.sum()) <= np.finfo(float).tiny:
        return float(np.median(x))
    order = np.argsort(x, kind="stable")
    xs = x[order]
    ws = w[order]
    cutoff = 0.5 * float(ws.sum())
    return float(xs[np.searchsorted(np.cumsum(ws), cutoff, side="left")])


def _raw_relation_edge_dissimilarity(sample: PairedSpeciesSample, nodes: np.ndarray) -> np.ndarray:
    """Raw edge distance in the declared relation state, before TTF rank scaling."""

    relation = np.asarray(pointwise_relative_state(sample), dtype=float)
    left = relation[nodes[:, 0]]
    right = relation[nodes[:, 1]]
    delta = left - right
    if delta.ndim == 1:
        raw = np.abs(delta)
    else:
        raw = np.linalg.norm(delta.reshape(len(delta), -1), axis=1)
    if not np.isfinite(raw).all():
        raise ValueError("non-finite raw relation edge dissimilarity")
    return np.asarray(raw, dtype=float)


def estimate_training_front_center(
    samples: Sequence[OrientedPairedSpeciesSample],
    *,
    train_species: Sequence[str],
    k: int,
    front_edge_fraction: float = 0.15,
) -> float:
    """Estimate one orientation coordinate for the transferable relation front.

    Localization uses training-system relation states only. Raw relation edge
    dissimilarity is used here deliberately rather than rank-standardized TTF-C
    turnover: tied zero-change edges must not acquire artificial localization
    weight. The inferential TTF-C statistic itself remains unchanged.
    """

    if not (0.0 < float(front_edge_fraction) <= 0.5):
        raise ValueError("front_edge_fraction must be in (0, 0.5]")
    sample_map = {x.species: x for x in samples}
    centers: list[float] = []
    for name in map(str, train_species):
        oriented = sample_map[name]
        nodes = knn_edges(oriented.sample.coordinates, k=int(k))
        raw = _raw_relation_edge_dissimilarity(oriented.sample, nodes)
        gmid = 0.5 * (
            oriented.orientation[nodes[:, 0]] + oriented.orientation[nodes[:, 1]]
        )
        n_top = max(3, int(np.ceil(float(front_edge_fraction) * len(nodes))))
        order = np.argsort(raw, kind="stable")[-n_top:]
        centers.append(_weighted_median(gmid[order], raw[order]))
    if len(centers) < 2:
        raise ValueError("at least two training systems required for front localization")
    return float(np.median(np.asarray(centers, dtype=float)))


def heldout_oriented_mismatch_deltas(
    samples: Sequence[OrientedPairedSpeciesSample],
    *,
    eval_species: Sequence[str],
    front_center: float,
    front_window: float,
    min_side_records: int = 4,
) -> dict[str, float]:
    """Compute positive-minus-negative mismatch around a frozen front center."""

    if front_window <= 0 or min_side_records < 2:
        raise ValueError("front_window must be positive and min_side_records >= 2")
    sample_map = {x.species: x for x in samples}
    out: dict[str, float] = {}
    c = float(front_center)
    w = float(front_window)
    for name in map(str, eval_species):
        oriented = sample_map[name]
        g = oriented.orientation
        mismatch = pointwise_mismatch(oriented.sample)
        neg = (g >= c - w) & (g < c)
        pos = (g > c) & (g <= c + w)
        if int(np.count_nonzero(neg)) < int(min_side_records):
            continue
        if int(np.count_nonzero(pos)) < int(min_side_records):
            continue
        out[name] = float(mismatch[pos].mean() - mismatch[neg].mean())
    return out


def directional_coupling_test(
    samples: Sequence[OrientedPairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    k: int,
    bandwidth: float,
    front_edge_fraction: float,
    front_window: float,
    sign_consistency_fraction: float = 0.75,
    alpha: float = 0.05,
    n_bootstrap: int = 1999,
    seed: int = 0,
) -> DirectionalCouplingResult:
    """Qualify directional breakdown/recoupling after direction-free TTF-C.

    A directional label requires BOTH a significant held-out TTF-C relation
    front and a significant held-out oriented mismatch contrast with at least
    the predeclared fraction of species sharing the same sign.
    """

    if not (0.5 < float(sign_consistency_fraction) <= 1.0):
        raise ValueError("sign_consistency_fraction must be in (0.5, 1]")
    if not (0.0 < float(alpha) < 0.5):
        raise ValueError("alpha must be in (0, 0.5)")

    paired = [x.sample for x in samples]
    coupling = paired_heldout_species_bootstrap_test(
        paired,
        train_species=train_species,
        eval_species=eval_species,
        k=int(k),
        bandwidth=float(bandwidth),
        n_bootstrap=int(n_bootstrap),
        seed=int(seed),
        edge_chunk_size=32,
        train_chunk_size=2048,
    )
    center = estimate_training_front_center(
        samples,
        train_species=train_species,
        k=int(k),
        front_edge_fraction=float(front_edge_fraction),
    )
    deltas = heldout_oriented_mismatch_deltas(
        samples,
        eval_species=eval_species,
        front_center=center,
        front_window=float(front_window),
    )
    vec = np.asarray([deltas[s] for s in map(str, eval_species) if s in deltas], dtype=float)
    if len(vec) < 6:
        raise ValueError("fewer than six held-out systems have usable oriented contrasts")

    breakdown = centered_species_bootstrap_mean_test(
        vec,
        n_bootstrap=int(n_bootstrap),
        seed=int(seed) + 10_000,
    )
    recoupling = centered_species_bootstrap_mean_test(
        -vec,
        n_bootstrap=int(n_bootstrap),
        seed=int(seed) + 20_000,
    )
    positive_fraction = float(np.mean(vec > 0.0))
    negative_fraction = float(np.mean(vec < 0.0))
    relation_detected = bool(coupling.coupling_bootstrap.p_value < float(alpha))
    breakdown_label = bool(
        relation_detected
        and breakdown.p_value < float(alpha)
        and positive_fraction >= float(sign_consistency_fraction)
    )
    recoupling_label = bool(
        relation_detected
        and recoupling.p_value < float(alpha)
        and negative_fraction >= float(sign_consistency_fraction)
    )
    if breakdown_label and recoupling_label:
        raise RuntimeError("breakdown and recoupling labels cannot both be true")

    return DirectionalCouplingResult(
        coupling=coupling,
        front_center=center,
        species_delta=deltas,
        breakdown_bootstrap=breakdown,
        recoupling_bootstrap=recoupling,
        positive_fraction=positive_fraction,
        negative_fraction=negative_fraction,
        relation_detected=relation_detected,
        breakdown_label=breakdown_label,
        recoupling_label=recoupling_label,
    )


__all__ = [
    "OrientedPairedSpeciesSample",
    "DirectionalCouplingResult",
    "estimate_training_front_center",
    "heldout_oriented_mismatch_deltas",
    "directional_coupling_test",
]
