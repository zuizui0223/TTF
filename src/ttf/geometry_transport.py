from __future__ import annotations

from typing import Literal, Sequence

import numpy as np

from .core import SpeciesEdges

GeometryTransportMode = Literal[
    "pooled_inverse_density",
    "self_inverse_density",
    "target_density_ratio",
]


def _stack_species_equal(
    edge_sets: Sequence[SpeciesEdges],
) -> tuple[np.ndarray, np.ndarray, dict[str, slice]]:
    data = tuple(edge_sets)
    if not data:
        raise ValueError("at least one edge set is required")
    names = [str(edges.species) for edges in data]
    if len(set(names)) != len(names):
        raise ValueError("system labels must be unique")
    positions: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    slices: dict[str, slice] = {}
    cursor = 0
    dimensionality: int | None = None
    for edges in data:
        midpoint = np.asarray(edges.midpoint, dtype=float)
        if midpoint.ndim != 2 or len(midpoint) < 1 or not np.isfinite(midpoint).all():
            raise ValueError("edge midpoints must be finite non-empty matrices")
        if dimensionality is None:
            dimensionality = int(midpoint.shape[1])
        elif midpoint.shape[1] != dimensionality:
            raise ValueError("edge midpoint dimensionality mismatch")
        n = int(len(midpoint))
        slices[str(edges.species)] = slice(cursor, cursor + n)
        cursor += n
        positions.append(midpoint)
        weights.append(np.full(n, 1.0 / n, dtype=float))
    return np.vstack(positions), np.concatenate(weights), slices


def _gaussian_density(
    query: np.ndarray,
    support: np.ndarray,
    support_weights: np.ndarray,
    *,
    bandwidth: float,
    chunk_size: int = 256,
) -> np.ndarray:
    q = np.asarray(query, dtype=float)
    s = np.asarray(support, dtype=float)
    w = np.asarray(support_weights, dtype=float)
    if q.ndim != 2 or s.ndim != 2 or q.shape[1] != s.shape[1]:
        raise ValueError("query/support dimensionality mismatch")
    if w.shape != (len(s),) or np.any(w <= 0) or not np.isfinite(w).all():
        raise ValueError("support weights must be finite and positive")
    if not np.isfinite(q).all() or not np.isfinite(s).all():
        raise ValueError("geometry must be finite")
    if bandwidth <= 0 or chunk_size < 1:
        raise ValueError("bandwidth and chunk_size must be positive")

    out = np.empty(len(q), dtype=float)
    h2 = float(bandwidth) * float(bandwidth)
    for start in range(0, len(q), int(chunk_size)):
        stop = min(start + int(chunk_size), len(q))
        delta = q[start:stop, None, :] - s[None, :, :]
        distance2 = np.sum(delta * delta, axis=2)
        kernel = np.exp(-0.5 * distance2 / h2)
        out[start:stop] = kernel @ w
    if not np.isfinite(out).all() or np.any(out < 0):
        raise RuntimeError("geometry density became invalid")
    return out


def _normalize_species_mass(
    raw: np.ndarray,
    slices: dict[str, slice],
) -> np.ndarray:
    values = np.asarray(raw, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("raw transport weights must be finite and positive")
    out = np.empty_like(values)
    for sl in slices.values():
        block = values[sl]
        total = float(block.sum())
        if total <= np.finfo(float).tiny:
            raise RuntimeError("transport weight normalization underflow")
        out[sl] = block / total
    return out


def geometry_transport_train_weights(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    *,
    mode: GeometryTransportMode,
    bandwidth: float,
    density_chunk_size: int = 256,
) -> np.ndarray:
    """Return response-blind training-edge masses for a Q5.2 transport candidate.

    Every training system retains total mass one, exactly as core TTF. Only the
    within-system distribution of that mass changes. No turnover or state value
    is accepted by this function, which makes outcome leakage impossible at the
    API level.

    ``pooled_inverse_density`` downweights training-edge locations that are
    overrepresented in the species-equal pooled training geometry.

    ``self_inverse_density`` equalizes each training system against its own
    midpoint sampling density.

    ``target_density_ratio`` is the more invasive transductive candidate. It
    aligns training-edge mass to the species-equal held-out midpoint geometry
    using only held-out coordinates/graphs, never held-out responses.
    """
    allowed = {
        "pooled_inverse_density",
        "self_inverse_density",
        "target_density_ratio",
    }
    if mode not in allowed:
        raise ValueError(f"unknown geometry transport mode: {mode}")
    train = tuple(train_edges)
    evaluation = tuple(eval_edges)
    train_pos, train_base, train_slices = _stack_species_equal(train)
    eval_pos, eval_base, _ = _stack_species_equal(evaluation)
    if train_pos.shape[1] != eval_pos.shape[1]:
        raise ValueError("train/evaluation geometry dimensionality mismatch")
    tiny = np.finfo(float).tiny

    if mode == "self_inverse_density":
        raw = np.empty(len(train_pos), dtype=float)
        for edges in train:
            sl = train_slices[str(edges.species)]
            midpoint = train_pos[sl]
            base = np.full(len(midpoint), 1.0 / len(midpoint), dtype=float)
            density = _gaussian_density(
                midpoint,
                midpoint,
                base,
                bandwidth=float(bandwidth),
                chunk_size=int(density_chunk_size),
            )
            raw[sl] = 1.0 / np.maximum(density, tiny)
        return _normalize_species_mass(raw, train_slices)

    train_density = _gaussian_density(
        train_pos,
        train_pos,
        train_base,
        bandwidth=float(bandwidth),
        chunk_size=int(density_chunk_size),
    )
    if mode == "pooled_inverse_density":
        raw = 1.0 / np.maximum(train_density, tiny)
        return _normalize_species_mass(raw, train_slices)

    eval_density = _gaussian_density(
        train_pos,
        eval_pos,
        eval_base,
        bandwidth=float(bandwidth),
        chunk_size=int(density_chunk_size),
    )
    raw = np.maximum(eval_density, tiny) / np.maximum(train_density, tiny)
    return _normalize_species_mass(raw, train_slices)


def transport_effective_edge_count(
    weights: np.ndarray,
    train_edges: Sequence[SpeciesEdges],
) -> dict[str, float]:
    """Kish-style effective edge count for audit; never used for selection."""
    _, _, slices = _stack_species_equal(tuple(train_edges))
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1:
        raise ValueError("weights must be one-dimensional")
    expected = max(sl.stop for sl in slices.values())
    if len(w) != expected:
        raise ValueError("weight length mismatch")
    out: dict[str, float] = {}
    for name, sl in slices.items():
        block = w[sl]
        out[name] = float(1.0 / np.sum(block * block))
    return out


__all__ = [
    "GeometryTransportMode",
    "geometry_transport_train_weights",
    "transport_effective_edge_count",
]
