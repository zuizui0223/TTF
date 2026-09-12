from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .core import Dissimilarity, SpeciesEdges, SpeciesSample, build_species_edges

MismatchFunction = Callable[[np.ndarray, np.ndarray], float]
RelativeStateFunction = Callable[[np.ndarray, np.ndarray], np.ndarray]


@dataclass(frozen=True)
class PairedSpeciesSample:
    """Two co-located state fields observed on the same within-system support.

    TTF-M treats ``state_a`` and ``state_b`` as paired measurements at the
    same observation units.  The pairing is part of the estimand: rows may not
    be independently re-matched after outcomes are inspected.
    """

    species: str
    coordinates: np.ndarray
    state_a: np.ndarray
    state_b: np.ndarray
    blocks: np.ndarray | None = None

    def __post_init__(self) -> None:
        coords = np.asarray(self.coordinates, dtype=float)
        a = np.asarray(self.state_a)
        b = np.asarray(self.state_b)
        if coords.ndim != 2 or coords.shape[0] < 2:
            raise ValueError("coordinates must be n x d with n >= 2")
        if a.shape[0] != coords.shape[0] or b.shape[0] != coords.shape[0]:
            raise ValueError("state_a, state_b and coordinates must share the first dimension")
        if not np.isfinite(coords).all():
            raise ValueError("coordinates must be finite")
        blocks = None if self.blocks is None else np.asarray(self.blocks)
        if blocks is not None and blocks.shape != (coords.shape[0],):
            raise ValueError("blocks must have one value per observation")
        if not str(self.species):
            raise ValueError("species must be non-empty")
        object.__setattr__(self, "coordinates", coords)
        object.__setattr__(self, "state_a", a)
        object.__setattr__(self, "state_b", b)
        object.__setattr__(self, "blocks", blocks)


def default_mismatch(a: np.ndarray, b: np.ndarray) -> float:
    """Pointwise mismatch magnitude for commensurable numeric states.

    Scalars use absolute difference and vectors/tensors use Euclidean distance.
    Non-commensurable state spaces must supply an explicit mismatch function.
    """

    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    if aa.shape != bb.shape:
        raise ValueError("default_mismatch requires equal-shaped commensurable states")
    if aa.ndim == 0:
        return float(abs(float(aa) - float(bb)))
    return float(np.linalg.norm(aa - bb))


def default_relative_state(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Relative coupling state ``A - B`` for a declared common vector space."""

    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    if aa.shape != bb.shape:
        raise ValueError("default_relative_state requires equal-shaped states")
    return np.asarray(aa - bb, dtype=float)


def pointwise_mismatch(
    sample: PairedSpeciesSample,
    *,
    mismatch: MismatchFunction = default_mismatch,
) -> np.ndarray:
    """Compute the scalar mismatch state M_i = m(A_i, B_i)."""

    values = np.empty(len(sample.coordinates), dtype=float)
    for i in range(len(values)):
        values[i] = float(mismatch(sample.state_a[i], sample.state_b[i]))
    if not np.isfinite(values).all():
        raise ValueError("mismatch function produced non-finite values")
    return values


def pointwise_relative_state(
    sample: PairedSpeciesSample,
    *,
    relative_state: RelativeStateFunction = default_relative_state,
) -> np.ndarray:
    """Compute a declared relation state R_i = r(A_i, B_i).

    The default relation state is A-B and therefore requires A and B to live in
    the same numeric vector space.  A custom adapter may map heterogeneous
    paired states into any fixed numeric relation representation.
    """

    rows = [np.asarray(relative_state(sample.state_a[i], sample.state_b[i]), dtype=float) for i in range(len(sample.coordinates))]
    try:
        values = np.stack(rows, axis=0)
    except ValueError as exc:
        raise ValueError("relative_state must return equal-shaped numeric values") from exc
    if not np.isfinite(values).all():
        raise ValueError("relative_state produced non-finite values")
    return values


def mismatch_species_sample(
    sample: PairedSpeciesSample,
    *,
    mismatch: MismatchFunction = default_mismatch,
) -> SpeciesSample:
    """Project paired states to a scalar mismatch state for core TTF."""

    return SpeciesSample(
        species=sample.species,
        coordinates=sample.coordinates,
        trait=pointwise_mismatch(sample, mismatch=mismatch),
        blocks=sample.blocks,
    )


def relative_species_sample(
    sample: PairedSpeciesSample,
    *,
    relative_state: RelativeStateFunction = default_relative_state,
) -> SpeciesSample:
    """Project paired states to a relation/coupling state for core TTF."""

    return SpeciesSample(
        species=sample.species,
        coordinates=sample.coordinates,
        trait=pointwise_relative_state(sample, relative_state=relative_state),
        blocks=sample.blocks,
    )


def build_mismatch_edges(
    sample: PairedSpeciesSample,
    *,
    mismatch: MismatchFunction = default_mismatch,
    mismatch_dissimilarity: Dissimilarity | None = None,
    k: int = 4,
    max_distance: float | None = None,
    edge_nodes: np.ndarray | None = None,
) -> SpeciesEdges:
    """Build TTF edges for turnover of mismatch magnitude M(x).

    This asks where the *amount* of mismatch changes sharply.  Because the
    graph is undirected, the statistic identifies a mismatch-transition front,
    not which side represents increasing versus decreasing mismatch.
    """

    projected = mismatch_species_sample(sample, mismatch=mismatch)
    return build_species_edges(
        projected,
        k=k,
        max_distance=max_distance,
        dissimilarity=mismatch_dissimilarity,
        edge_nodes=edge_nodes,
    )


def build_coupling_edges(
    sample: PairedSpeciesSample,
    *,
    relative_state: RelativeStateFunction = default_relative_state,
    relation_dissimilarity: Dissimilarity | None = None,
    k: int = 4,
    max_distance: float | None = None,
    edge_nodes: np.ndarray | None = None,
) -> SpeciesEdges:
    """Build TTF edges for turnover of the A--B coupling relation R(x).

    With the default R=A-B, a common shift of A and B leaves R unchanged,
    whereas asymmetric shifts, opposite shifts, or rotations of their relative
    vector create coupling turnover.  This is a stricter relational estimand
    than mismatch magnitude alone: two sites can have equal ||A-B|| but
    different relative vectors.

    The output is direction-free coupling-change turnover.  Calling a boundary
    specifically a *breakdown* additionally requires a predeclared side-label
    or mismatch-level rule identifying the less-coupled side.
    """

    projected = relative_species_sample(sample, relative_state=relative_state)
    return build_species_edges(
        projected,
        k=k,
        max_distance=max_distance,
        dissimilarity=relation_dissimilarity,
        edge_nodes=edge_nodes,
    )
