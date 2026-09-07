from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import SpeciesEdges, average_ranks
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test
from .transfer import PreparedTransfer, prepare_transfer


@dataclass(frozen=True)
class RankOrthogonalization:
    residual: np.ndarray
    fitted: np.ndarray
    nuisance_rank_r2: float
    n_active_covariates: int


@dataclass(frozen=True)
class OrthogonalTransferResult:
    """Held-out rank transfer after removing geometry-only predictor components."""

    statistic: float
    species_scores: Mapping[str, float]
    nuisance_rank_r2: Mapping[str, float]
    n_eval_species: int


@dataclass(frozen=True)
class OrthogonalBootstrapResult:
    observed: OrthogonalTransferResult
    species_scores: np.ndarray
    bootstrap: MeanBootstrapResult

    @property
    def p_value(self) -> float:
        return float(self.bootstrap.p_value)


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.shape != yy.shape or xx.ndim != 1 or len(xx) < 3:
        raise ValueError("x and y must be equal-length vectors with n >= 3")
    dx = xx - xx.mean()
    dy = yy - yy.mean()
    nx = float(np.linalg.norm(dx))
    ny = float(np.linalg.norm(dy))
    scale_x = max(1.0, float(np.linalg.norm(xx)))
    scale_y = max(1.0, float(np.linalg.norm(yy)))
    # Rank-space least squares can leave residuals at ~1e-14 when the
    # predictor is exactly explained by geometry.  Treat those as constant
    # rather than turning round-off into an arbitrary correlation.
    if nx <= 1e-12 * scale_x or ny <= 1e-12 * scale_y:
        return 0.0
    return float(np.dot(dx, dy) / (nx * ny))


def orthogonalize_rank_predictor(
    predictor: np.ndarray,
    geometry_covariates: np.ndarray,
) -> RankOrthogonalization:
    """Remove geometry-only components from a predictor in rank space.

    The outcome/turnover vector is deliberately absent from this function.  Each
    predictor and nuisance covariate is converted to average ranks.  OLS then
    projects the ranked predictor onto an intercept plus non-constant ranked
    geometry covariates.  The returned residual is therefore exactly orthogonal
    (up to numerical precision) to the active nuisance columns in rank space.
    """
    y = np.asarray(predictor, dtype=float)
    z = np.asarray(geometry_covariates, dtype=float)
    if y.ndim != 1 or len(y) < 3 or not np.isfinite(y).all():
        raise ValueError("predictor must be a finite vector with n >= 3")
    if z.ndim == 1:
        z = z[:, None]
    if z.ndim != 2 or z.shape[0] != len(y) or not np.isfinite(z).all():
        raise ValueError("geometry_covariates must be finite n x p")

    ranked_y = average_ranks(y)
    active: list[np.ndarray] = []
    for j in range(z.shape[1]):
        ranked = average_ranks(z[:, j])
        centered = ranked - ranked.mean()
        if float(np.dot(centered, centered)) > np.finfo(float).tiny:
            active.append(centered)

    columns = [np.ones(len(y), dtype=float), *active]
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, ranked_y, rcond=None)
    fitted = design @ beta
    residual = ranked_y - fitted

    centered_y = ranked_y - ranked_y.mean()
    ss_total = float(np.dot(centered_y, centered_y))
    ss_resid = float(np.dot(residual, residual))
    r2 = 0.0 if ss_total <= np.finfo(float).tiny else 1.0 - ss_resid / ss_total
    return RankOrthogonalization(
        residual=np.asarray(residual, dtype=float),
        fitted=np.asarray(fitted, dtype=float),
        nuisance_rank_r2=float(np.clip(r2, 0.0, 1.0)),
        n_active_covariates=int(len(active)),
    )


def semi_partial_rank_score(
    predictor: np.ndarray,
    target: np.ndarray,
    geometry_covariates: np.ndarray,
) -> tuple[float, RankOrthogonalization]:
    """Correlation of geometry-orthogonalized predictor ranks with target ranks."""
    target = np.asarray(target, dtype=float)
    predictor = np.asarray(predictor, dtype=float)
    if target.shape != predictor.shape or target.ndim != 1:
        raise ValueError("predictor and target must be equal-length vectors")
    if not np.isfinite(target).all():
        raise ValueError("target must be finite")
    orth = orthogonalize_rank_predictor(predictor, geometry_covariates)
    target_rank = average_ranks(target)
    return _pearson(orth.residual, target_rank), orth


def _train_vector(
    prepared: PreparedTransfer,
    train_turnover: Mapping[str, np.ndarray],
) -> np.ndarray:
    out = np.empty(max(sl.stop for sl in prepared.train_slices.values()), dtype=float)
    for species in prepared.train_species:
        sl = prepared.train_slices[species]
        values = np.asarray(train_turnover[species], dtype=float)
        if values.shape != (sl.stop - sl.start,):
            raise ValueError(f"train turnover shape drift for {species}")
        out[sl] = values
    return out


def score_orthogonalized_transfer(
    prepared: PreparedTransfer,
    *,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    eval_edge_length: Mapping[str, np.ndarray] | None = None,
    include_edge_length: bool = False,
) -> OrthogonalTransferResult:
    """Score held-out transfer after predictor-side geometry orthogonalization.

    Nuisance covariates contain no evaluation turnover information.  The base
    covariate is training-field opportunity/support.  Optionally, held-out edge
    length is added as a second geometry-only nuisance.  The target is used only
    in the final rank correlation.
    """
    if include_edge_length and eval_edge_length is None:
        raise ValueError("eval_edge_length is required when include_edge_length=True")
    train_values = _train_vector(prepared, train_turnover)
    scores: dict[str, float] = {}
    r2: dict[str, float] = {}
    for species in prepared.eval_species:
        target = np.asarray(eval_turnover[species], dtype=float)
        predicted = prepared.eval_projection[species] @ train_values + prepared.eval_prior_offset[species]
        opportunity = np.asarray(prepared.eval_opportunity[species], dtype=float)
        covars = [opportunity]
        if include_edge_length:
            length = np.asarray(eval_edge_length[species], dtype=float)
            if length.shape != predicted.shape:
                raise ValueError(f"evaluation edge-length shape drift for {species}")
            covars.append(length)
        geometry = np.column_stack(covars)
        if target.shape != predicted.shape:
            raise ValueError(f"evaluation turnover shape drift for {species}")
        score, orth = semi_partial_rank_score(predicted, target, geometry)
        scores[species] = float(score)
        r2[species] = float(orth.nuisance_rank_r2)

    finite = np.asarray([scores[s] for s in prepared.eval_species if np.isfinite(scores[s])], dtype=float)
    if len(finite) == 0:
        raise ValueError("no finite orthogonalized held-out species scores")
    return OrthogonalTransferResult(
        statistic=float(finite.mean()),
        species_scores=scores,
        nuisance_rank_r2=r2,
        n_eval_species=int(len(finite)),
    )


def orthogonalized_transfer_test(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    *,
    bandwidth: float,
    include_edge_length: bool = False,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
    n_bootstrap: int = 1999,
    seed: int = 0,
) -> OrthogonalBootstrapResult:
    """Held-out species bootstrap for geometry-orthogonalized transfer."""
    if len(eval_edges) < 6:
        raise ValueError("orthogonalized inference requires at least six held-out species")
    prepared = prepare_transfer(
        train_edges,
        eval_edges,
        bandwidth=bandwidth,
        prior_strength=prior_strength,
        prior_mean=prior_mean,
        segment_points=segment_points,
    )
    train_turnover = {edges.species: edges.turnover for edges in train_edges}
    eval_turnover = {edges.species: edges.turnover for edges in eval_edges}
    eval_length = {edges.species: edges.length for edges in eval_edges}
    observed = score_orthogonalized_transfer(
        prepared,
        train_turnover=train_turnover,
        eval_turnover=eval_turnover,
        eval_edge_length=eval_length,
        include_edge_length=include_edge_length,
    )
    scores = np.asarray(
        [observed.species_scores[name] for name in prepared.eval_species if np.isfinite(observed.species_scores[name])],
        dtype=float,
    )
    if len(scores) < 6:
        raise ValueError("fewer than six finite orthogonalized species scores")
    bootstrap = centered_species_bootstrap_mean_test(
        scores,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    if not np.isclose(observed.statistic, bootstrap.observed_mean, atol=1e-12, rtol=0.0):
        raise RuntimeError("orthogonalized macro-average drift")
    return OrthogonalBootstrapResult(
        observed=observed,
        species_scores=scores,
        bootstrap=bootstrap,
    )
