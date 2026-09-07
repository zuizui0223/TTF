from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import SpeciesEdges, spearman_rho


@dataclass(frozen=True)
class PredictorTransferResult:
    statistic: float
    species_scores: Mapping[str, float]
    coefficients: np.ndarray
    intercept: float


@dataclass(frozen=True)
class PredictorCompetitionResult:
    scores: Mapping[str, float]
    increments: Mapping[str, float]


def _validated_features(
    edges: Sequence[SpeciesEdges],
    features: Mapping[str, np.ndarray],
    columns: Sequence[int],
) -> list[np.ndarray]:
    out: list[np.ndarray] = []
    for item in edges:
        matrix = np.asarray(features[item.species], dtype=float)
        if matrix.ndim != 2 or matrix.shape[0] != item.n_edges:
            raise ValueError(f"feature shape drift for {item.species}")
        selected = matrix[:, list(columns)]
        if not np.isfinite(selected).all():
            raise ValueError(f"non-finite predictor feature for {item.species}")
        out.append(selected)
    return out


def ridge_predictor_transfer(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    features: Mapping[str, np.ndarray],
    *,
    columns: Sequence[int],
    ridge: float = 1.0,
) -> PredictorTransferResult:
    """Evaluate a predictor space by species-disjoint edge-level transfer.

    Each training species receives equal total weight, preventing record-rich
    species from defining the predictor model.
    """
    if ridge < 0:
        raise ValueError("ridge must be non-negative")
    if not columns:
        raise ValueError("at least one predictor column is required")
    train_names = {e.species for e in train_edges}
    eval_names = {e.species for e in eval_edges}
    if train_names & eval_names:
        raise ValueError("train and evaluation species must be disjoint")

    train_x_list = _validated_features(train_edges, features, columns)
    eval_x_list = _validated_features(eval_edges, features, columns)
    x = np.vstack(train_x_list)
    y = np.concatenate([e.turnover for e in train_edges])
    w = np.concatenate(
        [np.full(e.n_edges, 1.0 / e.n_edges, dtype=float) for e in train_edges]
    )
    w /= w.sum()

    mean_x = np.sum(w[:, None] * x, axis=0)
    centered_x = x - mean_x
    scale_x = np.sqrt(np.sum(w[:, None] * centered_x * centered_x, axis=0))
    scale_x = np.where(scale_x > 1e-12, scale_x, 1.0)
    z = centered_x / scale_x
    mean_y = float(np.dot(w, y))
    centered_y = y - mean_y

    weighted_z = z * w[:, None]
    gram = z.T @ weighted_z + float(ridge) * np.eye(z.shape[1])
    rhs = z.T @ (w * centered_y)
    beta = np.linalg.solve(gram, rhs)

    scores: dict[str, float] = {}
    for edges, matrix in zip(eval_edges, eval_x_list):
        z_eval = (matrix - mean_x) / scale_x
        predicted = mean_y + z_eval @ beta
        scores[edges.species] = spearman_rho(predicted, edges.turnover)
    finite = np.asarray([v for v in scores.values() if np.isfinite(v)], dtype=float)
    if len(finite) == 0:
        raise ValueError("no finite evaluation species scores")
    return PredictorTransferResult(
        statistic=float(finite.mean()),
        species_scores=scores,
        coefficients=beta / scale_x,
        intercept=float(mean_y - np.dot(mean_x / scale_x, beta)),
    )


def compete_predictor_spaces(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    features: Mapping[str, np.ndarray],
    *,
    spaces: Mapping[str, Sequence[int]],
    increments: Mapping[str, tuple[str, str]] | None = None,
    ridge: float = 1.0,
) -> PredictorCompetitionResult:
    """Compare D/G/E/R feature spaces on identical held-out species.

    ``increments`` maps a contrast name to ``(larger_space, baseline_space)``.
    For example ``{"E|G": ("GE", "G")}`` returns T_GE - T_G.
    """
    scores: dict[str, float] = {}
    for name, columns in spaces.items():
        result = ridge_predictor_transfer(
            train_edges,
            eval_edges,
            features,
            columns=columns,
            ridge=ridge,
        )
        scores[str(name)] = result.statistic

    delta: dict[str, float] = {}
    for name, pair in (increments or {}).items():
        larger, baseline = pair
        if larger not in scores or baseline not in scores:
            raise ValueError(f"unknown predictor space in increment {name}")
        delta[str(name)] = float(scores[larger] - scores[baseline])
    return PredictorCompetitionResult(scores=scores, increments=delta)
