from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class TwoWayFEAbsorber:
    source: np.ndarray
    target: np.ndarray
    n_source: int
    n_target: int
    target_reference: int
    cholesky: np.ndarray

    def residualize(self, values: np.ndarray) -> np.ndarray:
        a = np.asarray(values, dtype=float)
        one_dimensional = a.ndim == 1
        if one_dimensional:
            a = a[:, None]
        if a.ndim != 2 or a.shape[0] != len(self.source):
            raise ValueError("values must have one row per dyad")
        if not np.isfinite(a).all():
            raise ValueError("values must be finite")

        p = a.shape[1]
        source_rhs = np.zeros((self.n_source, p), dtype=float)
        np.add.at(source_rhs, self.source, a)
        keep = self.target != self.target_reference
        target_rhs = np.zeros((self.n_target - 1, p), dtype=float)
        np.add.at(target_rhs, self.target[keep], a[keep])
        rhs = np.vstack((source_rhs, target_rhs))

        y = np.linalg.solve(self.cholesky, rhs)
        coef = np.linalg.solve(self.cholesky.T, y)
        source_coef = coef[: self.n_source]
        target_coef = coef[self.n_source :]

        fitted = source_coef[self.source].copy()
        fitted[keep] += target_coef[self.target[keep]]
        out = a - fitted
        return out[:, 0] if one_dimensional else out


@dataclass(frozen=True)
class PreparedDyadicRegression:
    absorber: TwoWayFEAbsorber
    x_residual: np.ndarray
    bread: np.ndarray
    condition_number: float
    primary_index: int
    primary_weight: np.ndarray


@dataclass(frozen=True)
class BatchPrimaryResult:
    coefficient: np.ndarray
    standard_error: np.ndarray
    z_score: np.ndarray
    p_value_one_sided: np.ndarray


def prepare_two_way_absorber(source: np.ndarray, target: np.ndarray) -> TwoWayFEAbsorber:
    s = np.asarray(source, dtype=np.int64)
    t = np.asarray(target, dtype=np.int64)
    if s.ndim != 1 or t.ndim != 1 or len(s) != len(t) or len(s) == 0:
        raise ValueError("source and target must be equal-length non-empty vectors")
    if np.any(s < 0) or np.any(t < 0):
        raise ValueError("cluster labels must be non-negative")
    n_source = int(s.max()) + 1
    n_target = int(t.max()) + 1
    if set(map(int, np.unique(s))) != set(range(n_source)):
        raise ValueError("source labels must be contiguous from zero")
    if set(map(int, np.unique(t))) != set(range(n_target)):
        raise ValueError("target labels must be contiguous from zero")
    if n_source < 2 or n_target < 2:
        raise ValueError("two-way clustering requires at least two source and target clusters")

    reference = n_target - 1
    source_counts = np.bincount(s, minlength=n_source).astype(float)
    target_counts = np.bincount(t, minlength=n_target).astype(float)
    cross = np.zeros((n_source, n_target - 1), dtype=float)
    keep = t != reference
    np.add.at(cross, (s[keep], t[keep]), 1.0)

    normal = np.zeros((n_source + n_target - 1, n_source + n_target - 1), dtype=float)
    normal[:n_source, :n_source] = np.diag(source_counts)
    normal[n_source:, n_source:] = np.diag(target_counts[:-1])
    normal[:n_source, n_source:] = cross
    normal[n_source:, :n_source] = cross.T
    try:
        chol = np.linalg.cholesky(normal)
    except np.linalg.LinAlgError as exc:
        raise ValueError("source-target support graph is not connected enough for exact two-way FE absorption") from exc
    return TwoWayFEAbsorber(
        source=s,
        target=t,
        n_source=n_source,
        n_target=n_target,
        target_reference=reference,
        cholesky=chol,
    )


def prepare_dyadic_regression(
    source: np.ndarray,
    target: np.ndarray,
    predictors: np.ndarray,
    *,
    primary_index: int = 0,
) -> PreparedDyadicRegression:
    x = np.asarray(predictors, dtype=float)
    if x.ndim != 2 or x.shape[1] < 1:
        raise ValueError("predictors must be n x p")
    if not 0 <= int(primary_index) < x.shape[1]:
        raise ValueError("invalid primary_index")
    absorber = prepare_two_way_absorber(source, target)
    xr = absorber.residualize(x)
    xtx = xr.T @ xr
    condition = float(np.linalg.cond(xtx))
    if not np.isfinite(condition):
        raise ValueError("non-finite predictor condition number")
    bread = np.linalg.inv(xtx)
    primary_weight = xr @ bread[:, int(primary_index)]
    return PreparedDyadicRegression(
        absorber=absorber,
        x_residual=xr,
        bread=bread,
        condition_number=condition,
        primary_index=int(primary_index),
        primary_weight=primary_weight,
    )


def _normal_upper_tail(z: np.ndarray) -> np.ndarray:
    zz = np.asarray(z, dtype=float)
    return np.asarray([0.5 * math.erfc(float(v) / math.sqrt(2.0)) for v in zz], dtype=float)


def batch_primary_test(
    prepared: PreparedDyadicRegression,
    responses: np.ndarray,
) -> BatchPrimaryResult:
    y = np.asarray(responses, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    if y.ndim != 2 or y.shape[0] != len(prepared.absorber.source):
        raise ValueError("responses must have one row per dyad")
    if not np.isfinite(y).all():
        raise ValueError("responses must be finite")

    yr = prepared.absorber.residualize(y)
    beta = prepared.bread @ (prepared.x_residual.T @ yr)
    residual = yr - prepared.x_residual @ beta

    h = prepared.primary_weight[:, None] * residual
    source = prepared.absorber.source
    target = prepared.absorber.target
    ns = prepared.absorber.n_source
    nt = prepared.absorber.n_target
    n, p = prepared.x_residual.shape

    source_score = np.zeros((ns, h.shape[1]), dtype=float)
    target_score = np.zeros((nt, h.shape[1]), dtype=float)
    np.add.at(source_score, source, h)
    np.add.at(target_score, target, h)

    source_correction = (ns / (ns - 1.0)) * ((n - 1.0) / (n - p))
    target_correction = (nt / (nt - 1.0)) * ((n - 1.0) / (n - p))
    dyad_correction = n / (n - p)
    variance = (
        source_correction * np.sum(source_score * source_score, axis=0)
        + target_correction * np.sum(target_score * target_score, axis=0)
        - dyad_correction * np.sum(h * h, axis=0)
    )
    variance = np.maximum(variance, np.finfo(float).tiny)
    standard_error = np.sqrt(variance)
    coefficient = beta[prepared.primary_index]
    z_score = coefficient / standard_error
    p_value = _normal_upper_tail(z_score)
    return BatchPrimaryResult(
        coefficient=np.asarray(coefficient, dtype=float),
        standard_error=np.asarray(standard_error, dtype=float),
        z_score=np.asarray(z_score, dtype=float),
        p_value_one_sided=p_value,
    )


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid binomial counts")
    p = successes / trials
    den = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / den
    half = z * math.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / den
    return center - half, center + half
