from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Iterable, Sequence

import numpy as np

from .relational_dyadic import (
    batch_primary_test,
    prepare_dyadic_regression,
    prepare_two_way_absorber,
    wilson_interval,
)


@dataclass(frozen=True)
class RelationRepeatability:
    replicates: int
    dyads: int
    repeatability_icc: float
    between_dyad_variance: float
    within_dyad_variance: float
    median_within_dyad_sd: float


@dataclass(frozen=True)
class UniqueRelationInformation:
    dyads: int
    source_clusters: int
    target_clusters: int
    control_columns: int
    control_rank_after_two_way_fe: int
    source_target_fe_retained_variance_fraction: float
    control_unique_variance_fraction_after_fe: float
    total_unique_variance_fraction: float
    partial_sd_in_primary_sd_units: float
    vif_like_after_fe: float


@dataclass(frozen=True)
class DyadicSignalSupport:
    dyads: int
    source_clusters: int
    target_clusters: int
    signal_effective_dyads: float
    signal_effective_sources: float
    signal_effective_targets: float
    max_source_signal_share: float
    max_target_signal_share: float


def _as_finite_vector(values: np.ndarray | Sequence[float], *, name: str) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) == 0:
        raise ValueError(f"{name} must be a non-empty vector")
    if not np.isfinite(x).all():
        raise ValueError(f"{name} must be finite")
    return x


def _remap_labels(values: np.ndarray | Sequence[int]) -> np.ndarray:
    x = np.asarray(values)
    if x.ndim != 1 or len(x) == 0:
        raise ValueError("cluster labels must be a non-empty vector")
    unique = sorted(np.unique(x).tolist())
    lookup = {value: i for i, value in enumerate(unique)}
    return np.asarray([lookup[value] for value in x], dtype=np.int64)


def _zscore(values: np.ndarray | Sequence[float]) -> np.ndarray:
    x = _as_finite_vector(values, name="values")
    sd = float(x.std())
    if sd <= np.finfo(float).eps:
        raise ValueError("cannot standardize a constant vector")
    return (x - float(x.mean())) / sd


def relation_repeatability(replicates_by_dyad: np.ndarray) -> RelationRepeatability:
    """Estimate response-blind repeatability of a pairwise ecological relation.

    Input rows are independent resampling replicates and columns are fixed dyads.
    The returned ICC is the balanced one-way random-effects method-of-moments
    repeatability, clipped to [0, 1] because negative estimates represent no
    detectable between-dyad repeatability for qualification purposes.
    """
    y = np.asarray(replicates_by_dyad, dtype=float)
    if y.ndim != 2 or y.shape[0] < 2 or y.shape[1] < 2:
        raise ValueError("replicate relation matrix must be at least 2 x 2")
    if not np.isfinite(y).all():
        raise ValueError("replicate relation matrix must be finite")
    k, n = y.shape
    dyad_mean = y.mean(axis=0)
    grand = float(y.mean())
    ss_between = float(k * np.sum((dyad_mean - grand) ** 2))
    ss_within = float(np.sum((y - dyad_mean[None, :]) ** 2))
    ms_between = ss_between / (n - 1)
    ms_within = ss_within / (n * (k - 1))
    denom = ms_between + (k - 1) * ms_within
    raw_icc = 0.0 if denom <= np.finfo(float).eps else (ms_between - ms_within) / denom
    icc = float(np.clip(raw_icc, 0.0, 1.0))
    between_var = max((ms_between - ms_within) / k, 0.0)
    within_sd = np.std(y, axis=0, ddof=1)
    return RelationRepeatability(
        replicates=int(k),
        dyads=int(n),
        repeatability_icc=icc,
        between_dyad_variance=float(between_var),
        within_dyad_variance=float(ms_within),
        median_within_dyad_sd=float(np.median(within_sd)),
    )


def residualize_primary_relation(
    source: np.ndarray | Sequence[int],
    target: np.ndarray | Sequence[int],
    primary: np.ndarray | Sequence[float],
    controls: np.ndarray,
) -> tuple[np.ndarray, UniqueRelationInformation]:
    """Remove source/target identity and measured controls from a focal relation.

    `primary` should already be on the scale used by the intended model (usually
    z-scored). `controls` should contain the exact transformed nuisance columns
    used by that model. The residual is therefore the portion of the relation
    that can identify a distinct source-target effect under that specification.
    """
    s = _remap_labels(source)
    t = _remap_labels(target)
    p = _as_finite_vector(primary, name="primary")
    c = np.asarray(controls, dtype=float)
    if len(s) != len(t) or len(s) != len(p):
        raise ValueError("source, target, and primary must have equal length")
    if c.ndim == 1:
        c = c[:, None]
    if c.ndim != 2 or c.shape[0] != len(p):
        raise ValueError("controls must have one row per dyad")
    if not np.isfinite(c).all():
        raise ValueError("controls must be finite")

    absorber = prepare_two_way_absorber(s, t)
    p_fe = absorber.residualize(p)
    raw_ss = float(np.dot(p - p.mean(), p - p.mean()))
    fe_ss = float(np.dot(p_fe, p_fe))
    if raw_ss <= np.finfo(float).eps:
        raise ValueError("primary has no raw variation after centering")
    if fe_ss <= np.finfo(float).eps:
        residual = np.zeros_like(p_fe)
        rank = 0 if c.shape[1] == 0 else int(np.linalg.matrix_rank(absorber.residualize(c)))
    elif c.shape[1] == 0:
        residual = p_fe.copy()
        rank = 0
    else:
        c_fe = absorber.residualize(c)
        rank = int(np.linalg.matrix_rank(c_fe))
        coef, _, _, _ = np.linalg.lstsq(c_fe, p_fe, rcond=None)
        residual = p_fe - c_fe @ coef

    unique_ss = float(np.dot(residual, residual))
    fe_fraction = fe_ss / raw_ss
    conditional_fraction = 0.0 if fe_ss <= np.finfo(float).eps else unique_ss / fe_ss
    total_fraction = unique_ss / raw_ss
    partial_sd = float(np.sqrt(max(total_fraction, 0.0)))
    vif_like = float("inf") if conditional_fraction <= np.finfo(float).eps else 1.0 / conditional_fraction
    summary = UniqueRelationInformation(
        dyads=int(len(p)),
        source_clusters=int(len(np.unique(s))),
        target_clusters=int(len(np.unique(t))),
        control_columns=int(c.shape[1]),
        control_rank_after_two_way_fe=rank,
        source_target_fe_retained_variance_fraction=float(fe_fraction),
        control_unique_variance_fraction_after_fe=float(conditional_fraction),
        total_unique_variance_fraction=float(total_fraction),
        partial_sd_in_primary_sd_units=partial_sd,
        vif_like_after_fe=vif_like,
    )
    return residual, summary


def dyadic_signal_support(
    source: np.ndarray | Sequence[int],
    target: np.ndarray | Sequence[int],
    residual_primary: np.ndarray | Sequence[float],
) -> DyadicSignalSupport:
    """Describe where the identifiable focal-relation signal is concentrated.

    Effective counts are inverse-Herfindahl concentration measures of squared
    residual relation signal. They are diagnostics, not independent-sample-size
    estimates.
    """
    s = np.asarray(source)
    t = np.asarray(target)
    r = _as_finite_vector(residual_primary, name="residual_primary")
    if s.ndim != 1 or t.ndim != 1 or len(s) != len(t) or len(s) != len(r):
        raise ValueError("source, target, and residual_primary must align")
    signal = r * r
    total = float(signal.sum())
    if total <= np.finfo(float).eps:
        raise ValueError("residual relation carries no identifiable signal")

    dyad_share = signal / total
    source_values = np.unique(s)
    target_values = np.unique(t)
    source_share = np.asarray([signal[s == value].sum() / total for value in source_values], dtype=float)
    target_share = np.asarray([signal[t == value].sum() / total for value in target_values], dtype=float)

    def effective_number(shares: np.ndarray) -> float:
        return float(1.0 / np.sum(shares * shares))

    return DyadicSignalSupport(
        dyads=int(len(r)),
        source_clusters=int(len(source_values)),
        target_clusters=int(len(target_values)),
        signal_effective_dyads=effective_number(dyad_share),
        signal_effective_sources=effective_number(source_share),
        signal_effective_targets=effective_number(target_share),
        max_source_signal_share=float(source_share.max()),
        max_target_signal_share=float(target_share.max()),
    )


def _cell_seed(namespace: str, amplitude: float, effect: float) -> int:
    token = f"{namespace}|A={amplitude:.12g}|beta={effect:.12g}"
    return int(hashlib.sha256(token.encode()).hexdigest()[:16], 16)


def detectability_surface(
    source: np.ndarray | Sequence[int],
    target: np.ndarray | Sequence[int],
    predictors: np.ndarray,
    *,
    effects: Iterable[float],
    private_amplitudes: Iterable[float] = (0.0,),
    primary_index: int = 0,
    alpha: float = 0.05,
    worlds_per_cell: int = 1000,
    source_intercept_sd: float = 0.18,
    target_intercept_sd: float = 0.18,
    dyad_noise_sd: float = 0.32,
    fixed_coefficients: Sequence[float] | None = None,
    source_random_slope_index: int | None = None,
    target_random_slope_index: int | None = None,
    random_slope_sd_per_amplitude: float = 0.0,
    response_transform: str = "tanh",
    seed_namespace: str = "ttf-q-v0.1",
    block_size: int = 100,
) -> list[dict[str, float | int]]:
    """Map Type-I error and power over effect size x unmodelled-heterogeneity grids."""
    s = _remap_labels(source)
    t = _remap_labels(target)
    x = np.asarray(predictors, dtype=float)
    if x.ndim != 2 or x.shape[0] != len(s) or x.shape[1] < 1:
        raise ValueError("predictors must be n x p and align with dyads")
    if not np.isfinite(x).all():
        raise ValueError("predictors must be finite")
    if not 0 <= int(primary_index) < x.shape[1]:
        raise ValueError("invalid primary_index")
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0, 1)")
    if worlds_per_cell < 1 or block_size < 1:
        raise ValueError("worlds_per_cell and block_size must be positive")
    if response_transform not in {"tanh", "identity"}:
        raise ValueError("response_transform must be tanh or identity")

    prepared = prepare_dyadic_regression(s, t, x, primary_index=int(primary_index))
    p = x.shape[1]
    fixed = np.zeros(p, dtype=float) if fixed_coefficients is None else np.asarray(fixed_coefficients, dtype=float)
    if fixed.shape != (p,) or not np.isfinite(fixed).all():
        raise ValueError("fixed_coefficients must be a finite length-p vector")
    fixed = fixed.copy()
    fixed[int(primary_index)] = 0.0

    for idx in (source_random_slope_index, target_random_slope_index):
        if idx is not None and not 0 <= int(idx) < p:
            raise ValueError("random-slope index outside predictor matrix")
    if random_slope_sd_per_amplitude < 0:
        raise ValueError("random_slope_sd_per_amplitude must be non-negative")

    ns = prepared.absorber.n_source
    nt = prepared.absorber.n_target
    n = len(s)
    fixed_term = x @ fixed
    rows: list[dict[str, float | int]] = []

    for amplitude in tuple(float(v) for v in private_amplitudes):
        if amplitude < 0:
            raise ValueError("private amplitudes must be non-negative")
        slope_sd = random_slope_sd_per_amplitude * amplitude
        for effect in tuple(float(v) for v in effects):
            rng = np.random.default_rng(_cell_seed(seed_namespace, amplitude, effect))
            rejected = 0
            coef_sum = 0.0
            se_values: list[float] = []
            completed = 0
            while completed < worlds_per_cell:
                w = min(block_size, worlds_per_cell - completed)
                source_intercept = rng.normal(0.0, source_intercept_sd, size=(ns, w))
                target_intercept = rng.normal(0.0, target_intercept_sd, size=(nt, w))
                noise = rng.normal(0.0, dyad_noise_sd, size=(n, w))
                latent = (
                    source_intercept[s]
                    + target_intercept[t]
                    + fixed_term[:, None]
                    + effect * x[:, int(primary_index), None]
                    + noise
                )
                if source_random_slope_index is not None and slope_sd > 0:
                    slopes = rng.normal(0.0, slope_sd, size=(ns, w))
                    latent += slopes[s] * x[:, int(source_random_slope_index), None]
                if target_random_slope_index is not None and slope_sd > 0:
                    slopes = rng.normal(0.0, slope_sd, size=(nt, w))
                    latent += slopes[t] * x[:, int(target_random_slope_index), None]
                response = np.tanh(latent) if response_transform == "tanh" else latent
                fit = batch_primary_test(prepared, response)
                if not (
                    np.isfinite(fit.coefficient).all()
                    and np.isfinite(fit.standard_error).all()
                    and np.isfinite(fit.p_value_one_sided).all()
                ):
                    raise RuntimeError("non-finite synthetic qualification world")
                rejected += int(np.count_nonzero(fit.p_value_one_sided <= alpha))
                coef_sum += float(np.sum(fit.coefficient))
                se_values.extend(map(float, fit.standard_error))
                completed += w

            lo, hi = wilson_interval(rejected, worlds_per_cell)
            rows.append({
                "private_amplitude": amplitude,
                "random_slope_sd": float(slope_sd),
                "effect": effect,
                "worlds": int(worlds_per_cell),
                "rejections": int(rejected),
                "rejection_rate": float(rejected / worlds_per_cell),
                "wilson95_lower": float(lo),
                "wilson95_upper": float(hi),
                "coefficient_mean": float(coef_sum / worlds_per_cell),
                "standard_error_median": float(np.median(se_values)),
            })
    return rows


def minimum_detectable_effects(
    surface: Sequence[dict[str, float | int]],
    *,
    target_power: float = 0.80,
) -> dict[str, float | None]:
    """Return the smallest positive grid effect whose Wilson lower bound reaches target power."""
    if not 0 < target_power < 1:
        raise ValueError("target_power must lie in (0, 1)")
    by_amplitude: dict[float, list[dict[str, float | int]]] = {}
    for row in surface:
        amplitude = float(row["private_amplitude"])
        by_amplitude.setdefault(amplitude, []).append(row)
    out: dict[str, float | None] = {}
    for amplitude, rows in sorted(by_amplitude.items()):
        eligible = sorted(
            (
                float(row["effect"])
                for row in rows
                if float(row["effect"]) > 0
                and float(row["wilson95_lower"]) >= target_power
            )
        )
        out[f"{amplitude:g}"] = eligible[0] if eligible else None
    return out


def calibrated_detectability_envelope(
    surface: Sequence[dict[str, float | int]],
    *,
    target_power: float = 0.80,
    type1_wilson_upper_max: float = 0.05,
) -> dict[str, dict[str, float | bool | None]]:
    """Return per-heterogeneity detectability only where null calibration passes.

    A minimum detectable effect is meaningful only when the same inferential
    procedure controls false positives under the corresponding nuisance regime.
    This function therefore pairs each raw MDE with the beta=0 calibration cell
    at the same private-amplitude level and suppresses the MDE when the Wilson
    upper bound exceeds the declared Type-I ceiling.
    """
    if not 0 < target_power < 1:
        raise ValueError("target_power must lie in (0, 1)")
    if not 0 < type1_wilson_upper_max < 1:
        raise ValueError("type1_wilson_upper_max must lie in (0, 1)")

    by_amplitude: dict[float, list[dict[str, float | int]]] = {}
    for row in surface:
        amplitude = float(row["private_amplitude"])
        by_amplitude.setdefault(amplitude, []).append(row)

    out: dict[str, dict[str, float | bool | None]] = {}
    for amplitude, rows in sorted(by_amplitude.items()):
        null_rows = [row for row in rows if float(row["effect"]) == 0.0]
        if len(null_rows) != 1:
            raise ValueError(
                f"expected exactly one beta=0 calibration cell for A={amplitude:g}"
            )
        null = null_rows[0]
        calibration_pass = (
            float(null["wilson95_upper"]) <= type1_wilson_upper_max
        )
        eligible = sorted(
            float(row["effect"])
            for row in rows
            if float(row["effect"]) > 0
            and float(row["wilson95_lower"]) >= target_power
        )
        raw_mde = eligible[0] if eligible else None
        out[f"{amplitude:g}"] = {
            "null_rejection_rate": float(null["rejection_rate"]),
            "null_wilson95_upper": float(null["wilson95_upper"]),
            "calibration_pass": bool(calibration_pass),
            "raw_grid_mde": raw_mde,
            "evaluable_grid_mde": raw_mde if calibration_pass else None,
        }
    return out
