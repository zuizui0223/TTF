from __future__ import annotations

from dataclasses import asdict

import numpy as np

from .relational_qualification import (
    dyadic_signal_support,
    residualize_primary_relation,
)


def _zscore(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) == 0 or not np.isfinite(x).all():
        raise ValueError("expected finite non-empty vector")
    sd = float(x.std())
    if sd <= np.finfo(float).eps:
        raise ValueError("constant vector")
    return (x - float(x.mean())) / sd


def nested_historical_geographic_information(
    *,
    source_index: np.ndarray,
    target_index: np.ndarray,
    r_hist: np.ndarray,
    r_current: np.ndarray,
    geographic_coverage: np.ndarray,
    same_class: np.ndarray,
    same_order: np.ndarray,
    same_family: np.ndarray,
    occurrence_count_ratio: np.ndarray,
) -> dict[str, object]:
    """Localize historical-relation information loss attributable to geography.

    C1 and C2 are fit on exactly the same dyads. C2 adds only standardized
    directed external geographic co-opportunity to the C1 control set.
    Therefore the residual-sum-of-squares ratio is an exact nested information
    survival fraction, not a post-hoc comparison of unrelated models.
    """
    s = np.asarray(source_index, dtype=np.int64)
    t = np.asarray(target_index, dtype=np.int64)
    rh = np.asarray(r_hist, dtype=float)
    rc = np.asarray(r_current, dtype=float)
    geo = np.asarray(geographic_coverage, dtype=float)
    sc = np.asarray(same_class, dtype=float)
    so = np.asarray(same_order, dtype=float)
    sf = np.asarray(same_family, dtype=float)
    ratio = np.asarray(occurrence_count_ratio, dtype=float)

    vectors = (s, t, rh, rc, geo, sc, so, sf, ratio)
    if any(x.ndim != 1 for x in vectors):
        raise ValueError("all C2 inputs must be vectors")
    if not all(len(x) == len(s) for x in vectors) or len(s) == 0:
        raise ValueError("all C2 inputs must align and be non-empty")
    if not all(np.isfinite(x).all() for x in (rh, rc, geo, sc, so, sf, ratio)):
        raise ValueError("all numeric C2 inputs must be finite")
    if np.any(ratio <= 0):
        raise ValueError("occurrence count ratios must be positive")

    primary = _zscore(rh)
    baseline = np.column_stack((
        sc - sc.mean(),
        so - so.mean(),
        sf - sf.mean(),
        _zscore(np.abs(np.log(ratio))),
    ))
    c1_controls = np.column_stack((_zscore(rc), baseline))
    c2_controls = np.column_stack((_zscore(rc), _zscore(geo), baseline))

    c1_residual, c1 = residualize_primary_relation(
        s, t, primary, c1_controls
    )
    c2_residual, c2 = residualize_primary_relation(
        s, t, primary, c2_controls
    )
    c2_support = dyadic_signal_support(s, t, c2_residual)

    c1_ss = float(np.dot(c1_residual, c1_residual))
    c2_ss = float(np.dot(c2_residual, c2_residual))
    tolerance = 1e-10 * max(1.0, c1_ss)
    if c2_ss > c1_ss + tolerance:
        raise RuntimeError("nested C2 residual information exceeds C1")
    retained = 0.0 if c1_ss <= np.finfo(float).eps else c2_ss / c1_ss
    retained = float(np.clip(retained, 0.0, 1.0))

    return {
        "schema": "ttf_q_c2_nested_geographic_information_v0.1",
        "C1_recomputed": asdict(c1),
        "C2_geography_adjusted": asdict(c2),
        "C2_signal_support": asdict(c2_support),
        "nested_information": {
            "C1_residual_sum_squares": c1_ss,
            "C2_residual_sum_squares": c2_ss,
            "geography_retained_fraction": retained,
            "geography_attributable_fraction_of_C1_unique_information": (
                1.0 - retained
            ),
            "exact_nested_identity": (
                "C2_residual_SS = C1_residual_SS * geography_retained_fraction"
            ),
        },
        "conditioning": {
            "C1": [
                "current_climate_similarity",
                "same_class",
                "same_order",
                "same_family",
                "absolute_log_retained_environment_occurrence_count_ratio",
                "source_fixed_effect",
                "target_fixed_effect",
            ],
            "C2_adds_exactly": [
                "directed_external_geographic_co_opportunity_500km"
            ],
        },
        "interpretation": {
            "genetic_response_used": False,
            "legacy_C_reopened": False,
            "role": (
                "response-blind localization of historical ecological "
                "information loss attributable to present spatial co-opportunity"
            ),
        },
    }
