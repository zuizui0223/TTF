from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from typing import Iterable, Sequence

import numpy as np


PANEL_TAG = "butterfly-climate-release-independent-panel-v0.1"
PRIMARY_TAG = "butterfly-climate-release-partial-spearman-v0.1"
PRIMARY_TAG_V02 = "butterfly-climate-release-freedman-lane-v0.2"


@dataclass(frozen=True)
class ExpansionDescriptor:
    species: str
    host_family_count: float
    host_wgsrpd3_unit_count: int
    resolved_host_species: int
    wing_size_proxy: float
    voltinism: str


def host_breadth_stratum(host_family_count: float) -> str:
    value = float(host_family_count)
    if value == 1:
        return "1_family"
    if value == 2:
        return "2_families"
    if 3 <= value <= 5:
        return "3_to_5_families"
    if value >= 6:
        return "6plus_families"
    raise ValueError(f"invalid positive host-family count: {value}")


def _tie_key(species: str, tag: str) -> str:
    return hashlib.sha256(f"{tag}|{species}".encode("utf-8")).hexdigest()


def _zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = float(np.std(values))
    if sd <= np.sqrt(np.finfo(float).eps):
        return np.zeros_like(values)
    return (values - float(np.mean(values))) / sd


def farthest_point_species(
    descriptors: Sequence[ExpansionDescriptor],
    *,
    n: int,
    tag: str,
) -> tuple[str, ...]:
    if n < 1:
        raise ValueError("n must be positive")
    if len(descriptors) < n:
        raise ValueError(f"need at least {n} descriptors; got {len(descriptors)}")

    rows = sorted(
        descriptors,
        key=lambda d: (_tie_key(d.species, tag), d.species),
    )
    x1 = _zscore(
        np.log1p(np.asarray([d.host_wgsrpd3_unit_count for d in rows], dtype=float))
    )
    x2 = _zscore(
        np.log1p(np.asarray([d.resolved_host_species for d in rows], dtype=float))
    )
    x3 = _zscore(np.asarray([d.wing_size_proxy for d in rows], dtype=float))
    points = np.column_stack([x1, x2, x3])

    radius2 = np.sum(points * points, axis=1)
    first = min(
        range(len(rows)),
        key=lambda i: (-float(radius2[i]), _tie_key(rows[i].species, tag)),
    )
    selected = [first]
    remaining = set(range(len(rows))) - {first}

    while len(selected) < n:
        best = min(
            remaining,
            key=lambda i: (
                -min(
                    float(np.sum((points[i] - points[j]) ** 2))
                    for j in selected
                ),
                _tie_key(rows[i].species, tag),
            ),
        )
        selected.append(best)
        remaining.remove(best)

    return tuple(rows[i].species for i in selected)


def select_independent_panel(
    descriptors: Sequence[ExpansionDescriptor],
    *,
    excluded_species: Iterable[str],
    per_stratum: int = 8,
    minimum_native_units: int = 10,
    tag: str = PANEL_TAG,
) -> tuple[str, ...]:
    excluded = {str(x) for x in excluded_species}
    eligible = [
        d
        for d in descriptors
        if d.species not in excluded
        and math.isfinite(d.host_family_count)
        and d.host_family_count > 0
        and d.host_wgsrpd3_unit_count >= int(minimum_native_units)
        and d.resolved_host_species >= d.host_family_count
        and math.isfinite(d.wing_size_proxy)
    ]
    strata = {
        "1_family": [],
        "2_families": [],
        "3_to_5_families": [],
        "6plus_families": [],
    }
    for descriptor in eligible:
        strata[host_breadth_stratum(descriptor.host_family_count)].append(descriptor)

    selected: list[str] = []
    for label in ("1_family", "2_families", "3_to_5_families", "6plus_families"):
        rows = strata[label]
        if len(rows) < per_stratum:
            raise ValueError(
                f"stratum {label} has {len(rows)} eligible species; "
                f"need {per_stratum}"
            )
        selected.extend(
            farthest_point_species(
                rows,
                n=per_stratum,
                tag=f"{tag}|{label}",
            )
        )
    if len(selected) != 4 * per_stratum or len(set(selected)) != len(selected):
        raise RuntimeError("independent panel size/uniqueness drift")
    return tuple(selected)


def average_ranks(values: Sequence[float]) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and x[order[stop]] == x[order[start]]:
            stop += 1
        avg = 0.5 * ((start + 1) + stop)
        ranks[order[start:stop]] = avg
        start = stop
    return ranks


def _residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(x)), x])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ coef


def partial_spearman(
    response: Sequence[float],
    predictor: Sequence[float],
    control: Sequence[float],
) -> float:
    y = average_ranks(response)
    x = average_ranks(predictor)
    z = average_ranks(control)
    yr = _residualize(y, z)
    xr = _residualize(x, z)
    sy = float(np.std(yr))
    sx = float(np.std(xr))
    if sy <= np.sqrt(np.finfo(float).eps) or sx <= np.sqrt(np.finfo(float).eps):
        raise ValueError("partial Spearman residual has zero variance")
    return float(np.corrcoef(yr, xr)[0, 1])


def _fit_reduced(y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    design = np.column_stack([np.ones(len(z)), z])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coef
    return fitted, y - fitted


def deterministic_residual_donor_indices(
    species: Sequence[str],
    *,
    iteration: int,
    tag: str = PRIMARY_TAG_V02,
) -> np.ndarray:
    names = tuple(map(str, species))
    if len(set(names)) != len(names):
        raise ValueError("species identities must be unique")
    canonical_targets = tuple(sorted(names))
    donors = tuple(
        sorted(
            names,
            key=lambda name: (
                hashlib.sha256(
                    f"{tag}|{int(iteration)}|{name}".encode("utf-8")
                ).hexdigest(),
                name,
            ),
        )
    )
    donor_for_target = dict(zip(canonical_targets, donors))
    index_by_name = {name: i for i, name in enumerate(names)}
    return np.asarray(
        [index_by_name[donor_for_target[name]] for name in names],
        dtype=int,
    )


def freedman_lane_partial_spearman_permutation(
    response: Sequence[float],
    predictor: Sequence[float],
    control: Sequence[float],
    species: Sequence[str],
    *,
    iterations: int = 9999,
    tag: str = PRIMARY_TAG_V02,
) -> dict[str, float | int | str]:
    """One-sided conditional permutation test for ranked partial association.

    The response, predictor and control are rank transformed once. The predictor
    residual after removing the control is held fixed. Under the reduced null,
    response ranks are fitted on the control; only those reduced-model residuals
    are permuted, added back to the fixed reduced fitted values, re-residualized
    on the same control, and correlated with the fixed predictor residual.
    """
    if not (
        len(response) == len(predictor) == len(control) == len(species)
    ):
        raise ValueError("response, predictor, control and species lengths differ")
    if len(set(map(str, species))) != len(species):
        raise ValueError("species identities must be unique")
    if iterations < 1:
        raise ValueError("iterations must be positive")

    y = average_ranks(response)
    x = average_ranks(predictor)
    z = average_ranks(control)
    y_fitted, y_resid = _fit_reduced(y, z)
    x_resid = _residualize(x, z)

    sy = float(np.std(y_resid))
    sx = float(np.std(x_resid))
    if sy <= np.sqrt(np.finfo(float).eps) or sx <= np.sqrt(np.finfo(float).eps):
        raise ValueError("partial Spearman residual has zero variance")
    observed = float(np.corrcoef(y_resid, x_resid)[0, 1])

    names = tuple(map(str, species))
    extreme = 0
    for i in range(iterations):
        donor_indices = deterministic_residual_donor_indices(
            names,
            iteration=i,
            tag=tag,
        )
        pseudo_y = y_fitted + y_resid[donor_indices]
        pseudo_resid = _residualize(pseudo_y, z)
        stat = float(np.corrcoef(pseudo_resid, x_resid)[0, 1])
        extreme += int(stat <= observed)

    p_value = (extreme + 1) / (iterations + 1)
    return {
        "method": "Freedman-Lane-style reduced-response residual permutation on ranks",
        "observed_partial_spearman": observed,
        "iterations": int(iterations),
        "lower_tail_extreme_permutations": int(extreme),
        "one_sided_p_value": float(p_value),
    }


def deterministic_permutations(
    labels: Sequence[float],
    *,
    iterations: int,
    tag: str = PRIMARY_TAG,
) -> Iterable[np.ndarray]:
    values = np.asarray(labels, dtype=float)
    if iterations < 1:
        raise ValueError("iterations must be positive")
    for i in range(iterations):
        keys = [
            hashlib.sha256(f"{tag}|{i}|{j}".encode("utf-8")).hexdigest()
            for j in range(len(values))
        ]
        order = np.argsort(np.asarray(keys, dtype=object), kind="mergesort")
        yield values[order]


def one_sided_partial_spearman_permutation(
    response: Sequence[float],
    predictor: Sequence[float],
    control: Sequence[float],
    *,
    iterations: int = 9999,
    tag: str = PRIMARY_TAG,
) -> dict[str, float | int]:
    observed = partial_spearman(response, predictor, control)
    extreme = 0
    for permuted in deterministic_permutations(
        predictor,
        iterations=iterations,
        tag=tag,
    ):
        stat = partial_spearman(response, permuted, control)
        extreme += int(stat <= observed)
    p_value = (extreme + 1) / (iterations + 1)
    return {
        "observed_partial_spearman": float(observed),
        "iterations": int(iterations),
        "lower_tail_extreme_permutations": int(extreme),
        "one_sided_p_value": float(p_value),
    }


def species_digest(species: Sequence[str]) -> str:
    payload = "\n".join(map(str, species)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


__all__ = [
    "ExpansionDescriptor",
    "PANEL_TAG",
    "PRIMARY_TAG",
    "PRIMARY_TAG_V02",
    "average_ranks",
    "deterministic_residual_donor_indices",
    "farthest_point_species",
    "freedman_lane_partial_spearman_permutation",
    "host_breadth_stratum",
    "one_sided_partial_spearman_permutation",
    "partial_spearman",
    "select_independent_panel",
    "species_digest",
]
