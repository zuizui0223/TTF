from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .relational_qualification import (
    dyadic_signal_support,
    residualize_primary_relation,
)


@dataclass(frozen=True)
class SurvivorInformationGateResult:
    dyads: int
    source_clusters: int
    target_clusters: int
    control_rank_after_two_way_fe: int
    endpoint_retained_variance_fraction: float
    control_unique_variance_fraction_after_fe: float
    total_unique_variance_fraction: float
    partial_sd_in_primary_sd_units: float
    signal_effective_sources: float
    signal_effective_targets: float
    max_source_signal_share: float
    max_target_signal_share: float
    total_unique_fraction_minimum: float
    maximum_endpoint_signal_share_ceiling: float
    unique_fraction_pass: bool
    source_concentration_pass: bool
    target_concentration_pass: bool
    overall_pass: bool


def _zscore(values: np.ndarray, *, name: str) -> np.ndarray:
    x=np.asarray(values,dtype=float)
    if x.ndim!=1 or len(x)==0 or not np.isfinite(x).all():
        raise ValueError(f"{name} must be a finite non-empty vector")
    sd=float(np.std(x,ddof=0))
    if sd<=np.finfo(float).eps:
        raise ValueError(f"{name} is constant under survivor geometry")
    return (x-float(np.mean(x)))/sd


def evaluate_survivor_information_gate(
    source_index: np.ndarray,
    target_index: np.ndarray,
    relation_g: np.ndarray,
    controls: np.ndarray,
    *,
    total_unique_fraction_minimum: float=0.15,
    maximum_endpoint_signal_share_ceiling: float=0.15,
) -> SurvivorInformationGateResult:
    s=np.asarray(source_index)
    t=np.asarray(target_index)
    g=np.asarray(relation_g,dtype=float)
    c=np.asarray(controls,dtype=float)
    if s.ndim!=1 or t.ndim!=1 or g.ndim!=1 or len(s)!=len(t) or len(s)!=len(g) or len(g)==0:
        raise ValueError("source, target, and G must be aligned non-empty vectors")
    if c.ndim!=2 or c.shape!=(len(g),4) or not np.isfinite(c).all():
        raise ValueError("survivor controls must be a finite n x 4 matrix")
    if len(np.unique(s))<2 or len(np.unique(t))<2:
        raise ValueError("survivor geometry requires at least two source and target species")
    if not 0 < float(total_unique_fraction_minimum) < 1:
        raise ValueError("invalid total unique fraction threshold")
    if not 0 < float(maximum_endpoint_signal_share_ceiling) < 1:
        raise ValueError("invalid endpoint concentration ceiling")

    z_g=_zscore(g,name="G")
    z_controls=np.column_stack([
        _zscore(c[:,j],name=f"control_{j}")
        for j in range(c.shape[1])
    ])
    residual,summary=residualize_primary_relation(s,t,z_g,z_controls)
    support=dyadic_signal_support(s,t,residual)

    unique_pass=bool(
        summary.total_unique_variance_fraction >= float(total_unique_fraction_minimum)
    )
    source_pass=bool(
        support.max_source_signal_share <= float(maximum_endpoint_signal_share_ceiling)
    )
    target_pass=bool(
        support.max_target_signal_share <= float(maximum_endpoint_signal_share_ceiling)
    )
    return SurvivorInformationGateResult(
        dyads=int(summary.dyads),
        source_clusters=int(summary.source_clusters),
        target_clusters=int(summary.target_clusters),
        control_rank_after_two_way_fe=int(summary.control_rank_after_two_way_fe),
        endpoint_retained_variance_fraction=float(summary.source_target_fe_retained_variance_fraction),
        control_unique_variance_fraction_after_fe=float(summary.control_unique_variance_fraction_after_fe),
        total_unique_variance_fraction=float(summary.total_unique_variance_fraction),
        partial_sd_in_primary_sd_units=float(summary.partial_sd_in_primary_sd_units),
        signal_effective_sources=float(support.signal_effective_sources),
        signal_effective_targets=float(support.signal_effective_targets),
        max_source_signal_share=float(support.max_source_signal_share),
        max_target_signal_share=float(support.max_target_signal_share),
        total_unique_fraction_minimum=float(total_unique_fraction_minimum),
        maximum_endpoint_signal_share_ceiling=float(maximum_endpoint_signal_share_ceiling),
        unique_fraction_pass=unique_pass,
        source_concentration_pass=source_pass,
        target_concentration_pass=target_pass,
        overall_pass=bool(unique_pass and source_pass and target_pass),
    )


def survivor_information_gate_dict(*args,**kwargs) -> dict:
    return asdict(evaluate_survivor_information_gate(*args,**kwargs))


__all__=[
    "SurvivorInformationGateResult",
    "evaluate_survivor_information_gate",
    "survivor_information_gate_dict",
]
