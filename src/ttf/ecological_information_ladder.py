from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .relational_qualification import (
    dyadic_signal_support,
    residualize_primary_relation,
)


@dataclass(frozen=True)
class RelationLayer:
    label: str
    raw_spearman_with_other: float | None
    total_unique_variance_fraction: float
    endpoint_retained_variance_fraction: float
    control_unique_variance_fraction_after_fe: float
    partial_sd_in_primary_sd_units: float
    signal_effective_sources: float
    signal_effective_targets: float
    max_source_signal_share: float
    max_target_signal_share: float


def _zscore(values: np.ndarray) -> np.ndarray:
    x=np.asarray(values,float)
    if x.ndim!=1 or len(x)==0 or not np.isfinite(x).all():
        raise ValueError("expected finite non-empty vector")
    sd=float(x.std())
    if sd<=np.finfo(float).eps:
        raise ValueError("constant vector")
    return (x-float(x.mean()))/sd


def _ranks(values: np.ndarray) -> np.ndarray:
    x=np.asarray(values,float)
    order=np.argsort(x,kind="stable")
    out=np.empty(len(x),float)
    sx=x[order]
    starts=np.r_[0,1+np.flatnonzero(sx[1:]!=sx[:-1])]
    stops=np.r_[starts[1:],len(x)]
    for start,stop in zip(starts,stops):
        out[order[start:stop]]=0.5*((start+1)+stop)
    return out


def _spearman(a: np.ndarray,b: np.ndarray) -> float:
    x=_ranks(a); y=_ranks(b)
    x-=x.mean(); y-=y.mean()
    den=float(np.sqrt(np.dot(x,x)*np.dot(y,y)))
    return 0.0 if den<=np.finfo(float).eps else float(np.dot(x,y)/den)


def _baseline_controls(
    *,
    class_name: np.ndarray,
    order_name: np.ndarray,
    family_name: np.ndarray,
    occurrence_n: np.ndarray,
    source_index: np.ndarray,
    target_index: np.ndarray,
) -> np.ndarray:
    cls=np.asarray(class_name).astype(str)
    order=np.asarray(order_name).astype(str)
    family=np.asarray(family_name).astype(str)
    occ=np.asarray(occurrence_n,float)
    s=np.asarray(source_index,np.int64)
    t=np.asarray(target_index,np.int64)
    same_class=(cls[s]==cls[t]).astype(float)
    same_order=(order[s]==order[t]).astype(float)
    same_family=(family[s]==family[t]).astype(float)
    locality_ratio=occ[s]/occ[t]
    return np.column_stack((
        same_class-same_class.mean(),
        same_order-same_order.mean(),
        same_family-same_family.mean(),
        _zscore(np.abs(np.log(locality_ratio))),
    ))


def _layer(
    label: str,
    source: np.ndarray,
    target: np.ndarray,
    primary: np.ndarray,
    controls: np.ndarray,
    other: np.ndarray | None,
) -> tuple[RelationLayer,np.ndarray]:
    residual,summary=residualize_primary_relation(
        source,target,_zscore(primary),controls
    )
    support=dyadic_signal_support(source,target,residual)
    return RelationLayer(
        label=label,
        raw_spearman_with_other=(
            None if other is None else _spearman(primary,other)
        ),
        total_unique_variance_fraction=float(summary.total_unique_variance_fraction),
        endpoint_retained_variance_fraction=float(summary.source_target_fe_retained_variance_fraction),
        control_unique_variance_fraction_after_fe=float(summary.control_unique_variance_fraction_after_fe),
        partial_sd_in_primary_sd_units=float(summary.partial_sd_in_primary_sd_units),
        signal_effective_sources=float(support.signal_effective_sources),
        signal_effective_targets=float(support.signal_effective_targets),
        max_source_signal_share=float(support.max_source_signal_share),
        max_target_signal_share=float(support.max_target_signal_share),
    ),residual


def ecological_information_ladder(
    *,
    class_name: np.ndarray,
    order_name: np.ndarray,
    family_name: np.ndarray,
    occurrence_n: np.ndarray,
    source_index: np.ndarray,
    target_index: np.ndarray,
    r_current: np.ndarray,
    r_hist: np.ndarray,
) -> dict[str,object]:
    """Nested response-blind B1 -> C1 environmental-information decomposition.

    B1 asks how much present-climate relation survives endpoint identity,
    lineage and locality imbalance. C1 then asks how much historical relation
    survives the same baseline structure *and* the present-climate relation.
    The layers are descriptive information diagnostics, not empirical outcome
    tests and not additive variance partitions of a biological response.
    """
    s=np.asarray(source_index,np.int64)
    t=np.asarray(target_index,np.int64)
    rc=np.asarray(r_current,float)
    rh=np.asarray(r_hist,float)
    if not (
        s.ndim==t.ndim==rc.ndim==rh.ndim==1
        and len(s)==len(t)==len(rc)==len(rh)>0
        and np.isfinite(rc).all() and np.isfinite(rh).all()
    ):
        raise ValueError("aligned finite dyad vectors are required")

    base=_baseline_controls(
        class_name=class_name,
        order_name=order_name,
        family_name=family_name,
        occurrence_n=occurrence_n,
        source_index=s,
        target_index=t,
    )
    current,current_residual=_layer(
        "B1_current_climate",
        s,t,rc,base,rh,
    )
    hist_controls=np.column_stack((_zscore(rc),base))
    historical,historical_residual=_layer(
        "C1_history_given_current",
        s,t,rh,hist_controls,rc,
    )

    baseline_hist_residual,_=residualize_primary_relation(
        s,t,_zscore(rh),base
    )
    baseline_current_residual,_=residualize_primary_relation(
        s,t,_zscore(rc),base
    )
    den=float(np.sqrt(
        np.dot(baseline_hist_residual,baseline_hist_residual)
        *np.dot(baseline_current_residual,baseline_current_residual)
    ))
    partial_corr=(
        0.0 if den<=np.finfo(float).eps
        else float(np.dot(baseline_hist_residual,baseline_current_residual)/den)
    )

    return {
        "schema":"ttf_q_ecological_information_ladder_v0.1",
        "layers":{
            "B1_current_climate":asdict(current),
            "C1_history_given_current":asdict(historical),
        },
        "cross_layer":{
            "raw_spearman_hist_vs_current":_spearman(rh,rc),
            "partial_correlation_hist_current_after_endpoint_lineage_locality":partial_corr,
            "history_fraction_remaining_after_adding_current_to_baseline":(
                float(np.dot(historical_residual,historical_residual))
                / float(np.dot(baseline_hist_residual,baseline_hist_residual))
                if float(np.dot(baseline_hist_residual,baseline_hist_residual))>np.finfo(float).eps
                else 0.0
            ),
        },
        "interpretation":{
            "B1":"present-climate pairwise information after endpoint identity, lineage and locality imbalance",
            "C1":"historical-climate pairwise information added beyond the B1 present-climate relation and the same baseline structure",
            "layers_are_not_biological_response_variance_components":True,
            "legacy_B_or_C_reopened":False,
        },
    }
