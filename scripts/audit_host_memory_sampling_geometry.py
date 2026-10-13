#!/usr/bin/env python3
"""Exploratory predictor-vs-sampling-geometry audit on frozen HHM external data.

No genotype, nucleotide identity, edge genetic distance, genetic turnover or
species-specific exposed genetic coefficient is read. The frozen 259-species
one-shot was consumed/NOT_EVALUABLE; this audit does not reopen or rescue it.
The host predictor measures modern host-climate analogue distance projected
onto LGM rasters, NOT true LGM host occupancy or traversable host corridors.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ANALYSIS_STATUS = "POST_OUTCOME_RESPONSE_BLIND_PREDICTOR_GEOMETRY_AUDIT"
PROXY_SCOPE = "historical_host_climate_analogue_distance_not_ancient_host_occupancy_or_corridor"
PINS = {
    "predictors": "70b748d588d9ab10309ace080f85b8053c6e14aa96d7a84288724c1000272a5b",
    "roles": "4168e7572d3378e7f1595ae562d53c81bb2ef165cbe50eb2ed8ec4e7406fc6ae",
    "localities": "037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9",
    "edges": "ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9",
    "original_result": "5d994a76f7f20a778f0394f545c08a729182ce32298503c128d95c00da29f965",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def z(x):
    a = np.asarray(x, dtype=float)
    if a.ndim != 1 or len(a) == 0 or not np.isfinite(a).all():
        raise ValueError("nonfinite or empty one-dimensional metric")
    sd = a.std(ddof=0)
    if not sd > 1e-12:
        raise ValueError("invariant standardized metric")
    return (a - a.mean()) / sd


def sampling_edge_features(xyz, edge_nodes):
    """Separate contemporary sample spacing from frozen edge-to-edge geography.

    Site nearest-neighbor spacing is evaluated on ALL original localities,
    never the currently selected edge list, and never on insect genotypes.
    Graph degree is an additional *sampling topology* diagnostic only.
    """
    xyz = np.asarray(xyz, dtype=float)
    edge_nodes = np.asarray(edge_nodes)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 3 or not np.isfinite(xyz).all():
        raise ValueError("finite >=3 ECEF localities are required")
    if (edge_nodes.ndim != 2 or edge_nodes.shape[1] != 2 or len(edge_nodes) < 6
            or not np.issubdtype(edge_nodes.dtype, np.integer)):
        raise ValueError("frozen undirected graph edges required")
    if (np.any(edge_nodes < 0) or np.any(edge_nodes >= len(xyz))
            or np.any(edge_nodes[:, 0] == edge_nodes[:, 1])):
        raise ValueError("out-of-range locality or self-edge")
    if len(np.unique(np.sort(edge_nodes, axis=1), axis=0)) != len(edge_nodes):
        raise ValueError("duplicate undirected frozen edge")
    dist = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    np.fill_diagonal(dist, np.inf)
    nn = dist.min(axis=1)
    if not np.isfinite(nn).all() or np.any(nn <= 1e-9):
        raise ValueError("nonunique or nonfinite frozen geographic locality")
    degree = np.bincount(edge_nodes.ravel(), minlength=len(xyz))
    a, b = edge_nodes[:, 0], edge_nodes[:, 1]
    chord = np.linalg.norm(xyz[a] - xyz[b], axis=1)
    if not np.isfinite(chord).all() or np.any(chord <= 0):
        raise ValueError("invalid frozen chord distance")
    return {
        "chord_km": chord,
        "endpoint_nearest_neighbor_km": (nn[a] + nn[b]) / 2,
        "endpoint_graph_degree": (degree[a] + degree[b]) / 2,
    }


def partial_sampling_overlap(m_host, m_self, d_host_current, geography):
    """Finite-population partial R² of host-climate analogue versus sample isolation.

    Intercept, insect-own climatic history, present host mismatch, and direct
    endpoint chord distance are removed from BOTH the host predictor and
    endpoint nearest-neighbor spacing. Degrees add only a secondary R².
    No species' genetic outcome enters either fit. No p-value is computed.
    """
    mh = z(m_host)
    n = len(mh)
    if n < 10:
        raise ValueError("too few original frozen edge rows")
    other = np.asarray([m_self, d_host_current], dtype=float)
    if other.shape != (2, n):
        raise ValueError("control predictor shape differs from host edges")
    if any(np.asarray(geography[k]).shape != (n,) for k in (
        "chord_km", "endpoint_nearest_neighbor_km", "endpoint_graph_degree"
    )):
        raise ValueError("frozen edge geometry not aligned with host predictor")
    C = np.column_stack((np.ones(n), z(other[0]), z(other[1]),
                         z(np.log1p(geography["chord_km"]))))
    if np.linalg.matrix_rank(C, tol=1e-10) != C.shape[1]:
        raise ValueError("rank-deficient original climate/distance controls")
    rh = mh - C @ np.linalg.lstsq(C, mh, rcond=None)[0]
    if np.linalg.norm(rh) <= 1e-10:
        raise ValueError("history information absent after climate and distance")
    try:
        nn = z(np.log1p(geography["endpoint_nearest_neighbor_km"]))
    except ValueError as exc:
        raise ValueError("invariant sampling isolation") from exc
    rn = nn - C @ np.linalg.lstsq(C, nn, rcond=None)[0]
    if np.linalg.norm(rn) <= 1e-10:
        raise ValueError("sampling isolation fully explained by base controls")
    r = float((rh @ rn) / np.sqrt((rh @ rh) * (rn @ rn)))
    r = max(-1.0, min(1.0, r))
    var = rh @ rh
    # An additional graph-degree measure could be collinear with the NN proxy;
    # least squares projects to the span, without inventing a second degree of freedom.
    gd = np.asarray(geography["endpoint_graph_degree"], dtype=float)
    try:
        deg = z(np.log1p(gd))
        B = np.column_stack((C, nn, deg))
    except ValueError:
        B = np.column_stack((C, nn))
    remain = mh - B @ np.linalg.lstsq(B, mh, rcond=None)[0]
    r2joint = float(max(0.0, min(1.0, 1 - (remain @ remain) / var)))
    # Ordinary residual degrees-of-freedom adjustment is descriptive only:
    # edges share nodes and are not independent statistical observations.
    rank_base = int(np.linalg.matrix_rank(C))
    rank_joint = int(np.linalg.matrix_rank(B))
    if n <= rank_joint:
        raise ValueError("insufficient original edges for adjusted overlap diagnostic")
    nn_adjusted_r2 = float(1.0 - (1.0 - r*r) * (n-rank_base)/(n-(rank_base+1)))
    joint_adjusted_r2 = float(1.0 - (1.0 - r2joint) * (n-rank_base)/(n-rank_joint))
    if r2joint < r*r - 1e-9:
        raise RuntimeError("sampling geometry joint projection less than its nearest-neighbor subspace")
    return {
        "nearest_neighbor_partial_r": r,
        "nearest_neighbor_partial_r2": r*r,
        "sampling_geometry_joint_partial_r2": r2joint,
        "nearest_neighbor_df_adjusted_partial_r2": nn_adjusted_r2,
        "sampling_geometry_df_adjusted_partial_r2": joint_adjusted_r2,
        "postclimate_chord_unique_history_fraction": float(var/(mh @ mh)),
        "n_frozen_edges": n,
    }


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def describe(values):
    x=np.asarray(values,float)
    if len(x)==0 or not np.isfinite(x).all():
        raise RuntimeError("empty or nonfinite original predictor geometry subgroup")
    return {"species":int(len(x)),"median":float(np.median(x)),
            "mean":float(x.mean()),"p90":float(np.quantile(x,0.90)),
            "fraction_above_0_10":float(np.mean(x>0.10)),
            "species_above_0_10":int(np.sum(x>0.10)),
            "fraction_above_0_25":float(np.mean(x>0.25)),
            "species_above_0_25":int(np.sum(x>0.25))}


def audit(frozen_predictors, frozen_edges, frozen_localities, panel_roles):
    preds=defaultdict(dict);edges=defaultdict(dict);loc=defaultdict(dict)
    for r in frozen_predictors:
        sp=r["species"];idx=int(r["edge_index"])
        if idx in preds[sp]: raise RuntimeError("duplicate source predictor index")
        preds[sp][idx]=r
    for r in frozen_edges:
        sp=r["species"];idx=int(r["edge_index"])
        if sp not in preds: continue
        if idx in edges[sp]: raise RuntimeError("duplicate source edge index")
        edges[sp][idx]=(int(r["node_left"]),int(r["node_right"]))
    for r in frozen_localities:
        sp=r["species"];idx=int(r["locality_index"])
        if sp not in preds: continue
        if idx in loc[sp]: raise RuntimeError("duplicate source locality index")
        loc[sp][idx]=(float(r["x_km"]),float(r["y_km"]),float(r["z_km"]))
    roles={r["species"]:r["panel"] for r in panel_roles}
    if len(roles)!=len(panel_roles) or set(roles)!=set(preds):
        raise RuntimeError("frozen panel or external species changed")
    if set(roles.values()) != {"development", "confirmatory"}:
        raise RuntimeError("unexpected response-blind role value")
    if set(loc)!=set(preds) or set(edges)!=set(preds):
        raise RuntimeError("frozen locality/edge coverage absent")
    found=[]
    edge_count=0
    for sp in sorted(roles):
        idx=sorted(preds[sp])
        if idx!=list(range(len(idx))) or sorted(edges[sp])!=idx:
            raise RuntimeError(f"frozen host edge indexes differ from geometry: {sp}")
        nodes=np.asarray([edges[sp][i] for i in idx],dtype=np.int64)
        if sorted(loc[sp])!=list(range(len(loc[sp]))):
            raise RuntimeError("noncontiguous original ECEF sample localities")
        geo=sampling_edge_features(np.array([loc[sp][i] for i in sorted(loc[sp])]),nodes)
        rows=[preds[sp][i] for i in idx]
        vals=lambda key:np.array([float(r[key]) for r in rows],float)
        mh=vals("M_host");cur=vals("D_host_0BP");lgm=vals("D_host_LGM")
        if not np.allclose(mh,lgm-cur,rtol=0,atol=1e-10):
            raise RuntimeError("original climatic-memory edge identity violated")
        result=partial_sampling_overlap(mh,vals("M_self"),cur,geo)
        found.append({"panel":roles[sp],**result})
        edge_count+=len(idx)
    panels={}
    for panel in ("development", "confirmatory", "all"):
        subset=[x for x in found if panel=="all" or x["panel"]==panel]
        panels[panel]={
            "nearest_neighbor_partial_r2":describe([v["nearest_neighbor_partial_r2"] for v in subset]),
            "joint_sampling_geometry_partial_r2":describe([v["sampling_geometry_joint_partial_r2"] for v in subset]),
            "nearest_neighbor_df_adjusted_partial_r2":describe([v["nearest_neighbor_df_adjusted_partial_r2"] for v in subset]),
            "joint_sampling_geometry_df_adjusted_partial_r2":describe([v["sampling_geometry_df_adjusted_partial_r2"] for v in subset]),
            "host_information_after_climate_and_chord":describe([v["postclimate_chord_unique_history_fraction"] for v in subset]),
        }
    strata=[]
    for lo,hi in ((0,30),(30,100),(100,None)):
        subset=[x for x in found if lo<=x["n_frozen_edges"] and (hi is None or x["n_frozen_edges"]<hi)]
        strata.append({
            "min_original_edges_inclusive":lo,"max_original_edges_exclusive":hi,
            "species":len(subset),
            "median_nearest_neighbor_partial_r2":float(np.median([v["nearest_neighbor_partial_r2"] for v in subset])) if subset else None,
            "median_nearest_neighbor_df_adjusted_partial_r2":float(np.median([v["nearest_neighbor_df_adjusted_partial_r2"] for v in subset])) if subset else None,
            "median_joint_sampling_geometry_partial_r2":float(np.median([v["sampling_geometry_joint_partial_r2"] for v in subset])) if subset else None,
        })
    return {
        "schema":"ttf_host_climate_memory_sampling_geometry_overlap_v0.1",
        "original_graph_edge_count_strata":strata,
        "status":ANALYSIS_STATUS,"source_bound_species":len(found),
        "response_blind_external_edge_rows":edge_count,
        "panels":panels,
        "interpretation":{
            "host_memory_proxy_scope":PROXY_SCOPE,
            "partial_r2_is_not_a_genetic_effect_explanatory_fraction":True,
            "coefficient_or_p_value_tested":False,
            "independently_confirmatory_genetic_panel_used":False,
            "original_259_empirical_status":"NOT_EVALUABLE_EMPIRICAL_RESPONSE",
            "original_259_genetic_response_not_reopened":True,
            "sampling_isolation_is_observed_sampling_spacing_not_historical_geographic_barrier":True,
            "edge_rows_within_species_not_independent_ecological_replicates":True,
            "this_response_exposed_source_audit_is_exploratory":True,
            "no_per_species_genotype_or_coefficient_exported":True,
        },
    }


def main():
    parser=argparse.ArgumentParser()
    for name in ("predictors", "roles", "edges", "localities", "original_result"):
        parser.add_argument("--"+name.replace("_","-"), type=Path, required=True)
    parser.add_argument("--output",type=Path,required=True)
    a=parser.parse_args()
    for name,expected in PINS.items():
        if sha256(getattr(a,name)) != expected:
            raise RuntimeError(f"exact previously frozen source SHA256 mismatch: {name}")
    original=json.loads(a.original_result.read_text())
    if not (original.get("schema")=="ttf_hhm_259_one_shot_empirical_v0.1"
            and original.get("status")=="NOT_EVALUABLE_EMPIRICAL_RESPONSE"
            and original.get("empirical_identity_opened") is True):
        raise RuntimeError("original consumed one-shot outcome not preserved")
    report=audit(read_csv(a.predictors),read_csv(a.edges),read_csv(a.localities),read_csv(a.roles))
    if (report["source_bound_species"],report["response_blind_external_edge_rows"],
        report["panels"]["development"]["nearest_neighbor_partial_r2"]["species"],
        report["panels"]["confirmatory"]["nearest_neighbor_partial_r2"]["species"]) != (641,71520,320,321):
        raise RuntimeError("previous response-blind 641-species source universe drift")
    report["source_sha256"]=PINS
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"all_species":641,
        "nearest_neighbor_partial_r2":report["panels"]["all"]["nearest_neighbor_partial_r2"]},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
