#!/usr/bin/env python3
"""Frozen response-blind historical host-memory synthetic qualification v0.1.

No sequence characters, genetic distances, empirical post-IBD turnover or
empirical historical-host effects are read. The entire simulation is generated
from external predictor columns and exact response-blind graph geography.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

PREDICTOR_SHA = "70b748d588d9ab10309ace080f85b8053c6e14aa96d7a84288724c1000272a5b"
ROLE_SHA = "4168e7572d3378e7f1595ae562d53c81bb2ef165cbe50eb2ed8ec4e7406fc6ae"
GEOMETRY_LOCAL_SHA = "037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9"
GEOMETRY_EDGE_SHA = "ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9"
QUAL_RULE_SHA = "2ec2594510021da653598e22a3d0f4b9e4cb226398aadbf342967ee532909415"
SOURCE_SHA = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce61a5"
WORLD_SEEDS = {
    "reference": "hhm-synthetic-null-reference-v0.1",
    "evaluation": "hhm-synthetic-null-evaluation-v0.1",
    "positive": "hhm-synthetic-positive-v0.1",
}
BASIS_SEED = "hhm-synthetic-spatial-basis-v0.1"


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def seed(namespace: str, species: str) -> int:
    return int(hashlib.sha256(f"{namespace}|{SOURCE_SHA}|{species}".encode()).hexdigest()[:16], 16)


def z(v: np.ndarray, axis: int = 0) -> np.ndarray:
    a = np.asarray(v, dtype=float)
    mu = a.mean(axis=axis, keepdims=True)
    sd = a.std(axis=axis, keepdims=True, ddof=0)
    if np.any(~np.isfinite(sd)) or np.any(sd < 1e-12):
        raise RuntimeError("nonfinite or invariant within-species synthetic field")
    return (a - mu) / sd


def load_geometry(locations: Path, edges: Path, species: set[str]):
    loc = defaultdict(dict)
    with locations.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["species"] in species:
                loc[row["species"]][int(row["locality_index"])] = np.array(
                    [float(row["x_km"]), float(row["y_km"]), float(row["z_km"])]
                )
    graph = defaultdict(dict)
    with edges.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["species"] in species:
                ix = int(row["edge_index"])
                if ix in graph[row["species"]]:
                    raise RuntimeError("duplicate frozen graph edge")
                graph[row["species"]][ix] = (int(row["node_left"]), int(row["node_right"]))
    if set(loc) != species or set(graph) != species:
        raise RuntimeError("missing frozen locality or edge graph for confirmatory species")
    return loc, graph


def prepare_species(species: str, predictors: list[dict], loc: dict, graph: dict):
    """Precompute all geometry and Frisch-Waugh residualization operators."""
    node_keys = sorted(loc)
    node_map = {key: i for i, key in enumerate(node_keys)}
    xyz = np.stack([loc[n] for n in node_keys])
    order = sorted(graph)
    if [int(r["edge_index"]) for r in predictors] != order:
        raise RuntimeError("predictor rows do not align with frozen graph edge identities")
    endpoints = np.array([graph[k] for k in order], dtype=int)
    try:
        left = np.array([node_map[a] for a in endpoints[:, 0]], dtype=int)
        right = np.array([node_map[a] for a in endpoints[:, 1]], dtype=int)
    except KeyError as exc:
        raise RuntimeError("edge refers to a nonexistent frozen locality") from exc
    if np.any(left == right) or len(set(tuple(sorted(pair)) for pair in endpoints)) != len(order):
        raise RuntimeError("self or duplicate undirected frozen edge")
    cols = [z(np.array([float(r[k]) for r in predictors])) for k in
            ("M_host", "M_self", "D_host_0BP")]
    mh, ms, dc = [a.ravel() for a in cols]
    controls = np.column_stack((np.ones(len(order)), ms, dc))
    resid = mh - controls @ np.linalg.lstsq(controls, mh, rcond=None)[0]
    norm = float(resid @ mh)
    if not math.isfinite(norm) or norm < 1e-8:
        raise RuntimeError(f"host information collapses in synthetic geometry: {species}")

    chord = np.linalg.norm(xyz[left] - xyz[right], axis=1)
    x = np.log1p(chord)
    x = z(x).ravel()
    m, n = len(order), len(node_keys)
    incidence = csr_matrix(
        (np.ones(2*m), (np.r_[left, right], np.r_[np.arange(m), np.arange(m)])),
        shape=(n, m)
    )
    degree = np.asarray(incidence.sum(axis=1)).ravel()
    sx = np.asarray(incidence @ x).ravel()
    sx2 = np.asarray(incidence @ (x*x)).ravel()
    cnt = m - degree[left] - degree[right] + 1
    train_sum_x = np.sum(x) - sx[left] - sx[right] + x
    train_sum_xx = np.sum(x*x) - sx2[left] - sx2[right] + x*x
    denom = cnt*train_sum_xx - train_sum_x**2
    if np.any(cnt < 5) or np.any(denom <= 1e-9) or np.any(~np.isfinite(denom)):
        raise RuntimeError(f"endpoint-safe IBD geometry cannot be qualified: {species}")

    basis_rng = np.random.default_rng(seed(BASIS_SEED, species))
    freq = basis_rng.normal(size=(3, 8)) / 200.0
    phase = xyz @ freq
    basis = np.column_stack((np.sin(phase), np.cos(phase)))

    return dict(name=species, left=left, right=right, x=x, mh=mh, ms=ms, dc=dc,
                resid=resid, norm=norm, incidence=incidence, cnt=cnt,
                train_sum_x=train_sum_x, denom=denom, basis=basis)


def endpoint_safe_ibd(raw: np.ndarray, pre: dict) -> np.ndarray:
    """Leave two endpoints out of every edge-specific IBD fit, in O(edges*worlds)."""
    x = pre["x"]
    left = pre["left"]
    right = pre["right"]
    incidence = pre["incidence"]
    sy_nodes = incidence @ raw
    sxy_nodes = incidence @ (raw*x[:, None])
    sy = raw.sum(axis=0)[None, :] - sy_nodes[left] - sy_nodes[right] + raw
    sxy = np.sum(raw*x[:, None], axis=0)[None, :] - sxy_nodes[left] - sxy_nodes[right] + raw*x[:, None]
    cnt = pre["cnt"][:, None]
    sx = pre["train_sum_x"][:, None]
    slope = (cnt*sxy - sx*sy) / pre["denom"][:, None]
    intercept = (sy - slope*sx) / cnt
    return raw - (intercept + slope*x[:, None])


def simulate_betas(pre: dict, group: str, n_worlds: int, block_size: int = 64) -> np.ndarray:
    rng = np.random.default_rng(seed(WORLD_SEEDS[group], pre["name"]))
    out = np.empty(n_worlds, dtype=float)
    m = len(pre["x"])
    spatial_basis = pre["basis"]
    for start in range(0, n_worlds, block_size):
        stop = min(n_worlds, start+block_size)
        b = stop-start
        mean = 0.10 if group == "positive" else 0.0
        beta = rng.normal(loc=mean, scale=0.10, size=b)
        weights = rng.normal(size=(16, b))
        latent_nodes = spatial_basis @ weights / np.sqrt(16.0)
        spatial_edge = np.abs(latent_nodes[pre["left"]] - latent_nodes[pre["right"]])
        private = z(spatial_edge, axis=0)
        noise = rng.normal(size=(m, b))
        raw = (0.60*pre["x"][:, None] + 0.22*pre["ms"][:, None]
               -0.18*pre["dc"][:, None] + pre["mh"][:, None]*beta[None, :]
               +0.45*private + 0.65*noise)
        post = endpoint_safe_ibd(raw, pre)
        sd = post.std(axis=0, ddof=0)
        if np.any(sd <= 1e-12) or np.any(~np.isfinite(sd)):
            raise RuntimeError("nonfinite post-IBD synthetic response")
        out[start:stop] = (pre["resid"] @ post) / (pre["norm"]*sd)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("nonfinite synthetic species beta")
    return out


def statistic(beta_matrix: np.ndarray) -> np.ndarray:
    if beta_matrix.ndim != 2 or beta_matrix.shape[0] < 2:
        raise RuntimeError("studentized statistic needs multiple independent species")
    mu = beta_matrix.mean(axis=0)
    se = beta_matrix.std(axis=0, ddof=1) / np.sqrt(beta_matrix.shape[0])
    if np.any(se <= 0) or np.any(~np.isfinite(se)):
        raise RuntimeError("nonfinite cross-species synthetic SE")
    return mu / se


def wilson(count: int, n: int, z95: float = 1.959963984540054):
    phat = count/n
    d = 1.0 + z95*z95/n
    center = (phat + z95*z95/(2*n))/d
    spread = (z95*np.sqrt(phat*(1-phat)/n + z95*z95/(4*n*n)))/d
    return (float(center-spread), float(center+spread))


def assess(reference: np.ndarray, evaluation: np.ndarray, positive: np.ndarray):
    if len(reference) != 1999 or len(evaluation) != 500 or len(positive) != 500:
        raise RuntimeError("frozen Monte Carlo world counts changed")
    null_reference = np.sort(reference)
    def pvalues(stats):
        # Count reference values >= observed (one-sided) with exact finite simulation bound.
        greater_equal = len(reference) - np.searchsorted(null_reference, stats, side="left")
        return (1+greater_equal)/(1+len(reference))
    pe = pvalues(evaluation)
    pp = pvalues(positive)
    null_reject = int(np.sum(pe <= 0.05))
    positive_reject = int(np.sum(pp <= 0.05))
    nlo,nhi = wilson(null_reject, 500)
    plo,phi = wilson(positive_reject, 500)
    passed = nhi <= 0.10 and plo >= 0.80
    return {
        "null_rejections": null_reject, "null_rejection_rate": null_reject/500,
        "null_wilson_lower": nlo, "null_wilson_upper": nhi,
        "positive_rejections": positive_reject, "positive_power": positive_reject/500,
        "positive_wilson_lower": plo, "positive_wilson_upper": phi,
        "type1_gate_pass": bool(nhi <= 0.10),
        "positive_power_gate_pass": bool(plo >= 0.80),
        "decision": ("PASS_TO_CONFIRMATORY_CHARACTER_MASK_ONLY" if passed else
                     "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_SYNTHETIC_QUALIFICATION"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    for k in ("edges", "localities", "predictors", "roles", "qualification_rule", "world_contract"):
        ap.add_argument("--"+k.replace("_","-"), type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    pinned = {"edges": GEOMETRY_EDGE_SHA, "localities": GEOMETRY_LOCAL_SHA,
              "predictors": PREDICTOR_SHA, "roles": ROLE_SHA,
              "qualification_rule": QUAL_RULE_SHA}
    for k, digest in pinned.items():
        if file_sha(getattr(args, k)) != digest:
            raise RuntimeError(f"frozen synthetic source hash mismatch: {k}")
    contract = json.loads(args.world_contract.read_text())
    rule = json.loads(args.qualification_rule.read_text())
    if contract["schema"] != "ttf_historical_host_memory_synthetic_world_contract_v0.1":
        raise RuntimeError("wrong synthetic world contract")
    if contract["status"] != "DESIGN_FROZEN_BEFORE_ANY_SYNTHETIC_QUALIFICATION_OUTCOME":
        raise RuntimeError("world contract not prospectively frozen")
    if any([contract["null_reference_worlds"] != 1999, contract["null_evaluation_worlds"] != 500,
            contract["positive_worlds"] != 500, contract["alpha"] != 0.05]):
        raise RuntimeError("synthetic world counts or alpha drift")
    if rule["synthetic_qualification"]["positive_structure"].find("0.10") < 0:
        raise RuntimeError("frozen positive effect floor drift")
    with args.roles.open(newline="", encoding="utf-8") as f:
        assignments = list(csv.DictReader(f))
    confirmatory = {r["species"] for r in assignments if r["panel"] == "confirmatory"}
    development = {r["species"] for r in assignments if r["panel"] == "development"}
    if (len(confirmatory), len(development), len(confirmatory & development)) != (321, 320, 0):
        raise RuntimeError("frozen response-blind split drift")
    preds = defaultdict(list)
    with args.predictors.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["species"] in confirmatory:
                preds[r["species"]].append(r)
    if set(preds) != confirmatory:
        raise RuntimeError("incomplete confirmatory external predictor")
    loc, graph = load_geometry(args.localities, args.edges, confirmatory)
    groups = {"reference": 1999, "evaluation": 500, "positive": 500}
    betas = {k: np.empty((321, n), dtype=float) for k,n in groups.items()}
    for i,sp in enumerate(sorted(confirmatory)):
        pre = prepare_species(sp, preds[sp], loc[sp], graph[sp])
        for group, n_world in groups.items():
            betas[group][i, :] = simulate_betas(pre, group, n_world, block_size=64)
        if (i+1) % 25 == 0:
            print(json.dumps({"processed_confirmatory_species":i+1,"total":321}), flush=True)
    stats = {k:statistic(v) for k,v in betas.items()}
    result = assess(stats["reference"],stats["evaluation"],stats["positive"])
    receipt = {
        "schema": "ttf_historical_host_memory_synthetic_qualification_v0.1",
        "contract_sha256": file_sha(args.world_contract),
        "predictor_sha256": PREDICTOR_SHA, "roles_sha256": ROLE_SHA,
        "confirmatory_species": 321,
        "frozen_groups": groups, "calibration": result,
        "no_empirical_genetic_response_read": True,
        "development_genetic_response_permanently_closed": True,
        "nucleotide_identity_opened": False,
        "genetic_distance_opened": False,
        "post_IBD_empirical_turnover_opened": False,
        "empirical_beta_host_opened": False,
        "interpretation": "Qualification of synthetic detection/error only; NOT an empirical ecological genetic association."
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
    # A scientific failed gate is a terminal recorded receipt, not process failure.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
