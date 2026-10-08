#!/usr/bin/env python3
"""Build the response-blind historical host-memory predictor.

This program never reads FASTA sequence characters or any genetic response.
It combines frozen genetic geometry with external host/range/climate assets.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from ttf.historical_host_memory import (
    great_circle_quadrature_latlon,
    historical_analog_memory,
    predictor_design_diagnostics,
    standardize_from_reference,
)

ASSET_RE = re.compile(r"CHELSA_TraCE21K_(bio0?1|bio0?7|bio12|bio15)_(-?\d+)_V1\.0\.tif$", re.I)
VARIABLE_ORDER=("bio01","bio07","bio12","bio15")

# Exact response-blind sources, frozen before host/climate predictor values.
FROZEN_INPUT_SHA256 = {
    "candidates": "c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b",
    "localities": "037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9",
    "edges": "ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9",
    "host_pairs": "0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9",
    "native_units": "c731906315f7452f83ce302c1d36c75ad94afc24ce7e24d461360312244ac558",
    "wgsrpd_support": "d0fc12f635ec56a442dd37dec06cf3aa06684a00b14f1bb1f2db2f6fe3c4a21f",
}



def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def norm_var(x: str) -> str:
    x=x.lower()
    if x=="bio1": return "bio01"
    if x=="bio7": return "bio07"
    return x


def norm_id(x: str) -> str:
    s=str(x).strip()
    if s.endswith(".0") and s[:-2].isdigit():
        s=s[:-2]
    return s



def verify_frozen_input_hashes(
    candidates: Path, localities: Path, edges: Path, host_pairs: Path,
    native_units: Path, wgsrpd_support: Path,
) -> None:
    """Refuse any unbound original source before reading biological predictors."""
    paths = (("candidates", candidates), ("localities", localities), ("edges", edges),
             ("host_pairs", host_pairs), ("native_units", native_units),
             ("wgsrpd_support", wgsrpd_support))
    for key, path in paths:
        if sha256_path(path) != FROZEN_INPUT_SHA256[key]:
            raise RuntimeError(f"frozen historical-host-memory {key} SHA256 mismatch")


def verify_host_pair_identity(
    candidates: list[dict[str, str]], pair_rows: list[dict[str, str]]
) -> None:
    """Require the exact candidate insect x accepted-host mapping in HOSTS/WCVP.

    The full sidecar can contain other insects; only candidate associations are
    compared. All IDs are normalized using the same WCVP numeric-ID rule.
    """
    mapping: dict[str, set[str]] = defaultdict(set)
    for row in pair_rows:
        species = str(row["insect_species"]).strip()
        host_id = norm_id(row["accepted_plant_name_id"])
        if not species or not host_id:
            raise RuntimeError("empty insect or accepted-host ID in frozen sidecar")
        mapping[species].add(host_id)
    for cand in candidates:
        species = str(cand["species"]).strip()
        frozen_ids = [norm_id(part) for part in cand["accepted_host_ids"].split(";")]
        if not frozen_ids or not all(frozen_ids) or len(set(frozen_ids)) != len(frozen_ids):
            raise RuntimeError(f"invalid frozen host IDs for {species}")
        if mapping.get(species, set()) != set(frozen_ids):
            raise RuntimeError(f"HOSTS/WCVP accepted-host mapping drift for {species}")


def load_assets(paths: list[Path]) -> dict[tuple[str,int],Path]:
    out={}
    for path in paths:
        m=ASSET_RE.search(path.name)
        if not m:
            raise RuntimeError(f"unexpected CHELSA asset basename: {path.name}")
        key=(norm_var(m.group(1)),int(m.group(2)))
        if key in out:
            raise RuntimeError(f"duplicate CHELSA logical asset {key}")
        out[key]=path
    expected={(v,t) for v in VARIABLE_ORDER for t in (-190,20)}
    if set(out)!=expected:
        raise RuntimeError(f"CHELSA logical asset set drift: {sorted(out)}")
    return out


def sample_matrix(coords_lonlat: list[tuple[float,float]], paths: list[Path]) -> np.ndarray:
    import rasterio
    if not coords_lonlat:
        return np.empty((0,len(paths)),float)
    columns=[]
    for path in paths:
        with rasterio.open(path) as ds:
            if ds.crs is None or not ds.crs.is_geographic:
                raise RuntimeError(f"raster is not geographic: {path}")
            vals=np.asarray([float(x[0]) for x in ds.sample(coords_lonlat)],float)
            if ds.nodata is not None:
                vals[np.isclose(vals,float(ds.nodata),rtol=0.0,atol=0.0)]=np.nan
            columns.append(vals)
    return np.column_stack(columns)


def load_geometry(localities_path: Path, edges_path: Path):
    local=defaultdict(dict)
    with localities_path.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            local[row["species"]][int(row["locality_index"])]=(float(row["latitude"]),float(row["longitude"]))
    edges=defaultdict(list)
    with edges_path.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            edges[row["species"]].append((int(row["edge_index"]),int(row["node_left"]),int(row["node_right"])))
    for sp in edges:
        edges[sp].sort()
    return local,edges


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--localities",type=Path,required=True)
    ap.add_argument("--edges",type=Path,required=True)
    ap.add_argument("--host-pairs",type=Path,required=True)
    ap.add_argument("--native-units",type=Path,required=True)
    ap.add_argument("--wgsrpd-support",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--historical-asset",type=Path,action="append",required=True)
    ap.add_argument("--output-edges",type=Path,required=True)
    ap.add_argument("--output-species",type=Path,required=True)
    ap.add_argument("--output-receipt",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!="ttf_historical_host_memory_predictor_materialization_rule_v0.1":
        raise RuntimeError("unexpected host-memory materialization rule")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("host-memory response firewall is open")

    verify_frozen_input_hashes(args.candidates, args.localities, args.edges,
                               args.host_pairs, args.native_units, args.wgsrpd_support)
    candidates=list(csv.DictReader(args.candidates.open(newline="",encoding="utf-8")))
    if len(candidates)!=int(rule["candidate_species"]):
        raise RuntimeError("candidate count drift")
    names=[r["species"] for r in candidates]
    if len(set(names))!=len(names):
        raise RuntimeError("duplicate candidate species")

    with args.host_pairs.open(newline="", encoding="utf-8") as f:
        host_pairs = list(csv.DictReader(f))
    verify_host_pair_identity(candidates, host_pairs)

    assets=load_assets(args.historical_asset)
    present_paths=[assets[(v,20)] for v in VARIABLE_ORDER]
    lgm_paths=[assets[(v,-190)] for v in VARIABLE_ORDER]

    support=list(csv.DictReader(args.wgsrpd_support.open(newline="",encoding="utf-8")))
    support_coords=[(float(r["longitude"]),float(r["latitude"])) for r in support]
    support_current_raw=sample_matrix(support_coords,present_paths)
    finite_support=np.all(np.isfinite(support_current_raw),axis=1)
    if int(np.sum(finite_support))<100:
        raise RuntimeError("too few finite WGSRPD3 support climates")
    reference=support_current_raw[finite_support]
    support_current=np.full_like(support_current_raw,np.nan,dtype=float)
    (support_current_finite,)=standardize_from_reference(reference,support_current_raw[finite_support])
    support_current[finite_support]=support_current_finite
    support_by_unit=defaultdict(list)
    for row,val,keep in zip(support,support_current,finite_support):
        if keep:
            support_by_unit[str(row["area_code_l3"]).strip()].append(np.asarray(val,float))

    native_by_host=defaultdict(set)
    with args.native_units.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            native_by_host[norm_id(row["accepted_plant_name_id"])].add(str(row["area_code_l3"]).strip())

    local,edges=load_geometry(args.localities,args.edges)
    if set(names)-set(local) or set(names)-set(edges):
        raise RuntimeError("candidate species missing from frozen geometry")

    edge_fields=["species","edge_index","M_host","D_host_0BP","D_host_LGM","M_self","D_self_0BP","D_self_LGM"]
    species_fields=["species","edges","host_cloud_points","self_cloud_points","unique_fraction_M_host","predictor_condition_number","status"]
    args.output_edges.parent.mkdir(parents=True,exist_ok=True)
    args.output_species.parent.mkdir(parents=True,exist_ok=True)
    edges_handle=args.output_edges.open("w",newline="",encoding="utf-8")
    species_handle=args.output_species.open("w",newline="",encoding="utf-8")
    ew=csv.DictWriter(edges_handle,fieldnames=edge_fields,lineterminator="\n")
    sw=csv.DictWriter(species_handle,fieldnames=species_fields,lineterminator="\n")
    ew.writeheader(); sw.writeheader()

    complete=0
    status_counts=defaultdict(int)
    for idx,cand in enumerate(candidates):
        sp=cand["species"]
        host_ids=[norm_id(x) for x in cand["accepted_host_ids"].split(";") if norm_id(x)]
        units=sorted(set().union(*(native_by_host.get(h,set()) for h in host_ids)))
        host_vectors=[v for unit in units for v in support_by_unit.get(unit,())]
        if len(host_vectors)<int(rule["host_cloud"]["minimum_Q_host_points"]):
            status_counts["insufficient_host_cloud"]+=1
            sw.writerow({"species":sp,"edges":len(edges[sp]),"host_cloud_points":len(host_vectors),"self_cloud_points":0,"unique_fraction_M_host":"","predictor_condition_number":"","status":"insufficient_host_cloud"})
            continue
        host_cloud=np.vstack(host_vectors)

        loc_rows=local[sp]
        loc_idx=sorted(loc_rows)
        loc_coords=[(loc_rows[i][1],loc_rows[i][0]) for i in loc_idx]
        loc_current_raw=sample_matrix(loc_coords,present_paths)
        finite_loc=np.all(np.isfinite(loc_current_raw),axis=1)
        loc_current=np.full_like(loc_current_raw,np.nan,dtype=float)
        if np.any(finite_loc):
            (loc_current_finite,)=standardize_from_reference(reference,loc_current_raw[finite_loc])
            loc_current[finite_loc]=loc_current_finite
        self_cloud=loc_current[finite_loc]
        if len(self_cloud)<int(rule["insect_self_cloud"]["minimum_Q_self_points"]):
            status_counts["insufficient_self_cloud"]+=1
            sw.writerow({"species":sp,"edges":len(edges[sp]),"host_cloud_points":len(host_cloud),"self_cloud_points":len(self_cloud),"unique_fraction_M_host":"","predictor_condition_number":"","status":"insufficient_self_cloud"})
            continue

        starts=np.array([loc_rows[left] for _,left,_ in edges[sp]],float)
        ends=np.array([loc_rows[right] for _,_,right in edges[sp]],float)
        quad=great_circle_quadrature_latlon(starts,ends,segment_points=int(rule["edge_positions"]["required_finite_positions"]))
        coords=[(float(lon),float(lat)) for lat,lon in quad.reshape(-1,2)]
        cur_raw=sample_matrix(coords,present_paths)
        lgm_raw=sample_matrix(coords,lgm_paths)
        if np.any(~np.isfinite(cur_raw)) or np.any(~np.isfinite(lgm_raw)):
            cur=cur_raw.reshape(len(edges[sp]),quad.shape[1],-1)
            lgm=lgm_raw.reshape(len(edges[sp]),quad.shape[1],-1)
        else:
            cur,lgm=standardize_from_reference(reference,cur_raw,lgm_raw)
            cur=cur.reshape(len(edges[sp]),quad.shape[1],-1)
            lgm=lgm.reshape(len(edges[sp]),quad.shape[1],-1)
        if np.any(~np.isfinite(cur)) or np.any(~np.isfinite(lgm)):
            status_counts["nonfinite_edge_climate"]+=1
            sw.writerow({"species":sp,"edges":len(edges[sp]),"host_cloud_points":len(host_cloud),"self_cloud_points":len(self_cloud),"unique_fraction_M_host":"","predictor_condition_number":"","status":"nonfinite_edge_climate"})
            continue

        m_host,d_host_cur,d_host_lgm=historical_analog_memory(cur,lgm,host_cloud)
        m_self,d_self_cur,d_self_lgm=historical_analog_memory(cur,lgm,self_cloud)
        try:
            unique,condition=predictor_design_diagnostics(m_host,m_self,d_host_cur)
        except ValueError:
            status_counts["rank_deficient_predictor"]+=1
            sw.writerow({"species":sp,"edges":len(edges[sp]),"host_cloud_points":len(host_cloud),"self_cloud_points":len(self_cloud),"unique_fraction_M_host":"","predictor_condition_number":"","status":"rank_deficient_predictor"})
            continue
        for (edge_index,_,_),vals in zip(edges[sp],zip(m_host,d_host_cur,d_host_lgm,m_self,d_self_cur,d_self_lgm)):
            ew.writerow(dict(zip(edge_fields,[sp,edge_index,*[repr(float(v)) for v in vals]])))
        sw.writerow({"species":sp,"edges":len(edges[sp]),"host_cloud_points":len(host_cloud),"self_cloud_points":len(self_cloud),"unique_fraction_M_host":repr(float(unique)),"predictor_condition_number":repr(float(condition)),"status":"complete"})
        status_counts["complete"]+=1; complete+=1
        if (idx+1)%50==0:
            print(json.dumps({"processed":idx+1,"complete":complete},sort_keys=True),flush=True)

    edges_handle.close()
    species_handle.close()
    # Output hashes must be calculated only after both CSV handles are closed.
    payload={
      "schema":"ttf_historical_host_memory_predictor_result_v0.1",
      "status":"RESPONSE_BLIND_PREDICTOR_MATERIALIZED" if complete>=500 else "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_EXTERNAL_DATA",
      "candidate_species":len(candidates),
      "complete_species":complete,
      "status_counts":dict(sorted(status_counts.items())),
      "inputs_sha256":{
        "candidates":sha256_path(args.candidates),"localities":sha256_path(args.localities),"edges":sha256_path(args.edges),"host_pairs":sha256_path(args.host_pairs),"native_units":sha256_path(args.native_units),"wgsrpd_support":sha256_path(args.wgsrpd_support),"rule":sha256_path(args.rule)
      },
      "climate_assets_sha256":{path.name:sha256_path(path) for path in args.historical_asset},
      "outputs_sha256":{"edge_predictors":sha256_path(args.output_edges),
                        "species_diagnostics":sha256_path(args.output_species)},
      "response_firewall":{"sequence_identity_opened":False,"pairwise_genetic_distances_opened":False,"post_IBD_turnover_opened":False,"beta_host_opened":False}
    }
    args.output_receipt.parent.mkdir(parents=True,exist_ok=True)
    args.output_receipt.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps(payload,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())