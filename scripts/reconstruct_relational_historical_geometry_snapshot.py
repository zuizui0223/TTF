#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.phylogatr_confirmatory import canonical_latlon_for_geometry, latlon_to_ecef_km


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--snapshot",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--output-localities",type=Path,required=True)
    ap.add_argument("--output-edges",type=Path,required=True)
    ap.add_argument("--output-receipt",type=Path,required=True)
    args=ap.parse_args()

    contract=json.loads(args.contract.read_text())
    if contract.get("schema")!="ttf_relational_historical_geometry_snapshot_v0.1":
        raise RuntimeError("unexpected historical geometry snapshot contract")
    if contract.get("status")!="FROZEN_RESPONSE_BLIND_GEOMETRY_SNAPSHOT_FROM_EXACT_ARCHIVE":
        raise RuntimeError("historical geometry snapshot is not frozen")
    if any(bool(v) for v in contract["response_firewall"].values()):
        raise RuntimeError("historical geometry snapshot response firewall is open")
    if sha256_path(args.candidates)!=contract["candidate_csv_sha256"]:
        raise RuntimeError("historical candidate CSV SHA drift")

    encoded=b"".join(args.snapshot.read_bytes().split())
    if sha256_bytes(encoded+b"\n")!=contract["snapshot_base64_file_sha256"]:
        raise RuntimeError("historical geometry base64 snapshot SHA drift")
    compressed=base64.b64decode(encoded,validate=True)
    if sha256_bytes(compressed)!=contract["snapshot_gzip_sha256"]:
        raise RuntimeError("historical geometry gzip SHA drift")
    raw=gzip.decompress(compressed)
    if sha256_bytes(raw)!=contract["snapshot_json_sha256"]:
        raise RuntimeError("historical geometry JSON SHA drift")
    locality_map=json.loads(raw.decode("utf-8"))

    candidates=list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    if len(candidates)!=int(contract["species"]):
        raise RuntimeError("historical geometry candidate count drift")
    names=[str(r["species"]) for r in candidates]
    if set(locality_map)!=set(names):
        raise RuntimeError("historical geometry snapshot species universe drift")

    args.output_localities.parent.mkdir(parents=True,exist_ok=True)
    args.output_edges.parent.mkdir(parents=True,exist_ok=True)
    loc_fields=["species","locality_index","latitude","longitude","x_km","y_km","z_km"]
    edge_fields=["species","edge_index","node_left","node_right","mid_x_km","mid_y_km","mid_z_km"]
    nloc=nedge=0
    failures=[]

    with (
        args.output_localities.open("w",newline="",encoding="utf-8") as lh,
        args.output_edges.open("w",newline="",encoding="utf-8") as eh,
    ):
        lw=csv.DictWriter(lh,fieldnames=loc_fields,lineterminator="\n"); lw.writeheader()
        ew=csv.DictWriter(eh,fieldnames=edge_fields,lineterminator="\n"); ew.writeheader()
        for row in candidates:
            species=str(row["species"])
            source=np.asarray(locality_map[species],dtype=float)
            if source.ndim!=2 or source.shape[1]!=2:
                raise RuntimeError(f"{species}: invalid frozen locality matrix")
            if len(source)!=int(row["n_localities"]):
                raise RuntimeError(f"{species}: frozen locality count drift")
            # The snapshot stores canonical lat/lon in the exact geometry order.
            ecef=latlon_to_ecef_km(source)
            geometry=prepare_density_scaled_genetic_geometry(
                ecef,neighbor_fraction=float(contract["neighbor_fraction"])
            )
            canonical=canonical_latlon_for_geometry(source,geometry)
            expected={
                "n_localities":int(row["n_localities"]),
                "graph_k":int(row["graph_k"]),
                "edges":int(row["edges"]),
                "min_endpoint_disjoint_training_edges":int(row["min_endpoint_disjoint_training_edges"]),
            }
            observed={
                "n_localities":int(geometry.n_localities),
                "graph_k":int(geometry.graph_k),
                "edges":int(geometry.n_edges),
                "min_endpoint_disjoint_training_edges":int(geometry.min_endpoint_disjoint_training_edges),
            }
            if observed!=expected:
                failures.append({"species":species,"expected":expected,"observed":observed})
                continue

            xyz=np.asarray(geometry.coordinates,dtype=float)
            for i,(latlon,p) in enumerate(zip(canonical,xyz)):
                lw.writerow({
                    "species":species,"locality_index":i,
                    "latitude":repr(float(latlon[0])),"longitude":repr(float(latlon[1])),
                    "x_km":repr(float(p[0])),"y_km":repr(float(p[1])),"z_km":repr(float(p[2])),
                })
            for i,(left,right) in enumerate(np.asarray(geometry.edge_nodes,dtype=np.int64)):
                midpoint=.5*(xyz[int(left)]+xyz[int(right)])
                ew.writerow({
                    "species":species,"edge_index":i,
                    "node_left":int(left),"node_right":int(right),
                    "mid_x_km":repr(float(midpoint[0])),
                    "mid_y_km":repr(float(midpoint[1])),
                    "mid_z_km":repr(float(midpoint[2])),
                })
            nloc+=int(geometry.n_localities)
            nedge+=int(geometry.n_edges)

    if failures:
        raise RuntimeError(f"historical geometry verification failures: {failures[:3]}")
    if nloc!=int(contract["locality_rows"]) or nedge!=int(contract["edge_rows"]):
        raise RuntimeError(f"historical aggregate geometry drift: localities={nloc}, edges={nedge}")
    locality_sha=sha256_path(args.output_localities)
    edge_sha=sha256_path(args.output_edges)
    if locality_sha!=contract["locality_csv_sha256"]:
        raise RuntimeError("historical reconstructed locality CSV SHA drift")
    if edge_sha!=contract["edge_csv_sha256"]:
        raise RuntimeError("historical reconstructed edge CSV SHA drift")

    receipt={
        "schema":"ttf_relational_historical_geometry_snapshot_reproduction_v0.1",
        "status":"PASS_EXACT_RESPONSE_BLIND_HISTORICAL_GEOMETRY_REPRODUCTION",
        "source_archive_sha256":contract["source_archive_sha256"],
        "candidate_csv_sha256":contract["candidate_csv_sha256"],
        "snapshot_base64_file_sha256":contract["snapshot_base64_file_sha256"],
        "species":len(candidates),
        "locality_rows":nloc,
        "edge_rows":nedge,
        "locality_csv_sha256":locality_sha,
        "edge_csv_sha256":edge_sha,
        "verification_failures":0,
        "response_firewall":{
            "sequence_lines_decoded":False,
            "Study_C_sequence_identity_opened":False,
            "Study_C_pairwise_genetic_distances_opened":False,
            "Study_C_T_st_computed":False,
            "Study_C_beta_hist_computed":False,
        },
    }
    args.output_receipt.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
