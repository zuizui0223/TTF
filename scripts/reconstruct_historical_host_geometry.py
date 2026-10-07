#!/usr/bin/env python3
"""Reconstruct exact response-blind geometry for the historical-host panel."""
from __future__ import annotations

import argparse,csv,hashlib,json
from pathlib import Path

try:
    from scripts.reconstruct_relational_environment_geometry import (
        edge_rows,locality_rows,reconstruct_candidate,verify_candidate,sha256_path,
    )
except ModuleNotFoundError:
    from reconstruct_relational_environment_geometry import (
        edge_rows,locality_rows,reconstruct_candidate,verify_candidate,sha256_path,
    )

EXPECTED_CANDIDATE_SHA="1c8b35a7a005af639bbabc4d730fbe723df134c8371d04019a564740b095dbef"
EXPECTED_SOURCE_SHA="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
EXPECTED_SPECIES_DIGEST="6bec34a07c98c19a7b8c90bb89333fa7dde4212f899573fae95639d72163f558"

def digest_names(names):
    return hashlib.sha256(("\n".join(sorted(names))+"\n").encode()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--panel",type=Path,required=True)
    ap.add_argument("--source-sha256",required=True)
    ap.add_argument("--output-localities",type=Path,required=True)
    ap.add_argument("--output-edges",type=Path,required=True)
    ap.add_argument("--output-receipt",type=Path,required=True)
    args=ap.parse_args()

    if args.source_sha256!=EXPECTED_SOURCE_SHA:
        raise RuntimeError("historical-host source archive SHA drift")
    if sha256_path(args.candidates)!=EXPECTED_CANDIDATE_SHA:
        raise RuntimeError("historical-host precensus candidate SHA drift")

    candidates=list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    by_species={str(r["species"]).strip():r for r in candidates}
    if len(by_species)!=600:
        raise RuntimeError("historical-host precensus must contain 600 unique species")
    panel=list(csv.DictReader(args.panel.open(encoding="utf-8")))
    species=[str(r["species"]).strip() for r in panel]
    if len(species)!=140 or len(set(species))!=140:
        raise RuntimeError("historical-host panel must contain 140 unique species")
    if digest_names(species)!=EXPECTED_SPECIES_DIGEST:
        raise RuntimeError("historical-host panel species digest drift")
    missing=sorted(set(species)-set(by_species))
    if missing:
        raise RuntimeError(f"panel species missing from frozen precensus: {missing}")

    root=args.root.resolve()
    if not (root/"genes.txt").is_file() or not (root/"cite.txt").is_file():
        raise FileNotFoundError("root must contain genes.txt and cite.txt")

    lf=["species","locality_index","latitude","longitude","x_km","y_km","z_km"]
    ef=["species","edge_index","node_left","node_right","mid_x_km","mid_y_km","mid_z_km"]
    args.output_localities.parent.mkdir(parents=True,exist_ok=True)
    args.output_edges.parent.mkdir(parents=True,exist_ok=True)
    nloc=nedge=0
    with args.output_localities.open("w",newline="",encoding="utf-8") as lfh, args.output_edges.open("w",newline="",encoding="utf-8") as efh:
        lw=csv.DictWriter(lfh,fieldnames=lf,lineterminator="\n")
        ew=csv.DictWriter(efh,fieldnames=ef,lineterminator="\n")
        lw.writeheader(); ew.writeheader()
        for sp in species:
            row=by_species[sp]
            geometry,canonical_latlon=reconstruct_candidate(root,row,neighbor_fraction=0.15)
            verify_candidate(row,geometry)
            lr=list(locality_rows(sp,geometry,canonical_latlon))
            er=list(edge_rows(sp,geometry))
            lw.writerows(lr); ew.writerows(er)
            nloc+=len(lr); nedge+=len(er)

    receipt={
        "schema":"ttf_genetic_historical_host_geometry_reconstruction_v0.1",
        "status":"PASS_RESPONSE_BLIND_HISTORICAL_HOST_GEOMETRY_RECONSTRUCTION",
        "source_archive_sha256":args.source_sha256,
        "candidate_csv_sha256":sha256_path(args.candidates),
        "panel_species_digest_sha256":digest_names(species),
        "species":len(species),
        "neighbor_fraction":0.15,
        "locality_rows":nloc,
        "edge_rows":nedge,
        "locality_csv_sha256":sha256_path(args.output_localities),
        "edge_csv_sha256":sha256_path(args.output_edges),
        "verification_failures":0,
        "response_firewall":{
            "sequence_lines_decoded":False,
            "fresh_sequence_identity_opened":False,
            "fresh_pairwise_genetic_distances_opened":False,
            "fresh_post_ibd_turnover_computed":False,
            "fresh_beta_host_LGM_computed":False,
        },
    }
    args.output_receipt.parent.mkdir(parents=True,exist_ok=True)
    args.output_receipt.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
