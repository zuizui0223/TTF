#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path


GENES_HEADERS = (
    "gene","dir","proportion_retained","num_seqs_unaligned","num_seqs_aligned",
    "kingdom","phylum","class","order","family","genus","species","subspecies",
    "different_genbank_species",
)


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize_token(value: str) -> str:
    return re.sub(r"[_\-\s]+"," ",str(value).upper().strip()).strip()


def normalized_locus(raw_gene: str, species: str) -> str:
    gene=normalize_token(raw_gene)
    prefix=normalize_token(species)+" "
    if gene.startswith(prefix):
        gene=gene[len(prefix):]
    return gene


def headers_only(handle) -> list[str]:
    headers=[]
    for line in handle:
        if line.startswith(b">"):
            headers.append(line[1:].decode("utf-8").strip())
    return headers


def matched_unique_latlon(zf: zipfile.ZipFile, row: dict[str,str]) -> tuple[list[tuple[float,float]],int]:
    base="phylogatr-results/"+row["dir"]+"/"
    afa=base+row["gene"]+".afa"
    occ=base+"occurrences.txt"
    headers=headers_only(zf.open(afa))
    occ_rows=list(csv.DictReader(zf.read(occ).decode("utf-8").splitlines(),delimiter="\t"))
    by_id={r["phylogatr_id"]:r for r in occ_rows}
    coords=[]
    matched=0
    for header in headers:
        r=by_id.get(header)
        if r is None:
            continue
        matched+=1
        try:
            lat=float(r["latitude"]); lon=float(r["longitude"])
        except (TypeError,ValueError):
            continue
        if math.isfinite(lat) and math.isfinite(lon):
            coords.append((lat,lon))
    return sorted(set(coords)),matched


def core_fraction(coords: list[tuple[float,float]], rule: dict) -> float:
    if not coords:
        return 0.0
    inside=0
    for lat,lon in coords:
        if (
            float(rule["latitude_min"]) <= lat <= float(rule["latitude_max"])
            and float(rule["longitude_min"]) <= lon <= float(rule["longitude_max"])
        ):
            inside+=1
    return inside/len(coords)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--authorization",type=Path,required=True)
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    contract=json.loads(args.contract.read_text())
    if contract.get("schema")!="ttf_palearctic_lepidoptera_screen_rule_v0.1":
        raise RuntimeError("unexpected screen contract")
    if sha256_path(args.archive)!=contract["source"]["exact_phylogatr_archive_sha256"]:
        raise RuntimeError("exact phylogatR archive SHA256 mismatch")
    if sha256_path(args.authorization)!=contract["source"]["identity_opening_authorization_sha256"]:
        raise RuntimeError("identity-opening authorization SHA256 mismatch")

    authorization=json.loads(args.authorization.read_text())
    species_block=authorization.get("species",{})
    survivors=sorted(
        set(map(str,species_block.get("train_species",[])))
        | set(map(str,species_block.get("eval_species",[])))
    )
    if len(survivors)!=int(contract["source"]["expected_survivors"]):
        raise RuntimeError("survivor count drift")

    aliases={normalize_token(x) for x in contract["panel_reconstruction"]["marker_aliases"]}
    primary_rule=contract["primary_screen"]
    sensitivity_rule=contract["sensitivity_screen"]

    with zipfile.ZipFile(args.archive) as zf:
        genes=list(csv.DictReader(
            zf.read("phylogatr-results/genes.txt").decode("utf-8").splitlines(),
            delimiter="\t",
        ))
        if tuple(genes[0].keys())!=GENES_HEADERS:
            raise RuntimeError("genes.txt schema drift")

        by_species={}
        for row in genes:
            species=str(row["species"])
            if species not in survivors:
                continue
            if normalized_locus(row["gene"],species) not in aliases:
                continue
            coords,matched=matched_unique_latlon(zf,row)
            by_species.setdefault(species,[]).append({
                "row":row,
                "coords":coords,
                "nloc":len(coords),
                "headers":int(row["num_seqs_aligned"]),
                "matched":matched,
            })

    selected={}
    for species in survivors:
        candidates=by_species.get(species,[])
        if not candidates:
            raise RuntimeError(f"no COI-family panel for frozen survivor: {species}")
        candidates.sort(key=lambda x:(-x["nloc"],-x["headers"],x["row"]["gene"]))
        selected[species]=candidates[0]

    rows=[]
    for species,panel in sorted(selected.items()):
        row=panel["row"]
        coords=panel["coords"]
        frac=core_fraction(coords,primary_rule)
        lats=[x[0] for x in coords]; lons=[x[1] for x in coords]
        rows.append({
            "species":species,
            "class":row["class"],
            "order":row["order"],
            "family":row["family"],
            "nloc":len(coords),
            "palearctic_core_fraction":frac,
            "lat_min":min(lats),
            "lat_max":max(lats),
            "lon_min":min(lons),
            "lon_max":max(lons),
        })

    insects=[r for r in rows if r["class"]=="Insecta"]
    lepidoptera=[r for r in insects if r["order"]=="Lepidoptera"]
    primary=[
        r for r in lepidoptera
        if r["palearctic_core_fraction"] >= float(primary_rule["required_fraction_of_selected_localities_inside"])
    ]
    sensitivity=[
        r for r in lepidoptera
        if r["palearctic_core_fraction"] >= float(sensitivity_rule["required_fraction_of_selected_localities_inside"])
    ]
    primary_names={r["species"] for r in primary}
    status=(
        "PASS_TO_LGM_PREDICTOR_CONSTRUCTION"
        if len(primary)>=int(primary_rule["minimum_species_to_continue"])
        else "STOP_PANEL_TOO_SMALL_FOR_LGM_PROGRAM"
    )

    nloc=[int(r["nloc"]) for r in primary]
    nloc_sorted=sorted(nloc)
    median=(
        0.5*(nloc_sorted[len(nloc)//2-1]+nloc_sorted[len(nloc)//2])
        if len(nloc)%2==0 else float(nloc_sorted[len(nloc)//2])
    ) if nloc else None

    payload={
        "schema":"ttf_palearctic_lepidoptera_response_blind_screen_v0.1",
        "status":status,
        "source":{
            "exact_phylogatr_archive_sha256":sha256_path(args.archive),
            "identity_opening_authorization_sha256":sha256_path(args.authorization),
            "authorization_survivors":len(survivors),
        },
        "screen":contract["primary_screen"],
        "sensitivity_screen":contract["sensitivity_screen"],
        "census":{
            "survivors_total":len(rows),
            "insecta":len(insects),
            "lepidoptera":len(lepidoptera),
            "primary_strict_species":len(primary),
            "sensitivity_80pct_species":len(sensitivity),
            "primary_ordered_nonself_dyads":len(primary)*(len(primary)-1),
            "sensitivity_ordered_nonself_dyads":len(sensitivity)*(len(sensitivity)-1),
            "primary_localities_median":median,
            "primary_localities_min":min(nloc) if nloc else None,
            "primary_localities_max":max(nloc) if nloc else None,
            "primary_families":len({r["family"] for r in primary}),
        },
        "primary_species":primary,
        "sensitivity_only_species":[r for r in sensitivity if r["species"] not in primary_names],
        "response_firewall":{
            "sequence_identity_read_for_screen":False,
            "pairwise_genetic_distance_read_for_screen":False,
            "phase4_species_score_read_for_screen":False,
            "pairwise_transfer_response_constructed":False,
        },
        "next_authorized_step":(
            "Construct the frozen response-blind LGM-SDM predictor and run TTF-Q before any subpanel genetic response."
            if status=="PASS_TO_LGM_PREDICTOR_CONSTRUCTION"
            else "STOP; do not construct a subpanel genetic response."
        ),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "survivors":len(rows),
        "insecta":len(insects),
        "lepidoptera":len(lepidoptera),
        "primary":len(primary),
        "sensitivity":len(sensitivity),
        "primary_dyads":len(primary)*(len(primary)-1),
        "genetic_response_used":False,
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
