#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collapse(value: str) -> str:
    return " ".join(str(value).strip().split())


def normalize_locus(value: str) -> str:
    return collapse(re.sub(r"[_\-\s]+", " ", str(value).upper().strip()))


def locus_after_species(raw_gene: str, species: str) -> str:
    gene=normalize_locus(raw_gene)
    sp=normalize_locus(species)
    prefix=sp+" "
    return gene[len(prefix):] if gene.startswith(prefix) else gene


def fasta_headers_only(raw: bytes) -> tuple[str, ...]:
    headers=[]
    for line in raw.splitlines():
        if line.startswith(b">"):
            headers.append(line[1:].decode("utf-8").strip())
    return tuple(headers)


def read_tsv(raw: bytes) -> list[dict[str, str]]:
    text=raw.decode("utf-8")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def strict_terrestrial(rule: dict, *, class_name: str, order: str, family: str) -> bool:
    tax=rule["authoritative_subpanel_filter"]["taxon"]
    if collapse(class_name) != str(tax["class"]):
        return False
    if collapse(order) not in set(map(str, tax["allowed_orders"])):
        return False
    if collapse(order) == "Coleoptera" and collapse(family) in set(
        map(str, tax["excluded_coleoptera_families"])
    ):
        return False
    return True


def palearctic_fraction(realm_labels: list[str], accepted_label: str="palearctic") -> float:
    if not realm_labels:
        return 0.0
    accepted=accepted_label.casefold()
    return sum(str(x).casefold() == accepted for x in realm_labels) / len(realm_labels)


def role_key(namespace: str, archive_sha: str, species: str) -> str:
    return hashlib.sha256(f"{namespace}|{archive_sha}|{species}".encode()).hexdigest()


def assign_roles(species: list[str], archive_sha: str) -> tuple[list[str], list[str]]:
    ordered=sorted(
        map(str, species),
        key=lambda sp: (role_key("palearctic-lgm-role-v0.1", archive_sha, sp), sp),
    )
    cut=len(ordered)//2
    return sorted(ordered[:cut]), sorted(ordered[cut:])


def load_realm_polygons(path: Path, realm_field: str):
    try:
        import geopandas as gpd
    except ImportError as exc:
        raise RuntimeError(
            "geopandas is required for authoritative RESOLVE realm assignment"
        ) from exc
    gdf=gpd.read_file(path)
    field_lookup={str(c).casefold(): str(c) for c in gdf.columns}
    key=str(realm_field).casefold()
    if key not in field_lookup:
        raise RuntimeError(
            f"realm field {realm_field!r} not found; columns={list(gdf.columns)}"
        )
    field=field_lookup[key]
    if gdf.crs is None:
        raise RuntimeError("realm polygons have no CRS")
    return gdf[[field,"geometry"]].rename(columns={field:"realm"}).to_crs(4326)


def assign_realms(points, realms):
    try:
        import geopandas as gpd
        from shapely.geometry import Point
    except ImportError as exc:
        raise RuntimeError("geopandas/shapely are required for realm assignment") from exc

    if not points:
        return []
    frame=gpd.GeoDataFrame(
        {"point_index":list(range(len(points)))},
        geometry=[Point(float(lon), float(lat)) for lat,lon in points],
        crs=4326,
    )
    joined=gpd.sjoin(frame,realms,how="left",predicate="intersects")
    labels=[]
    for i in range(len(points)):
        vals=[
            str(v) for v in joined.loc[joined["point_index"]==i,"realm"].dropna().tolist()
        ]
        if not vals:
            labels.append("")
        elif any(v.casefold()=="palearctic" for v in vals):
            labels.append("Palearctic")
        else:
            labels.append(sorted(vals)[0])
    return labels


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--survivors",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--marker-protocol",type=Path,required=True)
    ap.add_argument("--ecoregions",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!="ttf_genetic_palearctic_lgm_subpanel_rule_v0.1":
        raise RuntimeError("unexpected Palearctic-LGM rule")
    survivors=json.loads(args.survivors.read_text())
    if survivors.get("schema")!="ttf_genetic_phase2_survivor_species_v0.1":
        raise RuntimeError("unexpected survivor manifest")
    protocol=json.loads(args.marker_protocol.read_text())
    aliases={normalize_locus(x) for x in protocol["marker_contract"]["normalized_aliases"]}

    expected=rule["source_domain"]["exact_phylogatr_archive"]
    if args.archive.stat().st_size != int(expected["size_bytes"]):
        raise RuntimeError("exact phylogatR archive size drift")
    archive_sha=sha256_path(args.archive)
    if archive_sha != str(expected["sha256"]):
        raise RuntimeError("exact phylogatR archive SHA256 drift")

    species_order=(
        list(map(str,survivors["species"]["train_species"]))
        + list(map(str,survivors["species"]["eval_species"]))
    )
    if len(species_order)!=211 or len(set(species_order))!=211:
        raise RuntimeError("survivor list is not the frozen 211-species set")

    realms=load_realm_polygons(
        args.ecoregions,
        rule["authoritative_subpanel_filter"]["palearctic"]["polygon_attribute"],
    )

    with zipfile.ZipFile(args.archive) as zf:
        root_candidates=[n for n in zf.namelist() if n.endswith("/genes.txt")]
        if len(root_candidates)!=1:
            raise RuntimeError("archive must contain one phylogatR genes.txt")
        root=root_candidates[0][:-len("genes.txt")]
        genes=read_tsv(zf.read(root+"genes.txt"))
        by_species={}
        for row in genes:
            sp=collapse(row.get("species",""))
            if sp in set(species_order):
                by_species.setdefault(sp,[]).append(row)

        rows=[]
        for species in sorted(species_order):
            candidates=[]
            for row in by_species.get(species,[]):
                if locus_after_species(row.get("gene",""),species) not in aliases:
                    continue
                base=root+str(row["dir"]).rstrip("/")+"/"
                aligned=base+str(row["gene"])+".afa"
                occurrence=base+"occurrences.txt"
                if aligned not in zf.namelist() or occurrence not in zf.namelist():
                    continue
                headers=fasta_headers_only(zf.read(aligned))
                header_set=set(headers)
                occurrence_rows=read_tsv(zf.read(occurrence))
                points=[]
                seen=set()
                for record in occurrence_rows:
                    if str(record.get("phylogatr_id","")) not in header_set:
                        continue
                    try:
                        lat=float(record["latitude"]); lon=float(record["longitude"])
                    except (TypeError,ValueError):
                        continue
                    if not (-90<=lat<=90 and -180<=lon<=180):
                        continue
                    key=(lat,lon)
                    if key not in seen:
                        seen.add(key); points.append(key)
                candidates.append((len(points),len(headers),str(row["gene"]),row,points))

            if not candidates:
                raise RuntimeError(f"no COI-family panel reconstructed for {species}")
            candidates.sort(key=lambda x:(-x[0],-x[1],x[2]))
            nloc,nhead,gene,row,points=candidates[0]

            labels=assign_realms(points,realms)
            fraction=palearctic_fraction(
                labels,
                rule["authoritative_subpanel_filter"]["palearctic"]["accepted_label_casefold"],
            )
            terrestrial=strict_terrestrial(
                rule,
                class_name=row.get("class",""),
                order=row.get("order",""),
                family=row.get("family",""),
            )
            eligible=bool(
                terrestrial
                and fraction >= float(
                    rule["authoritative_subpanel_filter"]["palearctic"][
                        "minimum_fraction_of_selected_localities_in_palearctic"
                    ]
                )
            )
            rows.append({
                "species":species,
                "class":row.get("class",""),
                "order":row.get("order",""),
                "family":row.get("family",""),
                "selected_gene":gene,
                "selected_coi_headers":nhead,
                "selected_unique_localities":nloc,
                "palearctic_localities":sum(x.casefold()=="palearctic" for x in labels),
                "palearctic_fraction":fraction,
                "strict_terrestrial_core":terrestrial,
                "eligible":eligible,
            })

    retained=[row["species"] for row in rows if row["eligible"]]
    source,target=assign_roles(retained,archive_sha)
    filt=rule["authoritative_subpanel_filter"]
    passed=(
        len(retained)>=int(filt["minimum_species"])
        and len(source)>=int(filt["minimum_source_clusters"])
        and len(target)>=int(filt["minimum_target_clusters"])
    )
    status=(
        "PASS_TO_RESPONSE_BLIND_LGM_EXTERNAL_DATA_CENSUS"
        if passed else str(filt["if_failed"])
    )

    role={sp:"source" for sp in source}
    role.update({sp:"target" for sp in target})
    for row in rows:
        row["role"]=role.get(row["species"],"")

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_csv.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    counts={}
    for row in rows:
        if row["eligible"]:
            counts[row["order"]]=counts.get(row["order"],0)+1

    payload={
        "schema":"ttf_genetic_palearctic_lgm_subpanel_census_v0.1",
        "status":status,
        "rule_sha256":sha256_path(args.rule),
        "survivor_manifest_sha256":sha256_path(args.survivors),
        "exact_archive_sha256":archive_sha,
        "ecoregions_source":str(args.ecoregions),
        "ecoregions_file_sha256":sha256_path(args.ecoregions),
        "counts":{
            "frozen_survivors":len(species_order),
            "eligible_species":len(retained),
            "source_clusters":len(source),
            "target_clusters":len(target),
            "eligible_by_order":dict(sorted(counts.items())),
        },
        "retained_species":sorted(retained),
        "source_species":source,
        "target_species":target,
        "next_step":(
            "Acquire frozen independent GBIF occurrences and build current/LGM climatic suitability before any subgroup genetic response."
            if passed else
            "STOP. Do not compute subgroup pairwise genetic transfer response."
        ),
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_subpanel_T_st_computed":False,
            "beta_LGM_computed":False,
        },
    }
    args.output_summary.parent.mkdir(parents=True,exist_ok=True)
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":status,"counts":payload["counts"]},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
