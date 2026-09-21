#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_host_resource import (
    build_insect_host_footprints,
    jaccard_units,
    select_and_split_species,
    species_list_sha256,
)
from ttf.phylogatr_confirmatory import (
    GENES_HEADERS,
    OCCURRENCE_HEADERS,
    collapse_whitespace,
    coordinates_for_headers,
    is_coi_family_locus,
    latlon_to_ecef_km,
)

SOURCE_SHA="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
COMPLETED_DESIGN_SHA="f2f2088c294ecfac3c49659439fd7858330452bd8fbea090513a6b3d590c98f8"
ALIASES=("COI","CO1","COX1","COXI","CYTOCHROME C OXIDASE SUBUNIT I","CYTOCHROME C OXIDASE SUBUNIT 1")


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def quantiles(values):
    x=np.asarray(list(values),float)
    if len(x)==0:
        return None
    q=np.quantile(x,[0,.25,.5,.75,1])
    return dict(zip(("min","q25","median","q75","max"),map(float,q)))


def load_sidecar(interactions:Path,distributions:Path):
    pairs=[]
    with interactions.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        required={"insect_species","accepted_plant_name_id"}
        if not required<=set(reader.fieldnames or ()):
            raise RuntimeError("insect-host sidecar schema drift")
        for row in reader:
            insect=str(row["insect_species"]).strip()
            host=str(row["accepted_plant_name_id"]).strip()
            if insect and host:
                pairs.append((insect,host))
    units=defaultdict(set)
    with distributions.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        required={"accepted_plant_name_id","area_code_l3"}
        if not required<=set(reader.fieldnames or ()):
            raise RuntimeError("host-distribution sidecar schema drift")
        for row in reader:
            host=str(row["accepted_plant_name_id"]).strip()
            unit=str(row["area_code_l3"]).strip()
            if host and unit:
                units[host].add(unit)
    return build_insect_host_footprints(pairs,units)


def _archive_root(z:zipfile.ZipFile)->tuple[str,list[dict[str,str]]]:
    names=set(z.namelist())
    genes_members=[
        name for name in names
        if name.endswith("genes.txt")
        and str(Path(name).parent/"cite.txt").replace("\\","/") in names
    ]
    if len(genes_members)!=1:
        raise RuntimeError("archive root drift")
    genes_member=genes_members[0]
    prefix=str(Path(genes_member).parent).replace("\\","/")
    reader=csv.DictReader(
        io.StringIO(z.read(genes_member).decode("utf-8")),
        delimiter="\t",
    )
    if tuple(reader.fieldnames or ())!=tuple(GENES_HEADERS):
        raise RuntimeError("genes.txt schema drift")
    return prefix,[{k:str(v or "") for k,v in row.items()} for row in reader]


def _headers_only(z:zipfile.ZipFile,member:str)->tuple[str,...]:
    headers=[]
    with z.open(member) as f:
        for raw in f:
            if not raw.startswith(b">"):
                continue
            header=raw[1:].strip().decode("utf-8")
            if header:
                headers.append(header)
    return tuple(headers)


def _occurrences(z:zipfile.ZipFile,member:str)->list[dict[str,str]]:
    reader=csv.DictReader(
        io.StringIO(z.read(member).decode("utf-8")),
        delimiter="\t",
    )
    if tuple(reader.fieldnames or ())!=tuple(OCCURRENCE_HEADERS):
        raise RuntimeError(f"occurrence schema drift: {member}")
    return [{k:str(v or "") for k,v in row.items()} for row in reader]


def scan_geometry(
    archive:Path,
    footprint_species:set[str],
    excluded:set[str],
):
    chosen={}
    family={}
    with zipfile.ZipFile(archive) as z:
        members=set(z.namelist())
        prefix,rows=_archive_root(z)
        candidates=[]
        for row in rows:
            species=collapse_whitespace(row.get("species",""))
            if (
                row.get("kingdom")!="Animalia"
                or row.get("order")!="Lepidoptera"
                or not species
                or species in excluded
                or species not in footprint_species
            ):
                continue
            if not is_coi_family_locus(row.get("gene",""),species,ALIASES):
                continue
            try:
                aligned=int(float(row.get("num_seqs_aligned","0")))
            except ValueError:
                aligned=0
            if aligned<12:
                continue
            candidates.append(row)

        scanned_species=set()
        for row in candidates:
            species=collapse_whitespace(row["species"])
            scanned_species.add(species)
            rel=Path(row["dir"])
            if rel.is_absolute() or ".." in rel.parts:
                raise RuntimeError("unsafe phylogatR dir")
            gene=str(row["gene"])
            fasta=str(Path(prefix)/rel/f"{gene}.afa").replace("\\","/")
            occ=str(Path(prefix)/rel/"occurrences.txt").replace("\\","/")
            if fasta not in members or occ not in members:
                continue
            headers=_headers_only(z,fasta)
            if not headers:
                continue
            latlon=coordinates_for_headers(headers,_occurrences(z,occ))
            if len(latlon)==0:
                continue
            exact=np.unique(latlon,axis=0)
            if len(exact)<12:
                continue
            geometry=prepare_density_scaled_genetic_geometry(
                latlon_to_ecef_km(exact),neighbor_fraction=.15
            )
            if geometry.min_endpoint_disjoint_training_edges<5:
                continue
            key=(-len(exact),-len(headers),gene)
            current=chosen.get(species)
            if current is None or key<current[0]:
                chosen[species]=(key,geometry)
                family[species]=str(row.get("family","")).strip()

    return (
        {species:value[1] for species,value in chosen.items()},
        family,
        len(scanned_species),
    )


@dataclass(frozen=True)
class SupportPools:
    source_pool: dict[str, tuple[str, ...]]
    pair_coverage: dict[str, dict[str, float]]
    eligible_eval_species: tuple[str, ...]
    unsupported_eval_species: tuple[str, ...]


def _edge_midpoints(geometry):
    nodes=np.asarray(geometry.edge_nodes,dtype=np.int64)
    coords=np.asarray(geometry.coordinates,dtype=float)
    return 0.5*(coords[nodes[:,0]]+coords[nodes[:,1]])


def _coverage_numpy(target_mid,source_mid,radius=500.0,chunk=128):
    target_mid=np.asarray(target_mid,float)
    source_mid=np.asarray(source_mid,float)
    within=0
    r2=float(radius)*float(radius)
    for start in range(0,len(target_mid),int(chunk)):
        stop=min(start+int(chunk),len(target_mid))
        delta=target_mid[start:stop,None,:]-source_mid[None,:,:]
        nearest2=np.min(np.sum(delta*delta,axis=2),axis=1)
        within+=int(np.count_nonzero(nearest2<=r2))
    return float(within/len(target_mid))


def prepare_support_only_pools(
    geos,
    train,
    evaluation,
    *,
    support_radius=500.0,
    minimum_target_coverage=0.5,
    minimum_source_species=5,
):
    train=tuple(map(str,train))
    evaluation=tuple(map(str,evaluation))
    mid={name:_edge_midpoints(geos[name]) for name in train+evaluation}
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        cKDTree=None
    trees=(
        {name:cKDTree(mid[name]) for name in train}
        if cKDTree is not None else {}
    )
    source_pool={}
    pair_coverage={}
    eligible=[]
    unsupported=[]
    for target in evaluation:
        supported=[]
        coverage={}
        target_mid=mid[target]
        for source in train:
            if cKDTree is not None:
                distance,_=trees[source].query(
                    target_mid,k=1,distance_upper_bound=float(support_radius)
                )
                value=float(np.mean(np.isfinite(distance)))
            else:
                value=_coverage_numpy(
                    target_mid,mid[source],radius=float(support_radius)
                )
            if value>=float(minimum_target_coverage):
                supported.append(source)
                coverage[source]=value
        source_pool[target]=tuple(supported)
        pair_coverage[target]=coverage
        if len(supported)>=int(minimum_source_species):
            eligible.append(target)
        else:
            unsupported.append(target)
    return SupportPools(
        source_pool=source_pool,
        pair_coverage=pair_coverage,
        eligible_eval_species=tuple(eligible),
        unsupported_eval_species=tuple(unsupported),
    )


def save_design_npz(
    path:Path,
    selected,
    train,
    evaluation,
    geos,
    family,
    footprints,
    diagnostics,
    pools,
):
    names=tuple(selected)
    index={name:i for i,name in enumerate(names)}
    coordinate_offsets=[0]
    coordinate_blocks=[]
    for name in names:
        block=np.asarray(geos[name].coordinates,float)
        coordinate_blocks.append(block)
        coordinate_offsets.append(coordinate_offsets[-1]+len(block))
    coordinates=np.vstack(coordinate_blocks)

    unit_names=tuple(sorted(set().union(*(footprints[name] for name in names))))
    unit_index={unit:i for i,unit in enumerate(unit_names)}
    host_presence=np.zeros((len(names),len(unit_names)),dtype=np.uint8)
    for i,name in enumerate(names):
        for unit in footprints[name]:
            host_presence[i,unit_index[unit]]=1

    target_idx=[]
    source_idx=[]
    coverage=[]
    host_jaccard=[]
    centroid=[]
    edge_ratio=[]
    locality_ratio=[]
    source_counts=[]
    for target in pools.eligible_eval_species:
        tmid=np.asarray(
            0.5*(
                geos[target].coordinates[geos[target].edge_nodes[:,0]]
                +geos[target].coordinates[geos[target].edge_nodes[:,1]]
            ),
            float,
        )
        tc=np.mean(tmid,axis=0)
        sources=pools.source_pool[target]
        source_counts.append(len(sources))
        for source in sources:
            smid=np.asarray(
                0.5*(
                    geos[source].coordinates[geos[source].edge_nodes[:,0]]
                    +geos[source].coordinates[geos[source].edge_nodes[:,1]]
                ),
                float,
            )
            target_idx.append(index[target])
            source_idx.append(index[source])
            coverage.append(float(pools.pair_coverage[target][source]))
            host_jaccard.append(jaccard_units(footprints[target],footprints[source]))
            centroid.append(float(np.linalg.norm(tc-np.mean(smid,axis=0))))
            edge_ratio.append(float(geos[source].n_edges/geos[target].n_edges))
            locality_ratio.append(float(geos[source].n_localities/geos[target].n_localities))

    np.savez_compressed(
        path,
        species_order=np.asarray(names,dtype="U"),
        family=np.asarray([family.get(name,"") for name in names],dtype="U"),
        coordinates=coordinates,
        coordinate_offsets=np.asarray(coordinate_offsets,dtype=np.int64),
        train_indices=np.asarray([index[name] for name in train],dtype=np.int64),
        eval_indices=np.asarray([index[name] for name in evaluation],dtype=np.int64),
        eligible_eval_indices=np.asarray(
            [index[name] for name in pools.eligible_eval_species],dtype=np.int64
        ),
        host_unit_names=np.asarray(unit_names,dtype="U"),
        host_presence=host_presence,
        resolved_host_species=np.asarray(
            [diagnostics[name]["resolved_host_species"] for name in names],
            dtype=np.int64,
        ),
        hosts_with_primary_native_units=np.asarray(
            [diagnostics[name]["hosts_with_primary_native_units"] for name in names],
            dtype=np.int64,
        ),
        target_index=np.asarray(target_idx,dtype=np.int64),
        source_index=np.asarray(source_idx,dtype=np.int64),
        coverage=np.asarray(coverage,float),
        host_resource_jaccard=np.asarray(host_jaccard,float),
        centroid_distance_km=np.asarray(centroid,float),
        edge_count_ratio=np.asarray(edge_ratio,float),
        locality_count_ratio=np.asarray(locality_ratio,float),
        source_count_per_supported_target=np.asarray(source_counts,dtype=np.int64),
    )


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--prior-exclusion-json",type=Path,required=True)
    ap.add_argument("--completed-design-json",type=Path,required=True)
    ap.add_argument("--insect-host-accepted-csv",type=Path,required=True)
    ap.add_argument("--native-wgsrpd3-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    ap.add_argument("--output-design-npz",type=Path,required=True)
    args=ap.parse_args()

    if sha256_path(args.archive)!=SOURCE_SHA:
        raise RuntimeError("source archive SHA drift")
    if sha256_path(args.completed_design_json)!=COMPLETED_DESIGN_SHA:
        raise RuntimeError("completed trait-gradient design SHA drift")

    prior=json.loads(args.prior_exclusion_json.read_text())
    completed=json.loads(args.completed_design_json.read_text())
    excluded=set(map(str,prior["species"]))
    excluded.update(map(str,prior.get("additional_operator_exposure_exclusion",[])))
    excluded.update(map(str,completed["species_order"]))

    footprints,diagnostics=load_sidecar(
        args.insect_host_accepted_csv,args.native_wgsrpd3_csv
    )
    geos,family,scanned_species=scan_geometry(
        args.archive,set(footprints),excluded
    )
    admissible=sorted(set(geos)&set(footprints))
    selected,train,evaluation=select_and_split_species(
        admissible,maximum_species=500
    )
    if len(selected)<6:
        raise RuntimeError("too few species for support census")

    selected_geos={name:geos[name] for name in selected}
    pools=prepare_support_only_pools(
        selected_geos,
        train,
        evaluation,
        support_radius=500.,
        minimum_target_coverage=.5,
        minimum_source_species=5,
    )

    args.output_design_npz.parent.mkdir(parents=True,exist_ok=True)
    save_design_npz(
        args.output_design_npz,
        selected,train,evaluation,selected_geos,family,
        footprints,diagnostics,pools,
    )
    z=np.load(args.output_design_npz,allow_pickle=False)
    jaccard=np.asarray(z["host_resource_jaccard"],float)
    source_counts=np.asarray(z["source_count_per_supported_target"],float)

    fam_counts=Counter(family.get(name,"") or "UNKNOWN" for name in selected)
    species_rows=[]
    eligible=set(pools.eligible_eval_species)
    for name in selected:
        role="train" if name in set(train) else "eval"
        species_rows.append({
            "species":name,
            "family":family.get(name,""),
            "role":role,
            "unique_localities":int(selected_geos[name].n_localities),
            "edge_count":int(selected_geos[name].n_edges),
            **diagnostics[name],
            "geographically_supported_eval":bool(name in eligible),
            "supported_source_species":(
                len(pools.source_pool[name]) if name in set(evaluation) else None
            ),
        })

    gate={
        "minimum_total_species":180,
        "minimum_train_species":80,
        "minimum_eval_species":80,
        "minimum_geographically_supported_eval_targets":90,
        "total_pass":len(selected)>=180,
        "train_pass":len(train)>=80,
        "eval_pass":len(evaluation)>=80,
        "supported_eval_pass":len(pools.eligible_eval_species)>=90,
    }
    gate["all_pass"]=all(
        gate[key] for key in ("total_pass","train_pass","eval_pass","supported_eval_pass")
    )

    out={
        "schema":"ttf_lepidoptera_host_resource_distribution_census_v0.2",
        "status":(
            "PASS_TO_CHARACTER_MASK_GATE"
            if gate["all_pass"]
            else "NOT_EVALUABLE_HOST_RESOURCE_DISTRIBUTION_CENSUS"
        ),
        "source_archive_sha256":SOURCE_SHA,
        "completed_trait_gradient_design_sha256":COMPLETED_DESIGN_SHA,
        "independence":{
            "prior_phase1_species_excluded":len(prior["species"]),
            "completed_trait_gradient_pre_mask_species_excluded":len(completed["species_order"]),
            "operator_exposure_species_excluded":prior.get("additional_operator_exposure_exclusion",[]),
        },
        "host_resource_admissibility":{
            "insect_species_with_at_least_one_primary_native_footprint_before_geometry":len(footprints),
            "aligned_count_ge_12_host_footprint_species_scanned":scanned_species,
            "geometry_and_host_footprint_eligible_species":len(admissible),
        },
        "panel":{
            "selected_species":len(selected),
            "train_species":len(train),
            "eval_species":len(evaluation),
            "geographically_supported_eval_species":len(pools.eligible_eval_species),
            "unsupported_eval_species":len(pools.unsupported_eval_species),
            "selected_species_sha256":species_list_sha256(selected),
            "train_species_sha256":species_list_sha256(train),
            "eval_species_sha256":species_list_sha256(evaluation),
            "design_npz_sha256":sha256_path(args.output_design_npz),
            "source_count_quantiles":quantiles(source_counts),
            "host_resource_jaccard_supported_pair_quantiles":quantiles(jaccard),
            "supported_pair_count":int(len(jaccard)),
            "supported_pairs_zero_host_overlap":int(np.count_nonzero(jaccard==0)),
            "supported_pairs_identical_host_footprint":int(np.count_nonzero(jaccard==1)),
        },
        "selected_trait_free_descriptors":{
            "primary_native_wgsrpd3_unit_count_quantiles":quantiles(
                diagnostics[name]["primary_native_wgsrpd3_units"] for name in selected
            ),
            "resolved_host_species_quantiles":quantiles(
                diagnostics[name]["resolved_host_species"] for name in selected
            ),
            "hosts_with_primary_native_units_quantiles":quantiles(
                diagnostics[name]["hosts_with_primary_native_units"] for name in selected
            ),
            "top_families":dict(fam_counts.most_common(20)),
        },
        "gate":gate,
        "species":species_rows,
        "outcome_firewall":{
            "character_masks_opened":False,
            "sequence_identity_opened":False,
            "pairwise_genetic_distances_opened":False,
            "transfer_statistic_computed":False,
        },
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":out["status"],
        "geometry_host_eligible":len(admissible),
        "selected":len(selected),
        "train":len(train),
        "eval":len(evaluation),
        "supported_eval":len(pools.eligible_eval_species),
        "supported_pairs":len(jaccard),
        "host_jaccard_median":(
            None if len(jaccard)==0 else float(np.median(jaccard))
        ),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
