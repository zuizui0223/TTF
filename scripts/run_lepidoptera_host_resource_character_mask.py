#!/usr/bin/env python3
from __future__ import annotations

import argparse,csv,hashlib,io,json,stat,tempfile,zipfile
from pathlib import Path
import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_host_resource_qualification import frozen_alignment_indices
from ttf.phylogatr_character_mask import (
    CharacterMaskError,
    edge_mask_support,
    masks_by_frozen_locality,
    read_canonical_mask_alignment,
)
from ttf.phylogatr_confirmatory import (
    choose_one_panel_per_species,
    read_genes_rows,
    read_occurrence_rows,
    scan_phylogatr_phase1,
)

ALIASES=("COI","CO1","COX1","COXI","CYTOCHROME C OXIDASE SUBUNIT I","CYTOCHROME C OXIDASE SUBUNIT 1")
RULE_SCHEMA="ttf_lepidoptera_host_resource_character_mask_rule_v0.1"


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def selective_extract(archive:Path,dest:Path,species:set[str])->Path:
    with zipfile.ZipFile(archive) as z:
        infos={i.filename:i for i in z.infolist()}
        for info in infos.values():
            p=Path(info.filename)
            if p.is_absolute() or ".." in p.parts:
                raise RuntimeError("unsafe archive member")
            mode=(int(info.external_attr)>>16)&0o170000
            if mode==stat.S_IFLNK:
                raise RuntimeError("symlink archive member forbidden")
        genes_members=[
            n for n in infos
            if n.endswith("genes.txt")
            and str(Path(n).parent/"cite.txt").replace("\\","/") in infos
        ]
        if len(genes_members)!=1:
            raise RuntimeError("archive root drift")
        genes_member=genes_members[0]
        prefix=str(Path(genes_member).parent).replace("\\","/")
        reader=csv.DictReader(io.StringIO(z.read(genes_member).decode("utf-8")),delimiter="\t")
        selected={genes_member,str(Path(prefix)/"cite.txt").replace("\\","/")}
        for row in reader:
            sp=str(row.get("species","")).strip()
            if sp not in species:
                continue
            rel=Path(str(row.get("dir","")))
            if rel.is_absolute() or ".." in rel.parts:
                raise RuntimeError("unsafe genes dir")
            gene=str(row.get("gene",""))
            selected.add(str(Path(prefix)/rel/f"{gene}.afa").replace("\\","/"))
            selected.add(str(Path(prefix)/rel/"occurrences.txt").replace("\\","/"))
        missing=sorted(n for n in selected if n not in infos)
        if missing:
            raise RuntimeError(f"panel source files missing: {missing[:5]}")
        for name in sorted(selected):
            z.extract(infos[name],dest)
    root=dest/Path(genes_member).parent
    return root


def q(values):
    x=np.asarray(values,float)
    z=np.quantile(x,[0,.25,.5,.75,1])
    return dict(zip(("min","q25","median","q75","max"),map(float,z)))


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--design-npz",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    ap.add_argument("--output-survivor-npz",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!=RULE_SCHEMA:
        raise RuntimeError("mask rule schema drift")
    if any(rule["outcome_firewall"].values()):
        raise RuntimeError("mask rule firewall open")
    if sha256_path(args.archive)!=rule["source_archive_sha256"]:
        raise RuntimeError("source archive SHA drift")
    if sha256_path(args.design_npz)!=rule["pre_mask_design_npz_sha256"]:
        raise RuntimeError("pre-mask design SHA drift")

    z=np.load(args.design_npz,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    if len(names)!=500:
        raise RuntimeError("pre-mask species count drift")
    train_idx=np.asarray(z["train_indices"],np.int64)
    eval_idx=np.asarray(z["eval_indices"],np.int64)
    if len(train_idx)!=250 or len(eval_idx)!=250:
        raise RuntimeError("pre-mask split drift")
    coords=np.asarray(z["coordinates"],float)
    offsets=np.asarray(z["coordinate_offsets"],np.int64)
    expected_coords={name:coords[offsets[i]:offsets[i+1]] for i,name in enumerate(names)}

    with tempfile.TemporaryDirectory(prefix="ttf_host_resource_mask_") as td:
        root=selective_extract(args.archive,Path(td),set(names))
        genes=read_genes_rows(root/"genes.txt")
        universe={str(r.get("species","")).strip() for r in genes}
        scan=scan_phylogatr_phase1(
            root,aliases=ALIASES,excluded_species=universe-set(names),
            min_localities=12,min_endpoint_training_edges=5,neighbor_fraction=.15,
        )
        panels={p.species:p for p in choose_one_panel_per_species(scan.candidates)}
        if set(panels)!=set(names):
            raise RuntimeError("pre-mask source-panel reconstruction drift")

        survives=np.zeros(len(names),dtype=bool)
        failed_edges=np.zeros(len(names),dtype=np.int64)
        min_best=np.zeros(len(names),dtype=np.int64)
        parser_failures={}
        for i,name in enumerate(names):
            panel=panels[name]
            if not np.array_equal(panel.geometry.coordinates,expected_coords[name]):
                raise RuntimeError(f"geometry drift for {name}")
            try:
                alignment=read_canonical_mask_alignment(panel.fasta_path)
                occ=read_occurrence_rows(panel.occurrence_path)
                grouped=masks_by_frozen_locality(alignment,occ,panel.canonical_latlon)
                support=edge_mask_support(
                    grouped,panel.geometry.edge_nodes,
                    alignment_length=alignment.alignment_length,
                    minimum_comparable_fraction=.5,
                )
            except CharacterMaskError as exc:
                survives[i]=False
                failed_edges[i]=int(panel.geometry.n_edges)
                min_best[i]=0
                parser_failures[name]=type(exc).__name__+":"+str(exc)
                continue
            survives[i]=support.all_edges_valid
            failed_edges[i]=int(np.count_nonzero(~support.valid_edges))
            min_best[i]=int(np.min(support.best_comparable_columns))

    survivor_idx=np.flatnonzero(survives)
    failed_idx=np.flatnonzero(~survives)
    train_survivor=np.asarray([i for i in train_idx if survives[i]],np.int64)
    eval_survivor=np.asarray([i for i in eval_idx if survives[i]],np.int64)

    # Recompute unchanged 500-km support from survivor training species only.
    geos={}
    midpoint={}
    for i in survivor_idx:
        name=names[int(i)]
        g=prepare_density_scaled_genetic_geometry(expected_coords[name],neighbor_fraction=.15)
        geos[int(i)]=g
        nodes=np.asarray(g.edge_nodes,np.int64)
        midpoint[int(i)]=0.5*(g.coordinates[nodes[:,0]]+g.coordinates[nodes[:,1]])
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        cKDTree=None

    supported_eval=[]
    source_counts=[]
    support_pairs=[]
    for ti in eval_survivor:
        tmid=midpoint[int(ti)]
        sources=[]
        for si in train_survivor:
            smid=midpoint[int(si)]
            if cKDTree is not None:
                dist,_=cKDTree(smid).query(tmid,k=1)
            else:
                d=tmid[:,None,:]-smid[None,:,:]
                dist=np.sqrt(np.min(np.sum(d*d,axis=2),axis=1))
            coverage=float(np.mean(dist<=500.0))
            if coverage>=.5:
                sources.append((int(si),coverage))
        if len(sources)>=5:
            supported_eval.append(int(ti))
            source_counts.append(len(sources))
            support_pairs.extend((int(ti),si,cov) for si,cov in sources)

    gate=rule["gate"]
    passed=(
        len(survivor_idx)>=int(gate["minimum_total_survivors"])
        and len(train_survivor)>=int(gate["minimum_train_survivors"])
        and len(eval_survivor)>=int(gate["minimum_eval_survivors"])
        and len(supported_eval)>=int(gate["minimum_geographically_supported_eval_targets"])
    )

    # Subset frozen response-blind host-resource objects exactly.
    old_to_new={int(old):i for i,old in enumerate(survivor_idx)}
    new_train=np.asarray([old_to_new[int(i)] for i in train_survivor],np.int64)
    new_eval=np.asarray([old_to_new[int(i)] for i in eval_survivor],np.int64)
    new_supported=np.asarray([old_to_new[int(i)] for i in supported_eval],np.int64)
    survivor_names=np.asarray([names[int(i)] for i in survivor_idx],dtype="U")
    survivor_coords=np.concatenate([expected_coords[names[int(i)]] for i in survivor_idx],axis=0)
    lengths=np.asarray([len(expected_coords[names[int(i)]]) for i in survivor_idx],np.int64)
    survivor_offsets=np.r_[0,np.cumsum(lengths)]

    support_target=np.asarray([old_to_new[t] for t,_,_ in support_pairs],np.int64)
    support_source=np.asarray([old_to_new[s] for _,s,_ in support_pairs],np.int64)
    support_coverage=np.asarray([cov for _,_,cov in support_pairs],float)

    args.output_survivor_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_survivor_npz,
        species_order=survivor_names,
        family=np.asarray(z["family"],dtype="U")[survivor_idx],
        coordinates=survivor_coords,
        coordinate_offsets=survivor_offsets,
        train_indices=new_train,
        eval_indices=new_eval,
        eligible_eval_indices=new_supported,
        host_unit_names=np.asarray(z["host_unit_names"],dtype="U"),
        host_presence=np.asarray(z["host_presence"],np.uint8)[survivor_idx],
        hosts_with_primary_native_units=np.asarray(z["hosts_with_primary_native_units"],np.int64)[survivor_idx],
        support_target_index=support_target,
        support_source_index=support_source,
        support_coverage=support_coverage,
    )

    out={
        "schema":"ttf_lepidoptera_host_resource_character_mask_census_v0.1",
        "status":"PASS_TO_EXACT_SURVIVOR_QUALIFICATION" if passed else "NOT_EVALUABLE_CHARACTER_SUPPORT",
        "executed_on":"2026-09-21",
        "parent_rule":str(args.rule),
        "source_archive_sha256":sha256_path(args.archive),
        "pre_mask_design_npz_sha256":sha256_path(args.design_npz),
        "post_mask":{
            "survivors":int(len(survivor_idx)),
            "failed_character_support":int(len(failed_idx)),
            "train_survivors":int(len(train_survivor)),
            "eval_survivors":int(len(eval_survivor)),
            "geographically_supported_eval_targets":int(len(supported_eval)),
            "supported_target_source_pairs":int(len(support_pairs)),
            "source_count_quantiles":q(source_counts) if source_counts else None,
        },
        "survivor_species_sha256":hashlib.sha256(("\n".join(sorted(map(str,survivor_names)))+"\n").encode()).hexdigest(),
        "failed_species_sha256":hashlib.sha256(("\n".join(sorted(names[int(i)] for i in failed_idx))+"\n").encode()).hexdigest(),
        "survivor_design_npz_sha256":sha256_path(args.output_survivor_npz),
        "character_support":{
            "parser_structural_failures":int(len(parser_failures)),
            "parser_structural_failure_species":sorted(parser_failures),
            "failed_edges_total":int(np.sum(failed_edges)),
            "failed_edges_per_failed_species_quantiles":q(failed_edges[failed_idx]) if len(failed_idx) else None,
            "minimum_best_comparable_columns_quantiles":q(min_best),
            "nucleotide_identity_persisted":False,
            "pairwise_nucleotide_differences_computed":False,
        },
        "gate":{
            "minimum_total_survivors":int(gate["minimum_total_survivors"]),
            "minimum_train_survivors":int(gate["minimum_train_survivors"]),
            "minimum_eval_survivors":int(gate["minimum_eval_survivors"]),
            "minimum_geographically_supported_eval_targets":int(gate["minimum_geographically_supported_eval_targets"]),
            "all_pass":bool(passed),
        },
        "outcome_firewall":{
            "sequence_identity_opened":False,
            "pairwise_genetic_distances_opened":False,
            "transfer_statistic_computed":False,
        }
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],**out["post_mask"],"survivor_design_sha":out["survivor_design_npz_sha256"]},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
