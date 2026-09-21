#!/usr/bin/env python3
from __future__ import annotations

import argparse,csv,hashlib,io,json,stat,tempfile,zipfile
from pathlib import Path
import numpy as np

from ttf.conditional_transfer import (
    prepare_cached_target_conditioned_transfer,
    prepare_fully_cached_target_conditioned_transfer,
    prepare_target_source_pools,
    score_fully_cached_target_conditioned_batch,
)
from ttf.geometry_control import length_orthogonalized_turnover
import ttf.lepidoptera_trait_gradient_empirical as empirical_module
from ttf.lepidoptera_trait_gradient_empirical import (
    post_ibd_edge_response,
    score_empirical_trait_gradient,
)
from ttf.phylogatr_character_mask import (
    edge_mask_support,
    masks_by_frozen_locality,
    read_canonical_mask_alignment,
)
from ttf.phylogatr_compact_execution import prepare_phylogatr_compact_ttf_design
from ttf.phylogatr_confirmatory import (
    choose_one_panel_per_species,
    read_genes_rows,
    read_occurrence_rows,
    scan_phylogatr_phase1,
)
from ttf.phylogatr_empirical import extract_species_frozen_edge_distances

ALIASES=("COI","CO1","COX1","COXI","CYTOCHROME C OXIDASE SUBUNIT I","CYTOCHROME C OXIDASE SUBUNIT 1")
RULE_SCHEMA="ttf_lepidoptera_trait_gradient_empirical_rule_v0.1"

def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""): h.update(chunk)
    return h.hexdigest()

def selective_extract(archive:Path,dest:Path,survivors:set[str])->Path:
    with zipfile.ZipFile(archive) as z:
        infos={info.filename:info for info in z.infolist()}
        for info in infos.values():
            p=Path(info.filename)
            if p.is_absolute() or ".." in p.parts:
                raise RuntimeError("unsafe archive member")
            mode=(int(info.external_attr)>>16)&0o170000
            if mode==stat.S_IFLNK:
                raise RuntimeError("symlink archive member forbidden")
        genes_members=[
            name for name in infos
            if name.endswith("genes.txt")
            and str(Path(name).parent/"cite.txt").replace("\\","/") in infos
        ]
        if len(genes_members)!=1:
            raise RuntimeError("archive root drift")
        genes_member=genes_members[0]
        prefix=str(Path(genes_member).parent).replace("\\","/")
        genes_bytes=z.read(genes_member)
        reader=csv.DictReader(io.StringIO(genes_bytes.decode("utf-8")),delimiter="\t")
        selected={genes_member,str(Path(prefix)/"cite.txt").replace("\\","/")}
        for row in reader:
            species=str(row.get("species","")).strip()
            if species not in survivors:
                continue
            rel=Path(str(row.get("dir","")))
            if rel.is_absolute() or ".." in rel.parts:
                raise RuntimeError("unsafe genes.txt dir")
            gene=str(row.get("gene",""))
            selected.add(str(Path(prefix)/rel/f"{gene}.afa").replace("\\","/"))
            selected.add(str(Path(prefix)/rel/"occurrences.txt").replace("\\","/"))
        missing=sorted(name for name in selected if name not in infos)
        if missing:
            raise RuntimeError(f"survivor source files missing: {missing[:5]}")
        for name in sorted(selected):
            z.extract(infos[name],dest)
    root=dest/Path(genes_member).parent
    if not (root/"genes.txt").is_file() or not (root/"cite.txt").is_file():
        raise RuntimeError("selective archive extraction drift")
    return root

def load_design(path:Path,expected_sha:str):
    if sha256_path(path)!=expected_sha: raise RuntimeError("survivor design SHA drift")
    z=np.load(path,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    coords=np.asarray(z["coordinates"],float); offsets=np.asarray(z["coordinate_offsets"],np.int64)
    expected={n:coords[offsets[i]:offsets[i+1]] for i,n in enumerate(names)}
    return z,names,expected

def baseline_score(panels,response,z):
    names=tuple(map(str,z["species_order"]))
    train=tuple(names[i] for i in np.asarray(z["train_indices"],int))
    evaluation=tuple(names[i] for i in np.asarray(z["eval_indices"],int))
    geos={n:panels[n].geometry for n in names}
    compact=prepare_phylogatr_compact_ttf_design(
        geos,train_species=train,eval_species=evaluation,
        bandwidth=500.,prior_strength=.25,segment_points=5,min_training_edges=5,
    )
    tm={n:compact.template_edges[n].midpoint for n in compact.train_species}
    em={n:compact.template_edges[n].midpoint for n in compact.eval_species}
    pools=prepare_target_source_pools(
        compact.prepared,tm,em,support_radius=500.,
        minimum_target_coverage=.5,minimum_source_species=5,
    )
    expected_eligible=tuple(names[i] for i in np.asarray(z["eligible_eval_indices"],int))
    if tuple(pools.eligible_eval_species)!=expected_eligible:
        raise RuntimeError("survivor geographic support drift")
    cache=prepare_cached_target_conditioned_transfer(compact.prepared,pools)
    full=prepare_fully_cached_target_conditioned_transfer(cache,eval_species=expected_eligible)
    train_response={
        n:length_orthogonalized_turnover(
            response[n],compact.template_edges[n].length
        )[:,None]
        for n in compact.train_species
    }
    eval_response={n:response[n][:,None] for n in compact.eval_species}
    scored=score_fully_cached_target_conditioned_batch(full,train_response,eval_response)
    return float(scored.statistics[0]),{
        n:float(scored.species_scores[n][0]) for n in expected_eligible
    }

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--survivor-design-npz",type=Path,required=True)
    ap.add_argument("--private-reference-npy",type=Path,required=True)
    ap.add_argument("--geometry-reference-npy",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--survivor-requalification",type=Path,required=True)
    ap.add_argument("--authorization",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    authorization=json.loads(args.authorization.read_text())
    if authorization.get("schema")!="ttf_lepidoptera_trait_gradient_empirical_authorization_v0.1" or authorization.get("status")!="AUTHORIZE_ONE_SHOT_LEPIDOPTERA_IDENTITY_OPENING":
        raise RuntimeError("empirical opening authorization missing")
    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!=RULE_SCHEMA or rule.get("status")!="FROZEN_BEFORE_NEW_LEPIDOPTERA_NUCLEOTIDE_IDENTITY_OPENING":
        raise RuntimeError("empirical rule is not frozen")
    if any(rule["outcome_firewall"].values()): raise RuntimeError("empirical rule firewall open")
    frozen=authorization["frozen_sha256"]
    checks={
        "empirical_rule":sha256_path(args.rule),
        "empirical_scorer":sha256_path(Path(empirical_module.__file__)),
        "empirical_runner":sha256_path(Path(__file__)),
        "survivor_requalification":sha256_path(args.survivor_requalification),
        "source_archive":sha256_path(args.archive),
        "survivor_design_npz":sha256_path(args.survivor_design_npz),
        "private_reference_npy":sha256_path(args.private_reference_npy),
        "geometry_reference_npy":sha256_path(args.geometry_reference_npy),
    }
    if checks!=frozen:
        drift={key:(frozen.get(key),value) for key,value in checks.items() if frozen.get(key)!=value}
        raise RuntimeError(f"empirical authorization hash drift: {drift}")
    rq=json.loads(args.survivor_requalification.read_text())
    if rq.get("status")!="PASS" or rq.get("empirical_identity_opening_authorized") is not True:
        raise RuntimeError("survivor requalification has not authorized opening")
    if sha256_path(args.archive)!=rule["source_archive"]["sha256"]:
        raise RuntimeError("source archive SHA drift")
    if sha256_path(args.private_reference_npy)!=rule["survivor_references"]["private_npy_sha256"]:
        raise RuntimeError("private reference SHA drift")
    if sha256_path(args.geometry_reference_npy)!=rule["survivor_references"]["geometry_confounded_npy_sha256"]:
        raise RuntimeError("geometry reference SHA drift")
    z,names,expected_coords=load_design(
        args.survivor_design_npz,rule["survivor_design"]["design_npz_sha256"]
    )
    if len(names)!=240: raise RuntimeError("survivor species count drift")

    # Everything before the marked transition below is response blind or mask only.
    with tempfile.TemporaryDirectory(prefix="ttf_lepidoptera_empirical_") as td:
        rootdir=Path(td)
        root=selective_extract(args.archive,rootdir,set(names))
        genes=read_genes_rows(root/"genes.txt")
        universe={str(row.get("species","")).strip() for row in genes}
        excluded=universe-set(names)
        scan=scan_phylogatr_phase1(
            root,aliases=ALIASES,excluded_species=excluded,
            min_localities=12,min_endpoint_training_edges=5,neighbor_fraction=.15,
        )
        selected=choose_one_panel_per_species(scan.candidates)
        panels={p.species:p for p in selected}
        if set(panels)!=set(names): raise RuntimeError("survivor source-panel reconstruction drift")
        for name in names:
            panel=panels[name]
            if not np.array_equal(panel.geometry.coordinates,expected_coords[name]):
                raise RuntimeError(f"survivor geometry drift for {name}")
            mask=read_canonical_mask_alignment(panel.fasta_path)
            occ=read_occurrence_rows(panel.occurrence_path)
            grouped=masks_by_frozen_locality(mask,occ,panel.canonical_latlon)
            support=edge_mask_support(
                grouped,panel.geometry.edge_nodes,
                alignment_length=mask.alignment_length,
                minimum_comparable_fraction=.5,
            )
            if not support.all_edges_valid:
                raise RuntimeError(f"survivor mask admissibility drift for {name}")

        # AUTHORIZED IDENTITY OPENING: all frozen integrity checks have passed.
        response={}
        source_checks={}
        for name in names:
            panel=panels[name]
            distances=extract_species_frozen_edge_distances(
                panel.fasta_path,panel.occurrence_path,
                panel.canonical_latlon,panel.geometry,
                minimum_comparable_fraction=.5,
            )
            response[name]=post_ibd_edge_response(
                distances.genetic_distance,panel.geometry,min_training_edges=5
            )
            source_checks[name]={
                "frozen_edges":int(panel.geometry.n_edges),
                "minimum_valid_sequence_pairs_per_edge":int(np.min(distances.valid_pair_counts)),
            }

        private=np.load(args.private_reference_npy,allow_pickle=False)
        geometry=np.load(args.geometry_reference_npy,allow_pickle=False)
        primary=score_empirical_trait_gradient(
            response,z,private_reference=private,
            geometry_reference=geometry,alpha=.05,
        )
        # Frozen primary decision is complete before the descriptive baseline.
        baseline,baseline_species=baseline_score(panels,response,z)

    out={
        "schema":"ttf_lepidoptera_trait_gradient_empirical_result_v0.1",
        "status":"EMPIRICAL_RESULT_OPENED_UNDER_FROZEN_RULE",
        "primary_trait_gradient":{
            "statistic":primary.statistic,
            "envelope_p_value":primary.p_value,
            "alpha":0.05,
            "positive":primary.positive,
            "finite_target_correlations":len(primary.per_target),
            "per_target_correlations":{str(k):float(v) for k,v in sorted(primary.per_target.items())},
        },
        "secondary_within_lepidoptera_baseline":{
            "role":"DESCRIPTIVE_ONLY",
            "statistic":baseline,
            "species_scores":baseline_species,
            "inferential_claim_authorized":False,
        },
        "source_integrity":{
            "archive_sha256_verified":True,
            "survivor_design_sha256_verified":True,
            "survivor_mask_admissibility_reverified":True,
            "species":source_checks,
        },
        "decision":"POSITIVE_TRAIT_SIMILARITY_GRADIENT" if primary.positive else "NO_DETECTED_POSITIVE_TRAIT_SIMILARITY_GRADIENT",
        "new_lepidoptera_empirical_sequence_identity_opened":True,
        "new_lepidoptera_empirical_pairwise_genetic_distances_opened":True,
        "new_lepidoptera_empirical_transfer_statistic_computed":True,
        "serialized_sequence_identity":False,
        "serialized_edge_genetic_distance_vectors":False,
        "claim_boundary":rule["one_shot_decision"][
            "interpretation_positive" if primary.positive else "interpretation_negative"
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({
        "decision":out["decision"],
        "primary_statistic":primary.statistic,
        "primary_p_value":primary.p_value,
        "descriptive_baseline":baseline,
    },sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
