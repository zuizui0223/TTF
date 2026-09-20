#!/usr/bin/env python3
from __future__ import annotations

import argparse,csv,hashlib,io,json,stat,tempfile,zipfile
from pathlib import Path

import numpy as np

from ttf.lepidoptera_species_transferability import (
    crossvalidated_species_property_gain,
    decompose_species_transferability,
)
from ttf.lepidoptera_trait_gradient_empirical import (
    frozen_pair_surface,
    pair_transfer_congruence,
    post_ibd_edge_response,
)
from ttf.phylogatr_confirmatory import (
    choose_one_panel_per_species,
    read_genes_rows,
    scan_phylogatr_phase1,
)
from ttf.phylogatr_empirical import extract_species_frozen_edge_distances

ALIASES=("COI","CO1","COX1","COXI","CYTOCHROME C OXIDASE SUBUNIT I","CYTOCHROME C OXIDASE SUBUNIT 1")


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
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
        reader=csv.DictReader(io.StringIO(z.read(genes_member).decode("utf-8")),delimiter="\t")
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


def geometry_covariates(design_npz,pairs:np.ndarray)->np.ndarray:
    target=pairs[:,0].astype(np.int64)
    source=pairs[:,1].astype(np.int64)
    g=np.asarray(pairs[:,3:7],float).copy()
    g[:,1]=np.log1p(g[:,1]/500.0)
    g[:,2]=np.log(g[:,2])
    g[:,3]=np.log(g[:,3])
    kernel=np.asarray(design_npz["geometry_kernel"],float)
    return np.column_stack([g,kernel[target,source]])


def main()->int:
    ap=argparse.ArgumentParser(
        description="Exploratory post-primary decomposition of pairwise Lepidoptera transfer congruence into source exportability and target receptivity."
    )
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--survivor-design-npz",type=Path,required=True)
    ap.add_argument("--parent-primary-result",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    parent=json.loads(args.parent_primary_result.read_text())
    if parent.get("schema")!="ttf_lepidoptera_trait_gradient_empirical_primary_result_v0.1":
        raise RuntimeError("parent primary result schema drift")
    if parent.get("status")!="PRIMARY_ONE_SHOT_DECISION_COMPLETE":
        raise RuntimeError("primary one-shot result is not complete")
    if parent.get("outcome_state",{}).get("new_lepidoptera_empirical_sequence_identity_opened") is not True:
        raise RuntimeError("species-property analysis may run only after the frozen primary opening")
    if sha256_path(args.archive)!=str(parent["source_archive_sha256"]):
        raise RuntimeError("source archive SHA drift")
    if sha256_path(args.survivor_design_npz)!=str(parent["survivor_design_npz_sha256"]):
        raise RuntimeError("survivor design SHA drift")

    z=np.load(args.survivor_design_npz,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    coords=np.asarray(z["coordinates"],float)
    offsets=np.asarray(z["coordinate_offsets"],np.int64)
    expected_coords={
        name:coords[offsets[i]:offsets[i+1]]
        for i,name in enumerate(names)
    }

    with tempfile.TemporaryDirectory(prefix="ttf_lepidoptera_species_property_") as td:
        root=selective_extract(args.archive,Path(td),set(names))
        genes=read_genes_rows(root/"genes.txt")
        universe={str(row.get("species","")).strip() for row in genes}
        scan=scan_phylogatr_phase1(
            root,aliases=ALIASES,excluded_species=universe-set(names),
            min_localities=12,min_endpoint_training_edges=5,neighbor_fraction=.15,
        )
        panels={p.species:p for p in choose_one_panel_per_species(scan.candidates)}
        if set(panels)!=set(names):
            raise RuntimeError("survivor panel reconstruction drift")

        response={}
        for name in names:
            panel=panels[name]
            if not np.array_equal(panel.geometry.coordinates,expected_coords[name]):
                raise RuntimeError(f"survivor geometry drift for {name}")
            distances=extract_species_frozen_edge_distances(
                panel.fasta_path,panel.occurrence_path,
                panel.canonical_latlon,panel.geometry,
                minimum_comparable_fraction=.5,
            )
            response[name]=post_ibd_edge_response(
                distances.genetic_distance,panel.geometry,min_training_edges=5
            )

        pairs,target_index,source_index,_,pair_id=frozen_pair_surface(z)
        pair_score=pair_transfer_congruence(
            response,names,pairs,pair_id,
            np.asarray(z["alignment_target"],dtype=np.int64),
            np.asarray(z["alignment_source"],dtype=np.int64),
        )

    geometry=geometry_covariates(z,pairs)
    source_names=[names[i] for i in source_index]
    target_names=[names[i] for i in target_index]
    decomposition=decompose_species_transferability(
        pair_score,geometry,source_names,target_names
    )
    cv=crossvalidated_species_property_gain(
        pair_score,geometry,source_names,target_names,folds=10
    )

    source_support={name:0 for name in decomposition.source_names}
    target_support={name:0 for name in decomposition.target_names}
    for s,t in zip(source_names,target_names):
        source_support[s]+=1
        target_support[t]+=1

    source_rows=[
        {
            "species":name,
            "exportability":float(effect),
            "posterior_se":float(se),
            "supported_targets":int(source_support[name]),
        }
        for name,effect,se in zip(
            decomposition.source_names,
            decomposition.source_effect,
            decomposition.source_se,
        )
    ]
    target_rows=[
        {
            "species":name,
            "receptivity":float(effect),
            "posterior_se":float(se),
            "supported_sources":int(target_support[name]),
        }
        for name,effect,se in zip(
            decomposition.target_names,
            decomposition.target_effect,
            decomposition.target_se,
        )
    ]

    out={
        "schema":"ttf_lepidoptera_species_transferability_exploratory_v0.1",
        "status":"EXPLORATORY_AFTER_PRIMARY_DECISION",
        "parent_primary_result":str(args.parent_primary_result),
        "parent_primary_decision":parent["decision"],
        "source_archive_sha256":sha256_path(args.archive),
        "survivor_design_npz_sha256":sha256_path(args.survivor_design_npz),
        "pairwise_surface":{
            "pairs":int(len(pair_score)),
            "source_species":len(decomposition.source_names),
            "target_species":len(decomposition.target_names),
            "response":"Pearson congruence of the exact frozen aligned post-IBD edge responses used by the primary study",
            "geometry_covariates":[
                "target_to_source_coverage",
                "log1p_centroid_distance_over_500km",
                "log_edge_count_ratio",
                "log_locality_count_ratio",
                "pairwise_geometry_kernel_similarity",
            ],
        },
        "variance_decomposition":{
            "source_exportability_variance":decomposition.source_variance,
            "target_receptivity_variance":decomposition.target_variance,
            "pair_residual_variance":decomposition.residual_variance,
            "source_fraction":decomposition.source_fraction,
            "target_fraction":decomposition.target_fraction,
            "residual_fraction":decomposition.residual_fraction,
        },
        "pair_cross_validation":{
            "folds":cv.folds,
            "geometry_only_rmse":cv.geometry_only_rmse,
            "species_property_rmse":cv.species_property_rmse,
            "rmse_improvement":cv.rmse_improvement,
            "geometry_only_correlation":cv.geometry_only_correlation,
            "species_property_correlation":cv.species_property_correlation,
        },
        "source_exportability":source_rows,
        "target_receptivity":target_rows,
        "inference":{
            "confirmatory":False,
            "p_value":None,
            "trait_predictors_tested":False,
            "reason":"This species-property estimand was defined after the frozen empirical trait-gradient result. It is exploratory and cannot rewrite the primary decision.",
        },
        "serialized_sequence_identity":False,
        "serialized_edge_genetic_distance_vectors":False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({
        "source_fraction":out["variance_decomposition"]["source_fraction"],
        "target_fraction":out["variance_decomposition"]["target_fraction"],
        "residual_fraction":out["variance_decomposition"]["residual_fraction"],
        "cv_rmse_improvement":out["pair_cross_validation"]["rmse_improvement"],
        "cv_correlation_gain":(
            out["pair_cross_validation"]["species_property_correlation"]
            -out["pair_cross_validation"]["geometry_only_correlation"]
        ),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
