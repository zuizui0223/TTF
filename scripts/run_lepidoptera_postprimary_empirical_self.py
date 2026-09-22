#!/usr/bin/env python3
from __future__ import annotations

import argparse,csv,hashlib,io,json,stat,tempfile,zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_self_detectability import prepare_genetic_self_detectability
from ttf.phylogatr_character_mask import (
    edge_mask_support,
    masks_by_frozen_locality,
    read_canonical_mask_alignment,
)
from ttf.phylogatr_compact_ibd import prepare_compact_crossfit_ibd_design
from ttf.phylogatr_compact_self_detectability import (
    score_phylogatr_compact_self_world_batch,
)
from ttf.phylogatr_confirmatory import (
    choose_one_panel_per_species,
    read_genes_rows,
    read_occurrence_rows,
    scan_phylogatr_phase1,
)
from ttf.phylogatr_empirical import extract_species_frozen_edge_distances
from ttf.private_null_inference import upper_monte_carlo_pvalue
from ttf.core import SpeciesEdges

ALIASES=("COI","CO1","COX1","COXI","CYTOCHROME C OXIDASE SUBUNIT I","CYTOCHROME C OXIDASE SUBUNIT 1")
RULE_SCHEMA="ttf_lepidoptera_postprimary_self_detectability_rule_v0.1"
SOURCE_ARCHIVE_SHA="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def template_edges(name:str,geometry):
    nodes=np.asarray(geometry.edge_nodes,np.int64)
    start=geometry.coordinates[nodes[:,0]]
    end=geometry.coordinates[nodes[:,1]]
    return SpeciesEdges(
        species=name,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=.5*(start+end),
        length=np.linalg.norm(end-start,axis=1),
        turnover=np.zeros(len(nodes),float),
    )


def selective_extract(archive:Path,dest:Path,species:set[str])->Path:
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
        reader=csv.DictReader(
            io.StringIO(z.read(genes_member).decode("utf-8")),
            delimiter="\t",
        )
        selected={genes_member,str(Path(prefix)/"cite.txt").replace("\\","/")}
        for row in reader:
            name=str(row.get("species","")).strip()
            if name not in species:
                continue
            rel=Path(str(row.get("dir","")))
            if rel.is_absolute() or ".." in rel.parts:
                raise RuntimeError("unsafe genes.txt dir")
            gene=str(row.get("gene",""))
            selected.add(str(Path(prefix)/rel/f"{gene}.afa").replace("\\","/"))
            selected.add(str(Path(prefix)/rel/"occurrences.txt").replace("\\","/"))
        missing=sorted(name for name in selected if name not in infos)
        if missing:
            raise RuntimeError(f"source files missing: {missing[:5]}")
        for name in sorted(selected):
            z.extract(infos[name],dest)
    root=dest/Path(genes_member).parent
    if not (root/"genes.txt").is_file():
        raise RuntimeError("selective extraction drift")
    return root


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--design-npz",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--qualification",type=Path,required=True)
    ap.add_argument("--reference-npy",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!=RULE_SCHEMA:
        raise RuntimeError("post-primary self rule schema drift")
    qual=json.loads(args.qualification.read_text())
    if qual.get("schema")!="ttf_lepidoptera_postprimary_self_qualification_v0.1":
        raise RuntimeError("self qualification schema drift")
    if qual.get("panel")!="host" or qual.get("status")!="PASS" or qual.get("empirical_self_eligible") is not True:
        raise RuntimeError("host self qualification has not passed")
    panel=rule["panels"]["host_resource_panel"]
    if sha256_path(args.design_npz)!=panel["survivor_design_sha256"]:
        raise RuntimeError("host survivor design SHA drift")
    if sha256_path(args.archive)!=SOURCE_ARCHIVE_SHA:
        raise RuntimeError("source archive SHA drift")
    reference=np.load(args.reference_npy,allow_pickle=False)
    if reference.shape!=(1000,):
        raise RuntimeError("host self reference length drift")

    z=np.load(args.design_npz,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    coords=np.asarray(z["coordinates"],float)
    offsets=np.asarray(z["coordinate_offsets"],np.int64)
    eval_idx=np.asarray(z["eval_indices"],np.int64)
    eval_names=tuple(names[int(i)] for i in eval_idx)
    if len(names)!=int(panel["species"]) or len(eval_names)!=int(panel["eval_species"]):
        raise RuntimeError("host self panel count drift")
    expected={
        name:coords[offsets[i]:offsets[i+1]]
        for i,name in enumerate(names)
        if name in set(eval_names)
    }

    with tempfile.TemporaryDirectory(prefix="ttf_postprimary_host_self_") as td:
        root=selective_extract(args.archive,Path(td),set(eval_names))
        genes=read_genes_rows(root/"genes.txt")
        universe={str(row.get("species","")).strip() for row in genes}
        scan=scan_phylogatr_phase1(
            root,
            aliases=ALIASES,
            excluded_species=universe-set(eval_names),
            min_localities=12,
            min_endpoint_training_edges=5,
            neighbor_fraction=.15,
        )
        panels={p.species:p for p in choose_one_panel_per_species(scan.candidates)}
        if set(panels)!=set(eval_names):
            raise RuntimeError("host self source-panel reconstruction drift")

        geos={}
        templates={}
        ibd={}
        distances={}
        for name in eval_names:
            panel_obj=panels[name]
            if not np.array_equal(panel_obj.geometry.coordinates,expected[name]):
                raise RuntimeError(f"host self geometry drift for {name}")
            mask=read_canonical_mask_alignment(panel_obj.fasta_path)
            occ=read_occurrence_rows(panel_obj.occurrence_path)
            grouped=masks_by_frozen_locality(mask,occ,panel_obj.canonical_latlon)
            support=edge_mask_support(
                grouped,
                panel_obj.geometry.edge_nodes,
                alignment_length=mask.alignment_length,
                minimum_comparable_fraction=.5,
            )
            if not support.all_edges_valid:
                raise RuntimeError(f"host self mask admissibility drift for {name}")
            edge=extract_species_frozen_edge_distances(
                panel_obj.fasta_path,
                panel_obj.occurrence_path,
                panel_obj.canonical_latlon,
                panel_obj.geometry,
                minimum_comparable_fraction=.5,
            )
            geos[name]=panel_obj.geometry
            templates[name]=template_edges(name,panel_obj.geometry)
            ibd[name]=prepare_compact_crossfit_ibd_design(
                templates[name].length,
                templates[name].nodes,
                min_training_edges=5,
            )
            distances[name]=np.asarray(edge.genetic_distance,float)

        design=SimpleNamespace(
            template_edges=templates,
            ibd_designs=ibd,
            eval_species=eval_names,
            min_training_edges=5,
        )
        self_design=prepare_genetic_self_detectability(
            design,
            species=eval_names,
            bandwidth=float(rule["method_inheritance"]["bandwidth_km"]),
            prior_strength=float(rule["method_inheritance"]["prior_strength"]),
            prior_mean=float(rule["method_inheritance"]["prior_mean"]),
            segment_points=int(rule["method_inheritance"]["segment_points"]),
        )
        world=SimpleNamespace(genetic_distance=distances)
        scored=score_phylogatr_compact_self_world_batch(
            design,self_design,[world]
        )
        if int(scored.n_species[0])!=len(eval_names):
            raise RuntimeError("non-finite host empirical self species count")
        statistic=float(scored.statistics[0])
        species_scores={
            name:float(scored.species_scores[name][0])
            for name in eval_names
        }

    p=float(upper_monte_carlo_pvalue(statistic,reference))
    alpha=float(rule["method_inheritance"]["alpha"])
    detected=bool(p<=alpha)
    values=np.asarray(list(species_scores.values()),float)
    out={
        "schema":"ttf_lepidoptera_postprimary_empirical_self_v0.1",
        "status":"EMPIRICAL_SELF_COMPLETE",
        "panel":"host",
        "parent_cross_species_decision":"NO_DETECTED_POSITIVE_HOST_RESOURCE_GEOGRAPHY_GRADIENT",
        "design_npz_sha256":sha256_path(args.design_npz),
        "source_archive_sha256":SOURCE_ARCHIVE_SHA,
        "statistic":statistic,
        "p_value":p,
        "alpha":alpha,
        "detected_against_structural_null":detected,
        "evaluation_species":len(eval_names),
        "species_score_summary":{
            "minimum":float(np.min(values)),
            "q25":float(np.quantile(values,.25)),
            "median":float(np.median(values)),
            "q75":float(np.quantile(values,.75)),
            "maximum":float(np.max(values)),
        },
        "species_scores":dict(sorted(species_scores.items())),
        "interpretation":(
            "Within-species spatial genetic structure is detectable relative to the inherited structural null, while the completed host-resource cross-species gradient remains undetected."
            if detected else
            "Within-species self-detectability is not detected under the inherited structural-null calibration; the completed cross-species primary remains unchanged."
        ),
        "raw_statistic_sign_interpreted_against_zero":False,
        "cross_species_primary_rewritten":False,
        "serialized_sequence_identity":False,
        "serialized_edge_genetic_distance_vectors":False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({
        "statistic":statistic,
        "p_value":p,
        "detected":detected,
        "evaluation_species":len(eval_names),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
