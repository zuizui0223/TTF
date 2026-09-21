#!/usr/bin/env python3
from __future__ import annotations

import argparse,csv,hashlib,io,json,stat,tempfile,zipfile
from pathlib import Path
import numpy as np

import ttf.lepidoptera_host_resource_empirical as empirical_module
from ttf.lepidoptera_host_resource_empirical import score_empirical_host_resource
from ttf.lepidoptera_trait_gradient_empirical import post_ibd_edge_response
from ttf.phylogatr_character_mask import (
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
from ttf.phylogatr_empirical import extract_species_frozen_edge_distances

ALIASES=("COI","CO1","COX1","COXI","CYTOCHROME C OXIDASE SUBUNIT I","CYTOCHROME C OXIDASE SUBUNIT 1")
RULE_SCHEMA="ttf_lepidoptera_host_resource_empirical_rule_v0.1"
AUTH_SCHEMA="ttf_lepidoptera_host_resource_empirical_authorization_v0.1"


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
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
        raise RuntimeError("selective extraction drift")
    return root


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--design-npz",type=Path,required=True)
    ap.add_argument("--private-reference-npy",type=Path,required=True)
    ap.add_argument("--geometry-reference-npy",type=Path,required=True)
    ap.add_argument("--breadth-reference-npy",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--qualification",type=Path,required=True)
    ap.add_argument("--authorization",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!=RULE_SCHEMA or rule.get("status")!="FROZEN_BEFORE_HOST_RESOURCE_NUCLEOTIDE_IDENTITY_OPENING":
        raise RuntimeError("empirical rule is not frozen")
    if any(rule["outcome_firewall"].values()):
        raise RuntimeError("empirical rule firewall open")
    qualification=json.loads(args.qualification.read_text())
    if qualification.get("status")!="PASS" or qualification.get("empirical_identity_opening_eligible") is not True:
        raise RuntimeError("v0.2 qualification has not authorized identity opening")
    authorization=json.loads(args.authorization.read_text())
    if authorization.get("schema")!=AUTH_SCHEMA or authorization.get("status")!="AUTHORIZE_ONE_SHOT_HOST_RESOURCE_IDENTITY_OPENING":
        raise RuntimeError("empirical authorization missing")

    frozen=authorization["frozen_sha256"]
    checks={
        "empirical_rule":sha256_path(args.rule),
        "empirical_scorer":sha256_path(Path(empirical_module.__file__)),
        "empirical_runner":sha256_path(Path(__file__)),
        "qualification":sha256_path(args.qualification),
        "source_archive":sha256_path(args.archive),
        "survivor_design_npz":sha256_path(args.design_npz),
        "private_reference_npy":sha256_path(args.private_reference_npy),
        "geometry_reference_npy":sha256_path(args.geometry_reference_npy),
        "host_breadth_reference_npy":sha256_path(args.breadth_reference_npy),
    }
    if checks!=frozen:
        drift={k:(frozen.get(k),v) for k,v in checks.items() if frozen.get(k)!=v}
        raise RuntimeError(f"empirical authorization hash drift: {drift}")

    if checks["source_archive"]!=rule["source_archive"]["sha256"]:
        raise RuntimeError("source archive SHA drift")
    if checks["survivor_design_npz"]!=rule["survivor_design"]["design_npz_sha256"]:
        raise RuntimeError("survivor design SHA drift")
    refs=rule["reference_distributions"]
    if checks["private_reference_npy"]!=refs["private_npy_sha256"]:
        raise RuntimeError("private reference SHA drift")
    if checks["geometry_reference_npy"]!=refs["geometry_confounded_npy_sha256"]:
        raise RuntimeError("geometry reference SHA drift")
    if checks["host_breadth_reference_npy"]!=refs["host_breadth_confounded_npy_sha256"]:
        raise RuntimeError("breadth reference SHA drift")

    z=np.load(args.design_npz,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    if len(names)!=int(rule["survivor_design"]["species"]):
        raise RuntimeError("survivor species count drift")
    coords=np.asarray(z["coordinates"],float)
    offsets=np.asarray(z["coordinate_offsets"],np.int64)
    expected_coords={
        name:coords[offsets[i]:offsets[i+1]]
        for i,name in enumerate(names)
    }

    with tempfile.TemporaryDirectory(prefix="ttf_host_resource_empirical_") as td:
        root=selective_extract(args.archive,Path(td),set(names))
        genes=read_genes_rows(root/"genes.txt")
        universe={str(row.get("species","")).strip() for row in genes}
        scan=scan_phylogatr_phase1(
            root,aliases=ALIASES,excluded_species=universe-set(names),
            min_localities=12,min_endpoint_training_edges=5,neighbor_fraction=.15,
        )
        panels={p.species:p for p in choose_one_panel_per_species(scan.candidates)}
        if set(panels)!=set(names):
            raise RuntimeError("survivor source-panel reconstruction drift")

        # Response-blind integrity gate, including character masks.
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
                raise RuntimeError(f"survivor character-mask admissibility drift for {name}")

        # AUTHORIZED IDENTITY OPENING: all frozen integrity checks passed.
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
        breadth=np.load(args.breadth_reference_npy,allow_pickle=False)
        primary=score_empirical_host_resource(
            response,z,
            private_reference=private,
            geometry_reference=geometry,
            breadth_reference=breadth,
            alpha=.05,
        )

    values=np.asarray(list(primary.per_target.values()),float)
    out={
        "schema":"ttf_lepidoptera_host_resource_empirical_result_v0.1",
        "status":"PRIMARY_ONE_SHOT_DECISION_COMPLETE",
        "authorization_sha256":sha256_path(args.authorization),
        "source_archive_sha256":checks["source_archive"],
        "survivor_design_npz_sha256":checks["survivor_design_npz"],
        "primary_host_resource_gradient":{
            "statistic":primary.statistic,
            "envelope_p_value":primary.p_value,
            "alpha":.05,
            "positive":primary.positive,
            "finite_target_correlations":len(primary.per_target),
            "target_correlation_summary":{
                "minimum":float(np.min(values)),
                "q25":float(np.quantile(values,.25)),
                "median":float(np.median(values)),
                "q75":float(np.quantile(values,.75)),
                "maximum":float(np.max(values)),
                "positive_targets":int(np.count_nonzero(values>0)),
                "negative_targets":int(np.count_nonzero(values<0)),
            },
        },
        "decision":"POSITIVE_HOST_RESOURCE_GEOGRAPHY_GRADIENT" if primary.positive else "NO_DETECTED_POSITIVE_HOST_RESOURCE_GEOGRAPHY_GRADIENT",
        "interpretation":rule["one_shot_decision"][
            "interpretation_positive" if primary.positive else "interpretation_negative"
        ],
        "source_integrity":{
            "archive_sha256_verified":True,
            "survivor_design_sha256_verified":True,
            "survivor_mask_admissibility_reverified":True,
            "species":source_checks,
        },
        "outcome_state":{
            "sequence_identity_opened":True,
            "pairwise_genetic_distances_opened":True,
            "host_resource_transfer_statistic_computed":True,
            "serialized_sequence_identity":False,
            "serialized_edge_genetic_distance_vectors":False,
        },
        "post_result_retuning_allowed":False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({
        "decision":out["decision"],
        "primary_statistic":primary.statistic,
        "primary_p_value":primary.p_value,
        "finite_targets":len(primary.per_target),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
