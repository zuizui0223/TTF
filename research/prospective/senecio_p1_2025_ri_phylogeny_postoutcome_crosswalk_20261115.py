"""Source-pinned post-outcome descriptive crosswalk of Senecio 2021 P1 and 2025 RI.

This is NOT a preregistered test, a genotype re-opening, or a new TTF
heldout prediction. It only joins frozen P1 aggregate pair IDs to original
author 2025 phylogeny tip labels and published averaged crossing results.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import re
import urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
BASE="https://raw.githubusercontent.com/MaddieEJames/James_et_al._2025_MolecularEcology"
COMMIT="fd8e8f075a72c1420267eecfd482dbc6a312caef"
PINNED={
  "data/Phylogeny.nex": "a7c79d4b0cccb4b628e46fdcce45f417dc487d38",
  "data/CrossesDH-ReciprocalsAveraged.txt": "ae6147df032baf492d692dea33aa03398b794d78"
}
P1_FILE=Path(__file__).with_name("senecio_p1_first_opening_aggregate_receipt_20261115.json")
P0_P1_FILE=Path(__file__).with_name("senecio_p0_p1_crosslayer_descriptive_20261115.json")
# Figure 2D group labels and Figure 2C phylogenetic brackets,
# transcribed as the *authors' source-defined clade labels*, not inferred
# from the P1 genotype source. Do not add geography for unmatched pairs.
PAPER_FIG2_CLADES={
  "D01-H01":"Eastern", "D03-H02":"Eastern", "D04-H05":"Eastern",
  "D14-H15":"Southeastern", "D32-H12":"Southeastern", "D23-H21":"Southeastern",
}
PAPER="https://doi.org/10.1111/mec.70090"


def git_blob_sha(raw: bytes) -> str:
    return hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()


def fetch_pinned(path: str) -> bytes:
    raw=urllib.request.urlopen(f"{BASE}/{COMMIT}/{path}",timeout=35).read(100_000)
    if git_blob_sha(raw)!=PINNED[path]:
        raise ValueError(f"original 2025 Git blob mismatch: {path}")
    return raw


def read_tree_tips(text: str) -> set[str]:
    # Only the declared NEXUS taxlabels block, not regex matches
    # of tip substrings in phylogeny annotation/comment text.
    m=re.search(r"begin\s+taxa\s*;.*?dimensions\s+ntax\s*=\s*(\d+)\s*;\s*taxlabels\s+(.*?)\s*;",text,re.I|re.S)
    if not m:
        raise ValueError("NEXUS taxa metadata block missing")
    tips=m.group(2).split()
    if len(tips)!=int(m.group(1)) or len(tips)!=len(set(tips)):
        raise ValueError("duplicate or incomplete phylogeny taxon labels")
    return set(tips)


def read_ri_pair_table(raw: bytes) -> dict:
    rd=csv.reader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t")
    hdr=next(rd)
    expected=["CrossType","RIseedSet-SobelChen","RIviability-SobelChen",
       "CrossType","Ecotype","Clade","Divergence time",
       "EnvDist-MDS-5PCs","PhenoDist-MDS-7PCs"]
    if hdr!=expected:
        raise ValueError("2025 RI source schema differs or duplicate-header positions changed")
    rows={}
    n=0
    for cells in rd:
        if not cells or not any(cells):
            continue
        n+=1
        if len(cells)!=len(hdr):
            raise ValueError("incomplete crossing-table row")
        pair=cells[0]
        if pair in rows:
            raise ValueError(f"duplicate source RI comparison {pair}")
        if cells[3] not in ("DH","DD","HH","HD") or cells[4] not in ("Same","Different"):
            raise ValueError(f"unrecognized source ecotype semantics: {pair}")
        rows[pair]={
           "ecotype_cross":cells[3], "ecotype_contrast":cells[4],
           "clade_comparison":cells[5],
           "ri_seed_set":float(cells[1]), "ri_viability":float(cells[2]),
        }
    if n!=18:
        raise ValueError(f"2025 source-crossing table row count drift: {n}")
    return rows


def join_sources(p1: dict, cross: dict, tips: set[str], p0p1: dict) -> dict:
    if p1.get("source_execution",{}).get("run_id")!=43693277926:
        raise ValueError("P1 original one-shot source not pinned")
    folds=p1["folds"]
    if len(folds)!=8 or len({f["pair"] for f in folds})!=8:
        raise ValueError("P1 frozen pair unit drift")
    pheno={r["pair"]:r for r in p0p1["rows"]}
    if set(pheno)!={f["pair"] for f in folds}:
        raise ValueError("P0/P1 source pair keys differ")
    items=[]
    for fold in folds:
        pair=fold["pair"]
        a,b=pair.split("-")
        tree_match=a in tips and b in tips
        ri=cross.get(pair)
        if ri is not None and not tree_match:
            raise ValueError(f"crossing row has no matching independent neutral tree tips: {pair}")
        if ri is not None and (ri["ecotype_cross"]!="DH" or ri["ecotype_contrast"]!="Different"):
            raise ValueError(f"wrong crossing role for P1 paired ecotypes: {pair}")
        items.append({
            "pair":pair,
            "p0_trait_directions_agree_of_9":pheno[pair]["morphology_signs_agree"],
            "p1_source_frozen_auc":fold["auc"],
            "p1_selected_snp_direction_fraction":fold["signed_snp_direction_fraction"],
            "2025_neutral_tree_both_tips":tree_match,
            "2025_paper_fig2_clade":PAPER_FIG2_CLADES.get(pair) if tree_match else None,
            "2025_direct_intrinsic_cross":ri,
            "not_same_individuals":True,
        })
    return {
      "status":"POSTOUTCOME_DESCRIPTIVE_2021_2025_LINKAGE_NOT_TTF_CONFIRMATORY",
      "original_p1_run_id":43693277926,
      "2025_author_commit":COMMIT,
      "2025_author_blob_shas":PINNED,
      "2025_paper":PAPER,
      "2025_ri_table_rows":len(cross),
      "2025_nexus_declared_taxa":len(tips),
      "p1_pair_count":len(items),
      "pairs_with_2025_tree_both_tips":sum(r["2025_neutral_tree_both_tips"] for r in items),
      "pairs_with_direct_2025_ri_cross":sum(r["2025_direct_intrinsic_cross"] is not None for r in items),
      "rows":items,
      "claims_not_licensed":[
         "causal_SNP_effect","cross_species_TTF","ancestral_mechanism_identification",
         "genomic_AUC_predicts_reproductive_isolation","2025_independent_genotype_replication",
         "sample_id_level_phenotype_genotype_cross"
      ],
      "no_new_original_vcf_or_genotype_rows_read":True,
    }


def main() -> None:
    p1=json.loads(P1_FILE.read_text())
    p0p1=json.loads(P0_P1_FILE.read_text())
    sources={p:fetch_pinned(p) for p in PINNED}
    tips=read_tree_tips(sources["data/Phylogeny.nex"].decode("utf-8"))
    cross=read_ri_pair_table(sources["data/CrossesDH-ReciprocalsAveraged.txt"])
    print(json.dumps(join_sources(p1,cross,tips,p0p1),indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":
    main()
