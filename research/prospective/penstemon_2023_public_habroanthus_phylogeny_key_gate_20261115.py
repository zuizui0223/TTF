"""Penstemon independent 2023 public Habroanthus tree/key-source coverage gate.

This public historical repo is not the deposited canonical 2026 manuscript tree.
Audits 2023 nuclear CDS and 10kb window species trees against author-written
hummingbird/bee group labels. Only tree tip labels + metadata keys are read;
no SNPs, gene sequences, floral measurements or genetic effect estimates.
"""
from __future__ import annotations
import collections
import csv
import hashlib
import io
import json
import re
import urllib.request
from pathlib import Path
import dendropy

REPO="lcwheeler/penstemon_novaseq"
COMMIT="fd93699387a8"
FILES={
 "CDS_nuclear_tree":{
   "path":"4_analysis/concordance_conflict_analyses/CDS_analyses/CDS_WithOutgroups_concord.cf.tree",
   "git_blob_sha":"ca1d2072e7df72649a8ec8699e31d8397fc46cae"},
 "10kb_window_tree":{
   "path":"4_analysis/concordance_conflict_analyses/genomic_windows_analyses/genomic_windows_WithOutgroups_concord.cf.tree",
   "git_blob_sha":"21f35d1490714d4584a0cd6f1af419243d3f3166"},
 "target_trios_groupmap":{
   "path":"4_analysis/twisst_analysis/target_trios/groups.tsv",
   "git_blob_sha":"360d904346db4ebc17b13d559aed5cd13dde2867"},
}
SOURCE_AUTHOR_REPO_LAST_COMMIT_DATE="2023-06-02"
PUB_2026_DOI="10.64898/2026.07.21.739878"
EXPECTED_ORIGINS={"barbatus","cardinalis","eatonii","labrosus"}
DECLARED_2026_HUMMINGBIRD_BEE_SISTER_CONTRAST={
    "barbatus":"virgatus",
    "cardinalis":"glaber",
    "eatonii":"laevis",
    "labrosus":"putus",
}


def git_blob_sha(raw:bytes)->str:
    return hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()


def fetch_source(info:dict)->bytes:
    url=f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{info['path']}"
    with urllib.request.urlopen(url,timeout=40) as response:
        raw=response.read(70_000)
    if git_blob_sha(raw)!=info["git_blob_sha"]:
        raise ValueError("pinned original 2023 Git blob mismatch")
    return raw


def strip_2023_internal_node_figure_annotations(newick: str) -> str:
    """Remove ONLY original IQ-TREE quote-annotation/slash-support after ')'.

    The 2023 published tree is nonstandard for DendroPy at
    )'[q1=...;q2=...;...]'/31.6:0.711; tip labels/branch lengths stay.
    This normalization authorizes TIP ID checks ONLY, never bootstrap inference.
    """
    pattern=r"(?<=\))'\[[^']+\]'/[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?=:)"
    return re.sub(pattern, "", newick)


def tip_labels(raw:bytes)->list[str]:
    newick=strip_2023_internal_node_figure_annotations(raw.decode("utf-8-sig"))
    phylo=dendropy.Tree.get(data=newick,schema="newick",preserve_underscores=True)
    labels=sorted(x.taxon.label.strip() for x in phylo.leaf_node_iter())
    if not labels or len(set(labels))!=len(labels):
        raise ValueError("repeated or empty 2023 original tree tip")
    return labels


def group_map(raw:bytes)->dict[str,str]:
    labels={}
    for row in csv.reader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t"):
        if not row:continue
        if len(row)!=2 or not all(row):
            raise ValueError("unexpected author 2023 target trios metadata columns")
        tip,group=row
        if tip in labels:raise ValueError("duplicate source tip mapping")
        labels[tip]=group
    if not labels:raise ValueError("empty original group mapping")
    return labels


def main():
    src={key:fetch_source(entry) for key,entry in FILES.items()}
    a=set(tip_labels(src["CDS_nuclear_tree"]))
    b=set(tip_labels(src["10kb_window_tree"]))
    groups=group_map(src["target_trios_groupmap"])
    mapped=set(groups)
    found=set(groups.values())
    namecount=collections.Counter(groups.values())
    sister_pairs={
        hb:{"hummingbird_group_sample_count":namecount.get(hb,0),
            "sister_bee_group":bee,
            "sister_bee_group_sample_count":namecount.get(bee,0),
            "both_source_groups_sampled_and_mapped":bool(namecount.get(hb,0) and namecount.get(bee,0))}
        for hb,bee in DECLARED_2026_HUMMINGBIRD_BEE_SISTER_CONTRAST.items()
    }
    print(json.dumps({
      "status":"PUBLIC_2023_HABROANTHUS_TREE_KEY_COVERAGE_NOT_CANONICAL_2026",
      "repository":REPO,"pinned_source_commit_prefix":COMMIT,
      "last_source_repo_commit_date":SOURCE_AUTHOR_REPO_LAST_COMMIT_DATE,
      "2026_manuscript":PUB_2026_DOI,
      "CDS_tree_tip_count":len(a),"10kb_window_tree_tip_count":len(b),
      "CDS_vs_window_tip_ids_equal":a==b,
      "tree_tips_shared":len(a&b),
      "CDS_only_tip_ids":sorted(a-b),"window_only_tip_ids":sorted(b-a),
      "2026_origin_sister_contrasts_on_older_2023_mapped_source":sister_pairs,
      "2026_origin_sister_contrasts_both_groups_mapped":sum(x["both_source_groups_sampled_and_mapped"] for x in sister_pairs.values()),
      "original_tree_may_contain_unmapped_sister_taxa_do_not_assume_taxon_absence":True,
      "group_map_rows":len(groups),
      "group_map_distinct_species_label_count":len(found),
      "group_map_species_labels":dict(sorted(collections.Counter(groups.values()).items())),
      "group_map_tips_mapped_to_both_tree_sets":len(mapped&a&b),
      "group_map_tips_not_in_CDS_tree":sorted(mapped-a),
      "group_map_tips_not_in_window_tree":sorted(mapped-b),
      "2026_four_hummingbird_origins_named_in_older_groupmap":sorted(EXPECTED_ORIGINS&found),
      "2023_source_sha":{k:v["git_blob_sha"] for k,v in FILES.items()},
      "same_2026_sample_specimens_not_verified":True,
      "no_floral_or_genomic_values_opened":True,
      "not_a_species_holdout_TTF_result":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
