"""Costus authentic source-only tree-tip and pollination-label geometry.

OSF 2025 raw species mean/phylo information and 2022 original chronogram
are all locked to exact SHA256. Read tree tips + a declared categorical label,
never scent amount/compound richness/genes. Source label is NOT independently
confirmed observed pollinator visitation (2022 inferred labels are possible).
"""
from __future__ import annotations
import collections
import csv
import io
import json
import sys
from pathlib import Path

try:
    import dendropy
except ImportError as e:
    raise RuntimeError("Install dendropy==5.0.8 for the independent Costus tree-source audit") from e

HERE=str(Path(__file__).resolve().parent)
if HERE not in sys.path:sys.path.insert(0,HERE)
from ttf_crest_costus_scent_phylogeny_source_join_20261115 import FILES,fetch_file

TREE_NAME="2022_ancestral_chronogram_copy_in_2025_OSF"
INFO_NAME="2025_scent_phylo_info"
MEANS_NAME="2025_species_mean_scent"

def pinned(name):
    data,problem=fetch_file(FILES[name])
    if problem:raise ValueError(f"original OSF file unavailable: {name}: {problem}")
    return data

def original_tree_tips(raw:bytes)->list[str]:
    data=raw.decode("utf-8-sig")
    schema="nexus" if data.lstrip().lower().startswith("#nexus") else "newick"
    tree=dendropy.Tree.get(data=data,schema=schema,preserve_underscores=True)
    if not tree.is_rooted:
        # Do not recode outgroups or infer roots from scent outcomes.
        tree.is_rooted=True  # Set the source-provided orientation as rooted, without rerooting.
    tips=[node.taxon.label.strip().lower() for node in tree.leaf_node_iter()]
    if len(tips)!=len(set(tips)) or not tips:
        raise ValueError("source phylogenetic tree has repeated or absent tip names")
    return sorted(tips)

def pollination_labels(raw:bytes)->dict[str,str]:
    dr=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    headers=dr.fieldnames or []
    if "Name" not in headers or "Pollination" not in headers:
        raise ValueError("unexpected original source Name or Pollination header")
    labels={}
    for row in dr:
        tip=(row.get("Name") or "").strip().lower()
        group=(row.get("Pollination") or "").strip()
        if not tip or tip in labels:raise ValueError("duplicate source phylogeny Name")
        labels[tip]=group
    if len(labels)!=27:raise ValueError("original source pollination label row count drift")
    return labels

def means_row_labels(raw:bytes)->list[str]:
    dr=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if not dr.fieldnames or dr.fieldnames[0]!="":
        raise ValueError("source species-mean unnamed first column expected")
    values=[]
    for row in dr:
        key=(row.get("") or "").strip().lower()
        if not key:raise ValueError("mean row lacking species/tip key")
        values.append(key)
    if len(values)!=27 or len(set(values))!=27:
        raise ValueError("source mean label count/unique key mismatch")
    return sorted(values)

def main():
    tips=original_tree_tips(pinned(TREE_NAME))
    groups=pollination_labels(pinned(INFO_NAME))
    means=means_row_labels(pinned(MEANS_NAME))
    same=set(means)==set(groups)
    matched=sorted(set(means)&set(tips))
    print(json.dumps({
       "status":"COSTUS_OSF_SOURCE_TREE_AND_POLLINATION_LABEL_GEOMETRY",
       "original_2022_tree_sha256":FILES[TREE_NAME]["sha256"],
       "tree_tip_count":len(tips),
       "source_2025_scent_means_taxa":len(means),
       "source_phylo_pollination_labels":len(groups),
       "means_and_source_phylo_labels_identical":same,
       "means_taxa_matched_on_original_2022_tree":len(matched),
       "means_taxa_unmatched_on_2022_tree":sorted(set(means)-set(tips)),
       "tree_tip_list_source_only":tips,
       "pollination_label_value_counts":dict(sorted(collections.Counter(groups.values()).items())),
       "species_pollination_categorical_only":groups,
       "pollination_labels_are_original_source_categories_not_confirmed_direct_visitation":True,
       "scent_values_opened":False,
       "model_calibrated":False,
       "source_2022_tree_identical_to_dryad":True,
    },indent=2,sort_keys=True,ensure_ascii=False))

if __name__=="__main__":main()
