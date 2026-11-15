"""Response-blind Costus fixed phylogeny-cluster holdout qualification.

Tree split algorithm and finite candidate order K=4,3,2 were fixed before
opening original 27 species' scent response columns. Pollination groups are
source labels, NOT authenticated field visitation labels. Original tree is
SHA-locked and only the original rooted topology is used.
"""
from __future__ import annotations
import collections
import json
import sys
from pathlib import Path
import dendropy

DIR=str(Path(__file__).resolve().parent)
if DIR not in sys.path:sys.path.insert(0,DIR)
from ttf_crest_costus_original_tree_label_geometry_20261115 import (
    pinned, TREE_NAME, INFO_NAME, MEANS_NAME, pollination_labels, means_row_labels
)

CANDIDATE_CLUSTERS=[4,3,2]
MIN_TEST_TIPS_PER_CLADE=3
MIN_TEST_TIPS_PER_SOURCE_POLLINATOR_CLASS=1
MIN_TRAIN_PER_SOURCE_POLLINATOR_CLASS=4

def name(node):
    return node.taxon.label.strip().lower()

def group_tip_labels(node):
    return sorted(name(tip) for tip in node.leaf_iter())

def source_tree(raw):
    s=raw.decode("utf-8-sig")
    schema="nexus" if s.lstrip().lower().startswith("#nexus") else "newick"
    return dendropy.Tree.get(data=s,schema=schema,preserve_underscores=True)

def split_largest_subclades(tree,k:int)->list[list[str]]|None:
    """Deterministically split largest currently retained subtree (no outcomes)."""
    groups=[tree.seed_node]
    while len(groups)<k:
        splittable=[n for n in groups if len(n.child_nodes())>=2]
        if not splittable:return None
        splittable.sort(key=lambda x:(-len(group_tip_labels(x)),tuple(group_tip_labels(x))))
        node=splittable[0]
        children=list(node.child_nodes())
        if len(groups)-1+len(children)>k:return None
        groups.remove(node)
        groups.extend(children)
    return sorted([group_tip_labels(x) for x in groups],key=lambda x:(-len(x),x))

def clade_folds(groups:list[list[str]], poll:dict[str,str])->dict:
    if not groups:return{"status":"NOT_EVALUABLE_TOPOLOGY_UNSPLITTABLE"}
    allnames=set(poll)
    if set().union(*map(set,groups))!=allnames or sum(map(len,groups))!=len(allnames):
        raise ValueError("source clades do not partition exact 27 scent taxa")
    fold=[]
    allgood=True
    for i,g in enumerate(groups):
        held=set(g)
        training=allnames-held
        test_counts=dict(sorted(collections.Counter(poll[x] for x in held).items()))
        train_counts=dict(sorted(collections.Counter(poll[x] for x in training).items()))
        qualifies=(len(held)>=MIN_TEST_TIPS_PER_CLADE and
            test_counts.get("Bee",0)>=MIN_TEST_TIPS_PER_SOURCE_POLLINATOR_CLASS and
            test_counts.get("Hummingbird",0)>=MIN_TEST_TIPS_PER_SOURCE_POLLINATOR_CLASS and
            train_counts.get("Bee",0)>=MIN_TRAIN_PER_SOURCE_POLLINATOR_CLASS and
            train_counts.get("Hummingbird",0)>=MIN_TRAIN_PER_SOURCE_POLLINATOR_CLASS)
        allgood &= qualifies
        fold.append({"fold":i+1,"test_taxa":g,"test_n":len(g),
                     "test_pollination_source_counts":test_counts,
                     "training_pollination_source_counts":train_counts,
                     "passes_source_only_support":qualifies})
    return {"status":"GEOMETRY_QUALIFIED" if allgood else "INSUFFICIENT_POLLINATOR_OR_CLADE_SUPPORT",
            "folds":fold,"n_folds":len(fold)}

def main():
    poll=pollination_labels(pinned(INFO_NAME))
    labels=means_row_labels(pinned(MEANS_NAME))
    if set(poll)!=set(labels):raise ValueError("source 2025 scent mean and pollination label taxon mismatch")
    tree=source_tree(pinned(TREE_NAME))
    original_tree_tip_count=sum(1 for _ in tree.leaf_node_iter())
    source_rooted=bool(tree.is_rooted)
    if not source_rooted:
        print(json.dumps({"status":"NOT_EVALUABLE_ORIGINAL_CHRONOGRAM_NOT_ROOTED",
              "original_tree_tips":original_tree_tip_count,"scent_values_read":False},indent=2))
        return
    original_case_labels=[tip.taxon.label for tip in tree.leaf_node_iter()
                          if name(tip) in poll]
    tree.retain_taxa_with_labels(original_case_labels)
    post=set(name(t) for t in tree.leaf_node_iter())
    if post!=set(poll):raise ValueError("source pruned tree tip identities mismatch")
    evaluated={}
    selected=None
    for k in CANDIDATE_CLUSTERS:
        groups=split_largest_subclades(tree,k)
        outcome=clade_folds(groups,poll)
        evaluated[str(k)]=outcome
        if selected is None and outcome["status"]=="GEOMETRY_QUALIFIED":
            selected=k
    print(json.dumps({
         "status":"COSTUS_PHYLOGENETIC_BLOCK_SOURCE_ONLY_QUALIFICATION",
         "source_2022_chronogram_original_tip_count":original_tree_tip_count,
         "source_rooted":source_rooted,
         "matched_scent_tree_taxa":len(post),
         "predeclared_candidate_k_order":CANDIDATE_CLUSTERS,
         "predeclared_test_min_class":MIN_TEST_TIPS_PER_SOURCE_POLLINATOR_CLASS,
         "predeclared_test_min_tips":MIN_TEST_TIPS_PER_CLADE,
         "predeclared_train_min_class":MIN_TRAIN_PER_SOURCE_POLLINATOR_CLASS,
         "geometry_checks":evaluated,
         "first_source_qualified_k":selected,
         "endpoint":"scent compounds (unopened)",
         "pollination_status":"author source categories, not verified directly observed",
         "no_scent_effect_magnitudes_read":True,
         "no_post_outcome_partition_selection":True,
       },indent=2,sort_keys=True,ensure_ascii=False))

if __name__=="__main__":main()
