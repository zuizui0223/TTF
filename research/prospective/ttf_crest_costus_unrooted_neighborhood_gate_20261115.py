"""Costus source-only unrooted phylogenetic-neighborhood heldout TTF-Q gate.

Algorithm fixed in COSTUS_UNROOTED_PHYLO_BLOCK_SOURCE_GATE_20261115.md
before any scent-response values were opened. Reads only original tree's
patristic distances, taxon keys and source categorical pollinator labels.
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
    pinned,TREE_NAME,INFO_NAME,MEANS_NAME,pollination_labels,means_row_labels
)
from ttf_crest_costus_clade_holdout_source_gate_20261115 import source_tree,clade_folds

K_ORDER=[4,3,2]
DISTANCE_ROUND_DIGITS=12
SOURCE_SPP_MIN=27


def exact_patristic_distances(tree, labels: list[str]) -> dict[tuple[str,str],float]:
    lookup={node.taxon.label.strip().lower():node.taxon for node in tree.leaf_node_iter()}
    if not set(labels).issubset(lookup):
        raise ValueError("unrooted source tree missing mean-scent tip labels")
    dm=tree.phylogenetic_distance_matrix()
    distances={}
    for i,a in enumerate(labels):
        for b in labels[i+1:]:
            val=float(dm.distance(lookup[a],lookup[b]))
            if val<0:
                raise ValueError("negative patristic distance unsupported")
            distances[(a,b)]=val
    return distances


def avg_linkage_partitions(labels:list[str],dist:dict[tuple[str,str],float]) -> dict[int,list[list[str]]]:
    clusters=[tuple([x]) for x in sorted(labels)]
    snapshots={}
    while len(clusters)>=2:
        k=len(clusters)
        if k in K_ORDER:
            snapshots[k]=[list(x) for x in sorted(clusters,key=lambda a:(-len(a),a))]
        ranked=[]
        for i,a in enumerate(clusters):
            for j in range(i+1,len(clusters)):
                b=clusters[j]
                cross=[dist[tuple(sorted((x,y)))] for x in a for y in b]
                average=sum(cross)/len(cross)
                ranked.append((round(average,DISTANCE_ROUND_DIGITS),a,b,i,j))
        if not ranked:break
        _,_,_,i,j=min(ranked)
        merged=tuple(sorted(clusters[i]+clusters[j]))
        clusters=[cluster for z,cluster in enumerate(clusters) if z not in (i,j)]
        clusters.append(merged)
        clusters.sort()
    return snapshots


def main():
    poll=pollination_labels(pinned(INFO_NAME))
    mean_labels=means_row_labels(pinned(MEANS_NAME))
    if set(poll)!=set(mean_labels) or len(mean_labels)!=SOURCE_SPP_MIN:
        raise ValueError("2025 source taxon-key mismatch")
    tree=source_tree(pinned(TREE_NAME))
    source_rooted=bool(tree.is_rooted)
    d=exact_patristic_distances(tree,mean_labels)
    groups=avg_linkage_partitions(mean_labels,d)
    candidates={}
    select=None
    for k in K_ORDER:
        p=clade_folds(groups.get(k),poll)
        candidates[str(k)]=p
        if select is None and p["status"]=="GEOMETRY_QUALIFIED":
            select=k
    print(json.dumps({
        "status":"COSTUS_UNROOTED_PHYLO_NEIGHBORHOOD_SOURCE_GATE",
        "source_original_tree_rooted":source_rooted,
        "source_original_tree_sha256":"24f1cbfaa82ce4d33aad1ba80787cc611a6730aafb38f8244954c4d0d50ae7a1",
        "source_scent_taxa":len(mean_labels),
        "source_tree_pairwise_distances":len(d),
        "distance_method":"DendroPy patristic sum, unrooted valid",
        "clustering":"average-linkage with rounded-12-decimal distance lexicographic ties",
        "candidate_K_order":K_ORDER,
        "source_support_checks":candidates,
        "first_qualified_K":select,
        "result":"GEOMETRY_ONLY" if select else "SOURCE_GEOMETRY_NOT_EVALUABLE",
        "not_true_independent_pollinator_transitions":True,
        "original_pollinator_labels_not_authenticated_observed_visitors":True,
        "scent_response_values_read":False,
    },sort_keys=True,indent=2,ensure_ascii=False))


if __name__=="__main__":main()
