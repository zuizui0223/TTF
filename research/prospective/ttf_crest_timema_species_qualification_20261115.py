"""TTF-Q source-only Timema dyad geometry and response-cell presence audit.

Source data are the SHA-locked publisher workbooks identified in the preceding
preflight. Decodes ONLY source identity columns and XML cell *presence* for
Figure3 habitat/sexual isolation outcomes; NEVER decodes response magnitudes.
"""
from __future__ import annotations
import collections
import io
import json
import statistics
import sys
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

# Imported directly by offline tests via spec_from_file_location, so sibling
# modules are not automatically on sys.path. Add only this checked source dir.
_SOURCE_DIR = str(Path(__file__).resolve().parent)
if _SOURCE_DIR not in sys.path:
    sys.path.insert(0, _SOURCE_DIR)

from ttf_crest_timema_publisher_species_holdout_geometry_20261115 import (
    M, PKG, REL, SOURCES, fetch, read_identity, split_geometry,
)
from posixpath import normpath

# Source-design admission thresholds: an evaluable host+sexual transfer contrast
# needs both within- and between-species dyads in each arm. Counts are only
# source geometry, not trial numbers, independent responses or statistical power.
MIN_WITHIN_PER_ARM=3
MIN_BETWEEN_PER_ARM=2
MIN_DISTINCT_SPECIES_PER_ARM=3


def connected_components(rows: list[dict]) -> list[dict]:
    nodes = sorted({x for r in rows for x in (r["species_1"],r["species_2"])})
    graph={n:set() for n in nodes}
    for r in rows:
        a,b=r["species_1"],r["species_2"]
        graph[a].add(b);graph[b].add(a)
    comps=[]
    seen=set()
    for n in nodes:
        if n in seen:continue
        group={n};pending=[n]
        while pending:
            for k in sorted(graph[pending.pop()]):
                if k not in group:group.add(k);pending.append(k)
        seen|=group
        inside=[r for r in rows if r["species_1"] in group and r["species_2"] in group]
        comps.append({"species":sorted(group),"species_count":len(group),
                      "dyads":len(inside),
                      "within_species_dyads":sum(r["species_1"]==r["species_2"] for r in inside),
                      "between_species_dyads":sum(r["species_1"]!=r["species_2"] for r in inside)})
    return sorted(comps,key=lambda x:(-x["dyads"],x["species"]))


def diversity_summary(rows: list[dict]) -> dict:
    counts=collections.Counter(x for r in rows for x in (r["species_1"],r["species_2"]))
    weight=sum(counts.values())
    simpson=weight**2/sum(k*k for k in counts.values()) if counts else None
    return {
        "dyads":len(rows),
        "species_involved":len(counts),
        "within_species":sum(r["species_1"]==r["species_2"] for r in rows),
        "between_species":sum(r["species_1"]!=r["species_2"] for r in rows),
        "effective_species_endpoint_diversity_not_effective_sample_size":simpson,
        "max_species_endpoint_incidence":max(counts.values(),default=0),
        "max_species_endpoint_incidence_fraction":max(counts.values(),default=0)/weight if weight else None,
    }


def split_support(rows: list[dict]) -> dict:
    geometry=split_geometry(rows)
    names=geometry["species_labels"]
    best=None;passing=0
    bestraw=geometry["best_balanced_species_split_source_only"]
    if bestraw is None:raise ValueError("no original publisher split")
    def info(train,test):
        tr=[r for r in rows if r["species_1"] in train and r["species_2"] in train]
        te=[r for r in rows if r["species_1"] in test and r["species_2"] in test]
        good=(all(len(x)>=MIN_DISTINCT_SPECIES_PER_ARM for x in (train,test))
              and all(sum(r["species_1"]==r["species_2"] for r in arm)>=MIN_WITHIN_PER_ARM for arm in (tr,te))
              and all(sum(r["species_1"]!=r["species_2"] for r in arm)>=MIN_BETWEEN_PER_ARM for arm in (tr,te)))
        return tr,te,good
    total=0
    for mask in range(1<<(len(names)-1)):
        train={names[0]}|{names[i+1] for i in range(len(names)-1) if mask&(1<<i)}
        test=set(names)-train
        if len(train)<MIN_DISTINCT_SPECIES_PER_ARM or len(test)<MIN_DISTINCT_SPECIES_PER_ARM:continue
        total+=1
        tr,te,ok=info(train,test)
        if not ok:continue
        passing+=1
        score=(min(len(tr),len(te)),min(sum(r["species_1"]!=r["species_2"] for r in arm) for arm in (tr,te)),
               len(tr)+len(te), -abs(len(tr)-len(te)))
        if best is None or score>best[0]:
            best=(score,{"train_species":sorted(train),"test_species":sorted(test),
                         "train":diversity_summary(tr),"test":diversity_summary(te),
                         "cross_boundary_dyads_excluded":len(rows)-len(tr)-len(te)})
    base_tr,base_te,_=info(set(bestraw["train_species"]),set(bestraw["test_species"]))
    return {
        "source_only_preselected_balanced_partition":{
            "train_species":bestraw["train_species"],"test_species":bestraw["test_species"],
            "train":diversity_summary(base_tr),"test":diversity_summary(base_te),
            "cross_boundary_dyads_excluded":len(rows)-len(base_tr)-len(base_te),
        },
        "minimum_within_per_arm":MIN_WITHIN_PER_ARM,
        "minimum_between_per_arm":MIN_BETWEEN_PER_ARM,
        "minimum_species_per_arm":MIN_DISTINCT_SPECIES_PER_ARM,
        "possible_species_partitions_examined":total,
        "source_geometry_qualified_partitions":passing,
        "best_both_types_supported_partition":best[1] if best else None,
        "status":"GEOMETRY_ONLY" if passing else "NOT_EVALUABLE_BOTH_ISOLATION_TYPES_UNDER_SOURCE_THRESHOLDS",
    }


def source_response_cell_presence(raw:bytes) -> dict:
    """Check only cells K/L are structurally populated on Figure3 rows."""
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap={r.get("Id"):r.get("Target") for r in rels.findall(f"{{{PKG}}}Relationship")}
        sh=next((x for x in wb.findall(f".//{{{M}}}sheets/{{{M}}}sheet")
                 if x.get("name")=="Figure3"),None)
        if sh is None:raise ValueError("Figure3 source sheet missing")
        path=relmap.get(sh.get(f"{{{REL}}}id"))
        if not path:raise ValueError("Figure3 source mapping missing")
        path=(path.lstrip("/") if path.startswith("/xl/") else normpath("xl/"+path))
        sheet=ET.fromstring(z.read(path))
        body=sheet.find(f"{{{M}}}sheetData")
        if body is None:raise ValueError("Figure3 source rows missing")
        rows=body.findall(f"{{{M}}}row")
        counts={"habitat_isolation":0,"sexual_isolation":0,
                "both_cells_structurally_present":0}
        for row in rows:
            idx=int(row.get("r","0"))
            if idx<2:continue
            indexed={c.get("r"):c for c in row.findall(f"{{{M}}}c")}
            present=[]
            for col in ("K","L"):
                cell=indexed.get(f"{col}{idx}")
                exists=cell is not None and (
                    cell.find(f"{{{M}}}v") is not None
                    or cell.find(f"{{{M}}}is") is not None
                    or cell.find(f"{{{M}}}f") is not None)
                present.append(exists)
            counts["habitat_isolation"]+=int(present[0])
            counts["sexual_isolation"]+=int(present[1])
            counts["both_cells_structurally_present"]+=int(all(present))
        counts["row_elements"]=len(rows)
        counts["meaning"]="serialized XML cells only; numeric values not read; not valid response trial count"
        return counts


def main():
    x=fetch("source_data");y=fetch("supplement_1")
    a=read_identity(x,"source_data");b=read_identity(y,"supplement_1")
    aa={r["pair_id"]:r for r in a};bb={r["pair_id"]:r for r in b}
    if set(aa)!=set(bb) or aa!=bb:
        raise ValueError("publisher independently sourced pair/species/population keys disagree")
    print(json.dumps({
        "status":"SOURCE_ONLY_TTF_Q_TIMEMA_SPECIES_SUPPORT_AND_RESPONSE_PRESENCE",
        "source_sha256":{k:p["sha256"] for k,p in SOURCES.items()},
        "source_observation":diversity_summary(a),
        "species_graph_components":connected_components(a),
        "holdout_support":split_support(a),
        "figure3_outcome_cell_presence_without_values":source_response_cell_presence(x),
        "response_values_decoded":False,
        "no_prediction_or_effect_estimated":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
