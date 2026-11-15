"""Strict original-XLSX metadata-only Timema mating availability by dyad.

No sexual isolation scores, Fst, host preference or mate trial counts decoded.
Only XLSX cell storage type and exact source pair identity are inspected.
"""
from __future__ import annotations
import collections
import hashlib
import io
import json
import sys
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile
from posixpath import normpath

DIR=str(Path(__file__).resolve().parent)
if DIR not in sys.path:sys.path.insert(0,DIR)
from ttf_crest_timema_publisher_species_holdout_geometry_20261115 import (
    M, PKG, REL, SOURCES, fetch, read_identity, split_geometry
)

def source_mating_cell_storage(original_xlsx:bytes, identity_rows:list[dict])->dict:
    if hashlib.sha256(original_xlsx).hexdigest()!=SOURCES["source_data"]["sha256"]:
        raise ValueError("original XLSX SHA256 mismatch")
    with zipfile.ZipFile(io.BytesIO(original_xlsx)) as z:
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap={r.get("Id"):r.get("Target") for r in rels.findall(f"{{{PKG}}}Relationship")}
        sh=next((x for x in wb.findall(f".//{{{M}}}sheets/{{{M}}}sheet") if x.get("name")=="Figure3"),None)
        if sh is None:raise ValueError("no Figure3")
        path=relmap.get(sh.get(f"{{{REL}}}id"))
        if not path:raise ValueError("broken Figure3 relationship")
        p=path.lstrip("/") if path.startswith("/xl/") else normpath("xl/"+path)
        rows=ET.fromstring(z.read(p)).find(f"{{{M}}}sheetData").findall(f"{{{M}}}row")
    actual=[r for r in rows if int(r.get("r","0"))>=2]
    if len(actual)!=42 or len(identity_rows)!=42:
        raise ValueError("identity source row count drift")
    if [int(r.get("r","0")) for r in actual]!=list(range(2,44)):
        raise ValueError("noncontiguous source identities prevent safe row join")
    kinds={}
    out=[]
    for identity,row in zip(identity_rows,actual):
        indexed={c.get("r"):c for c in row.findall(f"{{{M}}}c")}
        idx=int(row.get("r"))
        cell=indexed.get(f"L{idx}")
        if cell is None:
            klass="no_serialized_cell"
        elif cell.find(f"{{{M}}}f") is not None:
            klass="formula_not_evaluated"
        elif cell.get("t","n")=="n" and cell.find(f"{{{M}}}v") is not None:
            klass="numeric_literal_storage_not_authenticated_trials"
        elif cell.get("t")=="s":
            klass="shared_string_storage_not_decoded"
        elif cell.get("t") in ("inlineStr","str"):
            klass="text_storage_not_decoded"
        else:
            klass="unknown_storage_not_decoded"
        kinds[klass]=kinds.get(klass,0)+1
        out.append({**identity,"sexual_cell_storage_class":klass})
    return {"rows":out,"classes":dict(sorted(kinds.items()))}


def count_by_arm(rows,train,test):
    arms={}
    for name,sp in (("train",set(train)),("evaluation",set(test))):
        sample=[r for r in rows if r["species_1"] in sp and r["species_2"] in sp]
        numeric=[r for r in sample if r["sexual_cell_storage_class"]=="numeric_literal_storage_not_authenticated_trials"]
        arms[name]={
            "source_dyads":len(sample),
            "source_within":sum(r["species_1"]==r["species_2"] for r in sample),
            "source_between":sum(r["species_1"]!=r["species_2"] for r in sample),
            "numeric_literal_dyads":len(numeric),
            "numeric_literal_within":sum(r["species_1"]==r["species_2"] for r in numeric),
            "numeric_literal_between":sum(r["species_1"]!=r["species_2"] for r in numeric),
            "nonnumeric_or_unknown_storage_dyads":len(sample)-len(numeric),
            "species_labels":sorted(sp),
        }
    arms["cross_boundary_dyads"]=len(rows)-arms["train"]["source_dyads"]-arms["evaluation"]["source_dyads"]
    return arms


def evaluate_support(rows):
    raw=split_geometry(rows)
    source_split=raw["best_balanced_species_split_source_only"]
    assert source_split is not None
    checked=count_by_arm(rows,source_split["train_species"],source_split["test_species"])
    species=raw["species_labels"]
    strong=0
    best=None
    # EXACT source-only declared minimums, now restricted to stored numeric mating fields.
    for bits in range(1<<(len(species)-1)):
        train={species[0]}|{species[i+1] for i in range(len(species)-1) if bits & (1<<i)}
        test=set(species)-train
        if len(train)<3 or len(test)<3:continue
        arms=count_by_arm(rows,train,test)
        tr=arms["train"];te=arms["evaluation"]
        if min(tr["numeric_literal_within"],te["numeric_literal_within"])<3:continue
        if min(tr["numeric_literal_between"],te["numeric_literal_between"])<2:continue
        strong+=1
        score=(min(tr["numeric_literal_between"],te["numeric_literal_between"]),
               min(tr["numeric_literal_dyads"],te["numeric_literal_dyads"]),
               tr["numeric_literal_dyads"]+te["numeric_literal_dyads"])
        if best is None or score>best[0]:
            best=(score,arms)
    return {"original_21_21_balanced_split_storage_support":checked,
            "literal_numeric_cell_gate_partitions":strong,
            "best_literal_numeric_both_types_supported_split":best[1] if best else None,
            "numeric_storage_does_not_validate_an_original_mating_experiment":True,
            "no_response_values_decoded":True}


def main():
    raw=fetch("source_data")
    identities=read_identity(raw,"source_data")
    result=source_mating_cell_storage(raw,identities)
    report=evaluate_support(result["rows"])
    print(json.dumps({
        "status":"SOURCE_ONLY_TIMEMA_SEXUAL_CELL_STORAGE_AND_SPECIES_HOLDOUT",
        "source_sha256":SOURCES["source_data"]["sha256"],
        "total_source_dyads":len(identities),
        "sexual_storage_class_counts":result["classes"],
        "article_reported_mating_pair_count":30,
        "source_numeric_literal_cell_count":result["classes"].get("numeric_literal_storage_not_authenticated_trials",0),
        "source_vs_article_mating_count_unresolved":True,
        "species_heldout_numeric_storage_support":report,
        "no_effect_size_or_trial_count_read":True,
    },sort_keys=True,indent=2))


if __name__=="__main__":
    main()
