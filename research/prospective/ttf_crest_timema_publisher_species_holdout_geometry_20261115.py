"""Timema 2024 publisher source-only species/endpoint geometry gate.

Reads ONLY pair/species/population identity columns from two SHA256-locked XLSX
workbooks. Outcome columns (habitat/sexual RI and nucDNA/GBS magnitudes) are
never decoded, stored or reported. No species-heldout biological fit.
"""
from __future__ import annotations
import collections
import hashlib
import io
import json
import re
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from posixpath import normpath

URL="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-46294-9/MediaObjects/"
SOURCES={
    "source_data":{
      "filename":"41467_2024_46294_MOESM6_ESM.xlsx",
      "sha256":"c309402e1b1fa58ab0fb0763086a155aff028cdaa04b804c2c91796c5bed61e2",
      "sheet":"Figure3", "header_row":1,
    },
    "supplement_1":{
      "filename":"41467_2024_46294_MOESM4_ESM.xlsx",
      "sha256":"7e2b6210c182859d7b4d77ad8b0ca0858a1e437c0106a9c42e6c3787c24a33db",
      "sheet":"Supplementary Data 1", "header_row":2,
    },
}
M="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
PKG="http://schemas.openxmlformats.org/package/2006/relationships"
REL="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
IDENTITY_COLS={"A":"pair_id","C":"species_1","D":"species_2","E":"population_pair"}
HEADER={
  "pair_id":re.compile(r"^Pair\s*$",re.I),
  "species_1":re.compile(r"^Species\s*1\s*$",re.I),
  "species_2":re.compile(r"^Species\s*2\s*$",re.I),
  "population_pair":re.compile(r"^Populations?\s*$",re.I),
}


def fetch(name:str)->bytes:
    p=SOURCES[name]
    request=urllib.request.Request(URL+p["filename"],headers={"User-Agent":"TTF-source-identity/1.0"})
    with urllib.request.urlopen(request,timeout=45) as resp:
        raw=resp.read(1_000_001)
    if hashlib.sha256(raw).hexdigest()!=p["sha256"]:
        raise ValueError(f"{name}: exact XLSX SHA256 drift")
    return raw


def read_identity(raw:bytes,source_name:str)->list[dict]:
    p=SOURCES[source_name]
    if hashlib.sha256(raw).hexdigest()!=p["sha256"]:
        raise ValueError("original publisher XLSX SHA256 mismatch")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names=set(z.namelist())
        book=ET.fromstring(z.read("xl/workbook.xml"))
        rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap={r.get("Id"):r.get("Target") for r in rels.findall(f"{{{PKG}}}Relationship")}
        sst=[]
        if "xl/sharedStrings.xml" in names:
            ss=ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in ss.findall(f"{{{M}}}si"):
                sst.append("".join(t.text or "" for t in si.findall(f".//{{{M}}}t")))
        wanted=next((v for v in book.findall(f".//{{{M}}}sheets/{{{M}}}sheet")
                     if v.get("name")==p["sheet"]),None)
        if wanted is None:raise ValueError("expected source worksheet absent")
        target=relmap.get(wanted.get(f"{{{REL}}}id"))
        if not target:raise ValueError("source sheet relationship missing")
        path=(target.lstrip("/") if target.startswith("/xl/") else normpath("xl/"+target))
        if path not in names:raise ValueError("source sheet body absent")
        tree=ET.fromstring(z.read(path))
        sheetdata=tree.find(f"{{{M}}}sheetData")
        if sheetdata is None:raise ValueError("source sheet missing rows")
        rows=sheetdata.findall(f"{{{M}}}row")
        column_letters=re.compile(r"^([A-Z]+)\d+$")
        def identity_cells(row):
            out={}
            for cell in row.findall(f"{{{M}}}c"):
                match=column_letters.match(cell.get("r",""))
                if match is None or match.group(1) not in IDENTITY_COLS:continue
                col=match.group(1)
                if col in out:raise ValueError("duplicate identity column")
                if cell.get("t")=="s":
                    ix=cell.find(f"{{{M}}}v")
                    if ix is None or ix.text is None: continue
                    val=sst[int(ix.text)]
                elif cell.get("t")=="inlineStr":
                    val="".join(t.text or "" for t in cell.findall(f".//{{{M}}}t"))
                else:
                    # Numeric *pair ID* is an allowed identity, but no other
                    # numeric trait or response cell is inspected.
                    if col!="A":raise ValueError("nontext species/population identity column")
                    v=cell.find(f"{{{M}}}v")
                    val=v.text if v is not None and v.text else ""
                out[IDENTITY_COLS[col]]=val.strip()
            return out
        headers=[r for r in rows if int(r.get("r","0"))==p["header_row"]]
        if len(headers)!=1:raise ValueError("source header row missing or duplicate")
        hd=identity_cells(headers[0])
        if set(hd)!=set(IDENTITY_COLS.values()):
            raise ValueError(f"source identity header columns missing: {sorted(set(IDENTITY_COLS.values())-set(hd))}")
        for k,patt in HEADER.items():
            if not patt.fullmatch(hd[k]):
                raise ValueError(f"source header mismatch for {k}: {hd[k]!r}")
        data=[]
        for row in rows:
            if int(row.get("r","0"))<=p["header_row"]:continue
            info=identity_cells(row)
            if all(not info.get(k) for k in IDENTITY_COLS.values()):continue
            data.append(info)
    if len(data)!=42:raise ValueError(f"{source_name}: 42 source dyads expected, got {len(data)}")
    if len({d["pair_id"] for d in data})!=42:raise ValueError("duplicate pair identifiers")
    if any(not d.get(k) for d in data for k in IDENTITY_COLS.values()):
        raise ValueError("missing species or population identity source cell")
    return data


def split_geometry(rows:list[dict])->dict:
    species=sorted({x for row in rows for x in (row["species_1"],row["species_2"])})
    if len(species)>14 or len(species)<4:raise ValueError("unexpected number of source species")
    # Source species-interaction graph before any reproductive response is read.
    adjacency={sp:set() for sp in species}
    for row in rows:
        a,b=row["species_1"],row["species_2"]
        if a!=b:
            adjacency[a].add(b)
            adjacency[b].add(a)
    pending=set(species)
    components=[]
    while pending:
        start=min(pending)
        seen={start}
        todo=[start]
        while todo:
            node=todo.pop()
            for other in sorted(adjacency[node]):
                if other not in seen:
                    seen.add(other)
                    todo.append(other)
        components.append(sorted(seen))
        pending-=seen
    components.sort(key=lambda xs:(-len(xs),xs))
    comp_counts=[sum(r["species_1"] in group and r["species_2"] in group for r in rows)
                 for group in components]
    # Canonical species partitions: first species anchored to training.
    n=len(species)
    best=None
    viable=0
    for mask in range(1 << (n-1)):
        train={species[0]}|{species[i+1] for i in range(n-1) if mask & (1<<i)}
        test=set(species)-train
        if len(train)<2 or len(test)<2:continue
        tr=[r for r in rows if r["species_1"] in train and r["species_2"] in train]
        te=[r for r in rows if r["species_1"] in test and r["species_2"] in test]
        dropped=len(rows)-len(tr)-len(te)
        if tr and te:viable+=1
        # Deterministic geometry-only maximin score; NO response values.
        score=(min(len(tr),len(te)),min(len(train),len(test)),
               -abs(len(tr)-len(te)),len(tr)+len(te))
        if best is None or score>best[0]:
            best=(score,{"train_species":sorted(train),"test_species":sorted(test),
                         "train_dyads":len(tr),"test_dyads":len(te),
                         "cross_boundary_dyads_excluded":dropped})
    return {
       "species_count":n,
       "species_labels":species,
       "within_species_pair_count":sum(r["species_1"]==r["species_2"] for r in rows),
       "between_species_pair_count":sum(r["species_1"]!=r["species_2"] for r in rows),
       "per_species_endpoint_incidence":dict(sorted(collections.Counter(
           s for row in rows for s in (row["species_1"],row["species_2"])
       ).items())),
       "source_species_graph_components":components,
       "source_species_graph_component_pair_counts":comp_counts,
       "source_species_graph_component_count":len(components),
       "per_species_within_pair_counts":{
          sp:sum(r["species_1"]==sp and r["species_2"]==sp for r in rows)
          for sp in species
       },
       "partitions_with_nonempty_species_disjoint_train_and_test":viable,
       "best_balanced_species_split_source_only":best[1] if best else None,
       "comparisons_are_not_independent_evolutionary_replicates":True,
    }


def main():
    original=read_identity(fetch("source_data"),"source_data")
    supplement=read_identity(fetch("supplement_1"),"supplement_1")
    o={r["pair_id"]:r for r in original}
    s={r["pair_id"]:r for r in supplement}
    if set(o)!=set(s):
        raise ValueError("publisher source/supplement pair identifier mismatch")
    mismatches=[]
    for pid in sorted(o):
        if any(o[pid][field]!=s[pid][field] for field in ("species_1","species_2","population_pair")):
            mismatches.append(pid)
    out={"status":"ORIGINAL_PUBLISHER_SPECIES_GEOMETRY_NO_OUTCOMES",
         "sources":{key:{"filename":p["filename"],"sha256":p["sha256"],
                         "sheet":p["sheet"],"pair_count":42}
                   for key,p in SOURCES.items()},
         "cross_source_identity_mismatch_pair_ids":mismatches,
         "cross_source_identity_consistent":not bool(mismatches),
         "geometry":split_geometry(original),
         "published_identity_examples_first_5":original[:5],
         "n_literal_species_labels_vs_article_nine_species":"TAXON_CONCEPT_RECONCILIATION_PENDING",
         "pairs_compared":42,
         "original_response_columns_read":False,
         "no_TTF_transfer_effect_claim":True}
    print(json.dumps(out,sort_keys=True,indent=2,ensure_ascii=False))


if __name__=="__main__": main()
