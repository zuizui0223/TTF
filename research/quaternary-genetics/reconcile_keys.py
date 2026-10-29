"""Response-blind source-table key reconciliation; never reads He values."""
import collections, io, json, zipfile, xml.etree.ElementTree as ET
from pathlib import Path
from inspect_headers_only import fetch_verified
N="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def ix(ref):
    n=0
    for c in ref:
        if not c.isalpha(): break
        n=n*26+ord(c.upper())-64
    return n-1
def stream(blob, requested, presence=frozenset()):
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        strings=[]
        with z.open("xl/sharedStrings.xml") as f:
            for _,el in ET.iterparse(f,events=("end",)):
                if el.tag==N+"si":
                    strings.append("".join(x.text or "" for x in el.iter(N+"t")));el.clear()
        with z.open("xl/worksheets/sheet1.xml") as f:
            head=None
            for _,r in ET.iterparse(f,events=("end",)):
                if r.tag!=N+"row":continue
                out={}
                for c in r.findall(N+"c"):
                    i=ix(c.get("r",""))
                    if head is not None and head.get(i) not in requested and head.get(i) not in presence:continue
                    v=c.find(N+"v")
                    if head is not None and head.get(i) in presence:
                        out[i]=(v is not None and c.get("t")!="e")
                        continue
                    if v is None:continue
                    out[i]=strings[int(v.text)] if c.get("t")=="s" else v.text
                if head is None:head=out
                else:yield {name:out.get(i) for i,name in head.items() if name in requested or name in presence}
                r.clear()
def reconcile(pop,spec):
    meta={}
    for x in stream(spec,{"Spec_id","Study_id","N_pops"}):
        k=(str(x.get("Spec_id") or ""),str(x.get("Study_id") or ""))
        if k in meta:raise ValueError("nonunique metadata key")
        meta[k]=int(x["N_pops"])
    counts=collections.Counter()
    pop_ids=collections.defaultdict(list)
    for x in stream(pop,{"Spec_id","Study_id","Pop_id"}):
        k=(str(x.get("Spec_id") or ""),str(x.get("Study_id") or ""))
        counts[k]+=1
        if x.get("Pop_id") is not None:pop_ids[k].append(str(x["Pop_id"]))
    pm=set(counts);sm=set(meta)
    diffs=[(k,counts[k],meta[k],counts[k]-meta[k]) for k in pm&sm
           if counts[k]!=meta[k]]
    return {"status":"BLIND_KEY_RECONCILIATION",
        "population_rows":sum(counts.values()),
        "metadata_N_pops_sum":sum(meta.values()),
        "pop_group_count":len(pm),"metadata_group_count":len(sm),
        "groups_without_metadata":len(pm-sm),
        "groups_without_populations":len(sm-pm),
        "rows_without_metadata":sum(counts[k] for k in pm-sm),
        "orphan_groups":[{"species_id":k[0],"study_id":k[1],"rows":counts[k]} for k in sorted(pm-sm)],
        "expected_populations_in_unrepresented_metadata_groups":sum(meta[k] for k in sm-pm),
        "mismatched_joined_groups":len(diffs),
        "sum_positive_group_differences":sum(d for k,a,b,d in diffs if d>0),
        "sum_negative_group_differences":sum(d for k,a,b,d in diffs if d<0),
        "duplicate_Pop_id_within_group":sum(len(v)-len(set(v)) for v in pop_ids.values()),
        "duplicate_pop_id_groups":[{"species_id":k[0],"study_id":k[1],"duplicate_ids":len(v)-len(set(v))} for k,v in sorted(pop_ids.items()) if len(v)!=len(set(v))],
        "difference_histogram":dict(collections.Counter(d for k,a,b,d in diffs)),
        "top_mismatches":[{"species_id":k[0],"study_id":k[1],
           "population_rows":a,"N_pops":b,"difference":d}
           for k,a,b,d in sorted(diffs,key=lambda t:-abs(t[3]))[:15]],
        "genetic_values_opened":False}
if __name__=="__main__":
    import argparse
    a=argparse.ArgumentParser();a.add_argument("--metadata",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True);o=a.parse_args()
    m=json.loads(o.metadata.read_text());fs={x["filename"]:x for x in m["files"]}
    result=reconcile(fetch_verified(fs["pop_tab_v2025-03-31.xlsx"]),
                    fetch_verified(fs["spec_tab_v2025-03-31.xlsx"]))
    o.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
