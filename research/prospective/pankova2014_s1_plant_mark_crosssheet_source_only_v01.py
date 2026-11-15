"""Response-blind original Pankova2014 individual marker crosssheet metadata audit.

Only BIFF numeric/text cells in three predeclared columns (St,So,mark) in
survival and performance sheets are selected. NEVER index/read biological
surv2009, leaves, biomass, colonization or their numeric values. No raw
individual marker IDs or individual observation rows exported.

This checks 2014 S1 within-workbook plant ID consistency only. It cannot
establish individual tag overlap with 2007 without a 2007 individual source.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
import struct
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen

from pankova2014_s1_xls_schema_only_v01 import (
    MAGIC, records, sheet_info, sst_first_fragment_strings,
    extract_early_text_headers,
)

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"pankova2014_s1_original_plant_mark_join_contract_v01.json"
HEADER_CONTRACT=HERE/"pankova2014_s1_xls_schema_only_contract_v01.json"

def rk_number(bits:int) -> float:
    if bits & 2:
        n=bits>>2
        if n & (1<<29):n-=(1<<30)
        out=float(n)
    else:
        out=struct.unpack("<d",struct.pack("<Q",(bits & 0xFFFFFFFC)<<32))[0]
    return out/100 if bits & 1 else out

def scalar_id(x:object) -> str|None:
    if isinstance(x,str):
        s=x.strip()
        return s if s else None
    if isinstance(x,(int,float)) and not isinstance(x,bool):
        v=float(x)
        if not math.isfinite(v):return None
        if v.is_integer():return str(int(v))
    return None

def file_bytes(contract:dict) -> bytes:
    meta=contract["original"]
    url=meta["original_source_url"]
    if (urlparse(url).hostname!="journals.plos.org"
        or meta["paper_doi"]!="10.1371/journal.pone.0093967"):
        raise ValueError("Not the official original PLOS S1 source")
    req=Request(url,headers={"Accept":"*/*",
             "User-Agent":"TTF-original-2014-S1-plant-ID-only/2026-11-15"})
    with urlopen(req,timeout=60) as r:
        raw=r.read(meta["bytes"]+1)
    if (len(raw)!=meta["bytes"] or not raw.startswith(MAGIC)
        or hashlib.sha256(raw).hexdigest()!=meta["sha256"]):
        raise ValueError("Original XLS bytes do not match frozen published source")
    return raw

def workbook_stream(raw:bytes)->bytes:
    import olefile
    ole=olefile.OleFileIO(io.BytesIO(raw))
    try:
        names={"/".join(x) for x in ole.listdir()}
        name="Workbook" if "Workbook" in names else "Book" if "Book" in names else None
        if name is None:raise ValueError("No original OLE workbook binary")
        return ole.openstream(name).read()
    finally:
        ole.close()

def biff_sheet_locations(buf:bytes):
    sheets=[];sst=None
    for kind,rec in records(buf):
        if kind==0x0085:sheets.append(sheet_info(rec))
        if kind==0x00FC and sst is None:sst=rec
    if sst is None:raise ValueError("Original SST unavailable")
    return sheets,sst

def read_selected_id_cells(buf:bytes,offset:int,cols:dict[str,int],nlimit:int):
    if not set(cols)=={"St","So","mark"}:
        raise ValueError("Source-only individual marker metadata allowlist drift")
    reverse={col:key for key,col in cols.items()}
    rows=defaultdict(dict)
    sstrefs={}
    for kind,rec in records(buf,offset):
        if kind==0x00FD and len(rec)>=10: # LABELSST
            row,col=struct.unpack_from("<HH",rec,0)
            if row<1 or col not in reverse:continue
            if row>nlimit:raise ValueError("Metadata source rows outside frozen limit")
            idx=struct.unpack_from("<I",rec,6)[0]
            sstrefs[row,col]=idx
        elif kind==0x0204 and len(rec)>=8:
            row,col=struct.unpack_from("<HH",rec,0)
            if row<1 or col not in reverse:continue
            if row>nlimit:raise ValueError("Original source row count exceeds limit")
            n=struct.unpack_from("<H",rec,6)[0]
            if 8+n>len(rec):raise ValueError("Original direct LABEL truncated")
            rows[row][reverse[col]]=rec[8:8+n].decode("cp1252",errors="strict")
        elif kind==0x0203 and len(rec)>=14: # NUMBER
            row,col=struct.unpack_from("<HH",rec,0)
            if row<1 or col not in reverse:continue
            if row>nlimit:raise ValueError("Original source row count exceeds limit")
            rows[row][reverse[col]]=struct.unpack_from("<d",rec,6)[0]
        elif kind==0x027E and len(rec)>=10: # RK
            row,col=struct.unpack_from("<HH",rec,0)
            if row<1 or col not in reverse:continue
            if row>nlimit:raise ValueError("Original source row count exceeds limit")
            rows[row][reverse[col]]=rk_number(struct.unpack_from("<I",rec,6)[0])
        elif kind==0x00BD and len(rec)>=12: # MULRK
            row,first=struct.unpack_from("<HH",rec,0)
            last=struct.unpack_from("<H",rec,len(rec)-2)[0]
            if len(rec)!=6+(last-first+1)*6:raise ValueError("Original MULRK width drift")
            if row<1:continue
            if row>nlimit and any(c in reverse for c in range(first,last+1)):
                raise ValueError("Original source row count exceeds limit")
            for col in range(first,last+1):
                if col in reverse:
                    pos=4+(col-first)*6+2
                    rows[row][reverse[col]]=rk_number(struct.unpack_from("<I",rec,pos)[0])
        elif kind==0x0006 and len(rec)>=6:
            row,col=struct.unpack_from("<HH",rec,0)
            if row>=1 and col in reverse:
                raise ValueError("Source ID metadata column unexpectedly formula-based; abstain")
    return rows,sstrefs

def per_sheet_id_metadata(selected:dict,contract:dict)->dict:
    intended=contract["allow_row_cell_values_only_for_columns"]
    if intended!=["mark","St","So"]:
        raise ValueError("ID metadata fields changed after predata freeze")
    out={}
    keys={}
    for sheet in ("survival","performance"):
        rows=selected[sheet]
        composite=[];by_id=Counter();by_site=Counter()
        missing=Counter()
        for vals in rows.values():
            mark=scalar_id(vals.get("mark"))
            site=scalar_id(vals.get("St"))
            origin=scalar_id(vals.get("So"))
            if mark is None:
                missing["blank_or_invalid_mark"]+=1
            else:
                by_id[mark]+=1
            if site is None:missing["blank_or_invalid_target_site"]+=1
            if origin is None:missing["blank_or_invalid_origin_site"]+=1
            if mark is not None and site is not None:
                composite.append((site,mark,origin))
                by_site[site]+=1
        sitekeys={(s,m) for s,m,o in composite}
        pairs=Counter((s,m) for s,m,o in composite)
        origin_per_key=defaultdict(set)
        for s,m,o in composite:origin_per_key[s,m].add(o)
        keys[sheet]={"sitekeys":sitekeys,"origins":origin_per_key,
                     "full":{(s,m,o) for s,m,o in composite if o is not None}}
        out[sheet]={
            "metadata_row_count_with_any_selected_field":len(rows),
            "rows_with_target_site_and_mark":len(composite),
            "unique_mark_text_or_number_tokens":len(by_id),
            "distinct_site_by_mark_keys":len(sitekeys),
            "duplicate_site_by_mark_row_count":sum(n-1 for n in pairs.values() if n>1),
            "site_by_mark_keys_with_conflicting_origin_labels":sum(len(v)>1 for v in origin_per_key.values()),
            "unique_target_site_labels":len(by_site),
            "missingness":dict(sorted(missing.items())),
            "no_original_plant_mark_or_site_values_exported":True,
            "biological_numerical_values_read":False
        }
    a=keys["survival"];b=keys["performance"]
    common=a["sitekeys"]&b["sitekeys"]
    conflict=sum(a["origins"][k]!=b["origins"][k] for k in common)
    return {
        "per_sheet":out,
        "within_2014_crosssheet":{
            "shared_site_by_mark_identifiers":len(common),
            "performance_IDs_not_in_survival":len(b["sitekeys"]-a["sitekeys"]),
            "survival_IDs_not_in_performance":len(a["sitekeys"]-b["sitekeys"]),
            "matched_site_by_mark_with_disagreeing_origin_ID_sets":conflict,
            "same_site_mark_is_unique_in_both_sheets":all(out[s]["duplicate_site_by_mark_row_count"]==0 for s in out),
            "crosssheet_individual_2014_metadata_link_present":len(common)>0 and conflict==0,
            "crosspaper_2007_2014_individual_number_overlap_quantified":False
        }
    }

def audit_original_XLS()->dict:
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    hdr=json.loads(HEADER_CONTRACT.read_text(encoding="utf-8"))
    if (c["status"]!="FROZEN_AFTER_EXPLANATION_CODEBOOK_BEFORE_NUMERIC_PLANT_MARK_METADATA_CELLS"
        or c["allow_row_cell_values_only_for_columns"]!=["mark","St","So"]
        or c["original"]["sha256"]!="7a8de0c203aec2b2edf0a7b6bceea649663a17095471a475fdb8a0abd69665c5"):
        raise ValueError("Original source / allowed metadata fields not frozen")
    raw=file_bytes(c)
    wb=workbook_stream(raw)
    meta=extract_early_text_headers(wb,hdr)
    sheetmap={x["sheet_name"]:x for x in meta["sheets"]}
    if set(sheetmap)!={"survival","performance","explanation"}:
        raise ValueError("Unexpected published 2014 XLS sheet names")
    loc,sst=biff_sheet_locations(wb)
    offset={name:(pos,typ) for pos,typ,name in loc}
    stash={}
    selected={}
    all_ids=set()
    for sheet in c["allow_selected_original_sheets"]:
        orig_header=c["expected_sheet_headers"][sheet]
        cols={s["column_name"]:s["column_number"]-1 for s in sheetmap[sheet]["column_names_from_early_header"]}
        if list(cols)!=orig_header:
            raise ValueError("Original 2014 source sheet schema drifted")
        required={k:cols[k] for k in ("St","So","mark")}
        if offset[sheet][1]!=0:raise ValueError("Not an original XLS worksheet")
        rows,refs=read_selected_id_cells(wb,offset[sheet][0],required,c["row_limits"]["per_sheet_max"])
        stash[sheet]=(rows,refs,required)
        all_ids.update(refs.values())
    if all_ids:
        # Parse only shared-string prefix needed by allowlisted columns.
        strings=sst_first_fragment_strings(sst,max(all_ids))
    else:strings=[]
    for sheet,(rows,refs,cols) in stash.items():
        for (row,col),index in refs.items():
            rows[row][next(key for key,val in cols.items() if val==col)]=strings[index]
        selected[sheet]=rows
    result=per_sheet_id_metadata(selected,c)
    return {
        "schema":"ttf_Pankova2014_original_S1_individual_mark_within2014_crosssheet_metadata_v01",
        "status":"2014_MARK_SITE_METADATA_LINKAGE_ONLY_NOT_2007_INDIVIDUAL_PROOF",
        "original_S1_DOI":c["original"]["paper_doi"],
        "original_XLS_SHA256":c["original"]["sha256"],
        "original_header_explanation_mark":c["verified_dictionary"]["mark"],
        "source_field_allowlist":c["allow_row_cell_values_only_for_columns"],
        "summary":result,
        "raw_individual_mark_IDs_or_garden_or_origin_values_exported":False,
        "survival_or_performance_biological_values_selected":False,
        "publisher_2007_raw_individual_tag_table_located":False,
        "two_paper_adult_cohort_author_confirmed":True,
        "individual_level_2007_vs_2014_shared_plant_count":None,
        "independent_2004_adult_experiments_added_by_2014_paper":0,
        "frozen_PR154_Bontrager_44_of_50_source_stop_unchanged":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    out=audit_original_XLS()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("PANKOVA2014_ORIGINAL_PLANT_MARK_CROSSSHEET_SOURCE_ONLY="+json.dumps({
        "sheets":out["summary"]["per_sheet"],
        "within2014":out["summary"]["within_2014_crosssheet"],
        "crosspaper_individual_count":out["individual_level_2007_vs_2014_shared_plant_count"]
    },ensure_ascii=False,sort_keys=True))
if __name__=="__main__":main()
