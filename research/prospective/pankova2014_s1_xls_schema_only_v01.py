"""PLOS Pankova 2014 Table S1, legacy Excel header-only BIFF source audit.

Only official XLS workbook bytes are fetched; only text from the first
12 worksheet rows is decoded, numeric cell responses never accessed.
The parser is Python stdlib BIFF8 with generic OLE container reader.
If a shared string crosses a CONTINUE boundary, abstain.
"""
from __future__ import annotations
import argparse,hashlib,io,json,re,struct
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen

CONTRACT=Path(__file__).with_name("pankova2014_s1_xls_schema_only_contract_v01.json")
MAGIC=bytes.fromhex("D0CF11E0A1B11AE1")

def records(buf:bytes,start=0):
    pos=start
    while pos+4<=len(buf):
        sid,n=struct.unpack_from("<HH",buf,pos)
        stop=pos+4+n
        if stop>len(buf):raise ValueError("Truncated BIFF stream")
        yield sid,buf[pos+4:stop]
        pos=stop
        if sid==0x000A:break

def sheet_info(payload:bytes):
    if len(payload)<8:raise ValueError("Truncated BOUNDSHEET")
    offset=struct.unpack_from("<I",payload,0)[0]
    typ,count,flags=payload[5],payload[6],payload[7]
    width=2 if flags&1 else 1
    if len(payload)<8+count*width:raise ValueError("Truncated sheet name")
    name=payload[8:8+count*width].decode("utf-16le" if width==2 else "cp1252")
    return offset,typ,name

def sst_first_fragment_strings(payload:bytes,required_max:int):
    if required_max<0:return []
    if required_max>1000 or len(payload)<8:raise ValueError("Unsupported SST header indexes")
    unique=struct.unpack_from("<I",payload,4)[0]
    if required_max>=unique:raise ValueError("Header SST index outside original")
    values=[];pos=8
    for index in range(required_max+1):
        if pos+3>len(payload):raise ValueError("SST continuation needed; abstain")
        size,flags=struct.unpack_from("<HB",payload,pos);pos+=3
        runs=ext=0
        if flags&0x08:
            if pos+2>len(payload):raise ValueError("Rich header crosses CONTINUE")
            runs=struct.unpack_from("<H",payload,pos)[0];pos+=2
        if flags&0x04:
            if pos+4>len(payload):raise ValueError("Extended header crosses CONTINUE")
            ext=struct.unpack_from("<I",payload,pos)[0];pos+=4
        width=2 if flags&1 else 1
        nb=size*width
        if pos+nb+4*runs+ext>len(payload):
            raise ValueError("SST header string crosses CONTINUE; abstain")
        txt=payload[pos:pos+nb].decode("utf-16le" if width==2 else "cp1252")
        values.append(txt.strip())
        pos+=nb+4*runs+ext
    return values

def worksheet_early_labels(workbook:bytes,start:int,max_rows:int,max_cols:int):
    labels={};referenced=set()
    for kind,rec in records(workbook,start):
        if kind==0x00FD and len(rec)>=10:
            row,col=struct.unpack_from("<HH",rec,0)
            if row>=max_rows:break
            if col>=max_cols:continue
            idx=struct.unpack_from("<I",rec,6)[0]
            labels[row,col]=("sst",idx);referenced.add(idx)
        elif kind==0x0204 and len(rec)>=8:
            row,col=struct.unpack_from("<HH",rec,0)
            if row>=max_rows:break
            if col>=max_cols:continue
            n=struct.unpack_from("<H",rec,6)[0]
            if 8+n>len(rec):raise ValueError("BIFF LABEL text truncated")
            labels[row,col]=("inline",rec[8:8+n].decode("cp1252"))
        # No numeric/formula/value BIFF record is decoded or exported.
    return labels,referenced

def extract_early_text_headers(workbook:bytes,contract:dict)->dict:
    sheets=[];sst=None
    for kind,payload in records(workbook):
        if kind==0x0085:sheets.append(sheet_info(payload))
        if kind==0x00FC and sst is None:sst=payload
    if not 1<=len(sheets)<=contract["method"]["max_sheet_count"]:
        raise ValueError("Original worksheet count outside frozen limit")
    maxrow=contract["method"]["max_header_rows_per_sheet"]
    maxcol=contract["method"]["max_column_count"]
    rows=[];maxidx=-1
    for offset,typ,name in sheets:
        if typ!=0:
            rows.append((typ,name,{},set()));continue
        if offset>=len(workbook):
            raise ValueError("Invalid original worksheet stream offset")
        labels,indices=worksheet_early_labels(workbook,offset,maxrow,maxcol)
        if indices:maxidx=max(maxidx,max(indices))
        rows.append((typ,name,labels,indices))
    if maxidx>=0 and sst is None:raise ValueError("Missing original BIFF shared string table")
    stringtab=sst_first_fragment_strings(sst,maxidx) if maxidx>=0 else []
    output=[]
    for typ,name,labels,_ in rows:
        options={}
        for (row,col),(form,val) in labels.items():
            value=stringtab[val] if form=="sst" else val
            if value:options.setdefault(row,{})[col]=value
        sorted_rows=sorted(options.items(),key=lambda x:(-len(x[1]),x[0]))
        selected=sorted_rows[0] if sorted_rows else None
        chosenrow=selected[0]+1 if selected else None
        headings=[{"column_number":col+1,"column_name":value}
                  for col,value in sorted(selected[1].items())] if selected else []
        tokens={re.sub(r"[^a-z0-9]+","_",x["column_name"].lower()).strip("_") for x in headings}
        candidates=sorted(set(contract["method"]["candidate_plant_tag_column_names"])&tokens)
        output.append({"sheet_name":name,"original_sheet_type":typ,
          "chosen_early_header_row_1_based":chosenrow,
          "column_names_from_early_header":headings,
          "explicit_tag_id_field_name_candidates":candidates,
          "any_numeric_or_fitness_values_selected":False})
    return {"sheets":output,
            "any_explicit_tag_id_column_label_found":any(x["explicit_tag_id_field_name_candidates"] for x in output)}

def original_workbook_schema()->dict:
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (c["status"]!="FROZEN_AFTER_512BYTE_ORIGINAL_OLE_HEADER_BEFORE_ANY_WORKBOOK_ROWS"
       or c["method"]["max_header_rows_per_sheet"]!=12):
        raise ValueError("No unmodified original pre-data schema contract")
    url=c["original"]["official_url"]
    if urlparse(url).hostname!="journals.plos.org":
        raise ValueError("Source URL not original PLOS domain")
    with urlopen(Request(url,headers={"Accept":"*/*",
           "User-Agent":"TTF-public-original-PLOS-XLS-only-sheets-and-headings/2026-11-15"}),timeout=60) as r:
        raw=r.read(c["original"]["max_original_file_bytes"]+1)
    if len(raw)>c["original"]["max_original_file_bytes"] or not raw.startswith(MAGIC):
        raise ValueError("Published XLS unexpectedly large or not original OLE2")
    import olefile
    ole=olefile.OleFileIO(io.BytesIO(raw))
    try:
        streams={"/".join(x) for x in ole.listdir()}
        name="Workbook" if "Workbook" in streams else "Book" if "Book" in streams else None
        if name is None:raise ValueError("Published XLS OLE workbook stream not found")
        wbook=ole.openstream(name).read()
    finally:
        ole.close()
    m=extract_early_text_headers(wbook,c)
    return {
        "schema":"ttf_Pankova2014_original_PLOS_S1_XLS_early_worksheet_header_metadata_v01",
        "original_paper_doi":c["original"]["doi"],
        "source_url":url,
        "downloaded_original_size_bytes":len(raw),
        "first_received_XLS_SHA256":hashlib.sha256(raw).hexdigest(),
        "storage_format":"OLE2_BIFF8",
        "original_worksheet_count":len(m["sheets"]),
        "worksheets":m["sheets"],
        "explicit_individual_tag_column_candidate":m["any_explicit_tag_id_column_label_found"],
        "status":"EARLY_WORKSHEET_LABELS_ONLY_NO_INDIVIDUAL_DATA_CROSSJOIN",
        "all_individual_tag_ID_values_selected":False,
        "biological_fitness_measurement_cell_values_selected":False,
        "2007_individual_raw_tag_identifiers_found":False,
        "crosspublication_individual_ID_matching_completed":False,
        "author_confirmed_one_shared_2004_adult_cohort_unchanged":True,
        "Bontrager_44_of_50_frozen_gate_unchanged":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    out=original_workbook_schema()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("PANKOVA2014_PUBLISHED_XLS_SCHEMA_ONLY="+json.dumps({
      "bytes":out["downloaded_original_size_bytes"],
      "sha256":out["first_received_XLS_SHA256"],
      "explicit_individual_tag_column_candidate":out["explicit_individual_tag_column_candidate"],
      "worksheets":out["worksheets"],
      "crosspaper_ID_join":out["crosspublication_individual_ID_matching_completed"]
    },ensure_ascii=False,sort_keys=True))
if __name__=="__main__":main()
