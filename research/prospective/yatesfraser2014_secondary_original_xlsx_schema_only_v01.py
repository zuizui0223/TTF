"""Read only publisher-verified secondary XLSX sheet names and early header TEXT.

The Yates & Fraser 2014 Dryad 81k4v workbook may contain Raabova-supplied
meta-analysis data, not necessarily 2007 individual plant tags. The original
source must first match DOI + Zenodo preservation record + size + MD5.
Numeric cell values, formula cells, late rows, and response values are never
selected, decoded or exported. A lack of tag column in a short header window
is not evidence a 2007 raw tag table never existed.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import posixpath
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote,urlsplit
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"yatesfraser2014_raabova_secondary_xlsx_header_contract_v01.json"
NS={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "p":"http://schemas.openxmlformats.org/package/2006/relationships"}
def tag(local):return "{"+NS["m"]+"}"+local

def doi_norm(s):
    if not isinstance(s,str):return ""
    return re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:)","",s.strip().lower()).rstrip("/")

def verified_record_download(record:dict,c:dict)->str:
    org=c["original"]
    if str(record.get("id"))!=str(org["official_zenodo_record_id"]):
        raise ValueError("Not the original frozen Zenodo record")
    meta=record.get("metadata") or {}
    ids=[meta.get("doi"),record.get("doi")]
    ids += [v.get("identifier") for k in ("related_identifiers","identifiers",
                  "alternate_identifiers") for v in (meta.get(k) or []) if isinstance(v,dict)]
    pids=record.get("pids") or {}
    if isinstance(pids,dict) and isinstance(pids.get("doi"),dict):
        ids.append(pids["doi"].get("identifier"))
    if org["original_data_doi"] not in {doi_norm(v) for v in ids}:
        raise ValueError("No exact original published Dryad DOI")
    f=record.get("files") or []
    if isinstance(f,dict) and isinstance(f.get("entries"),dict):files=list(f["entries"].values())
    elif isinstance(f,list):files=f
    else:raise ValueError("Unexpected original preserved file list")
    selected=[x for x in files if isinstance(x,dict) and x.get("key")==org["original_filename"]]
    if len(selected)!=1:raise ValueError("Original Excel file absent or ambiguous")
    data=selected[0]
    if data.get("size")!=org["original_bytes"] or str(data.get("checksum","")).lower()!="md5:"+org["original_md5"]:
        raise ValueError("Original published workbook size/MD5 metadata changed")
    links=data.get("links") or {}
    link=links.get("content") or links.get("download") or links.get("self")
    if not isinstance(link,str):raise ValueError("Original file lacks preservation URL")
    uri=urlsplit(link);path=unquote(uri.path)
    allowed=(path.startswith("/api/files/") and path.endswith("/"+org["original_filename"])
         or path.startswith("/api/records/"+str(org["official_zenodo_record_id"])+"/files/")
            and path.endswith("/content"))
    if uri.scheme!="https" or uri.hostname!="zenodo.org" or not allowed:
        raise ValueError("Not an official Zenodo original file link")
    return link

def publisher_original_bytes(c:dict)->bytes:
    record=c["original"]
    url=record["official_record_api"]
    with urlopen(Request(url,headers={
      "Accept":"application/json","User-Agent":"TTF-YatesFraser-original-2014-XLSX-header-only/2026-11-16"
    }),timeout=30) as resp:
        b=resp.read(500001)
    if len(b)>500000:raise ValueError("Original preservation metadata too large")
    link=verified_record_download(json.loads(b),c)
    with urlopen(Request(link,headers={"Accept":"*/*",
             "User-Agent":"TTF-secondary-original-XLSX-sheetname-header-only/2026-11-16"}),timeout=60) as resp:
        raw=resp.read(record["original_bytes"]+1)
    if (len(raw)!=record["original_bytes"] or hashlib.md5(raw).hexdigest()!=record["original_md5"]
       or not raw.startswith(b"PK\x03\x04")):
        raise ValueError("Downloaded secondary original XLSX failed source checksum/signature")
    return raw

def sheet_links(z:zipfile.ZipFile,c:dict):
    book="xl/workbook.xml";rels="xl/_rels/workbook.xml.rels"
    if book not in z.namelist() or rels not in z.namelist():
        raise ValueError("Original XLSX workbook part missing")
    workbook=ET.fromstring(z.read(book))
    relationship=ET.fromstring(z.read(rels))
    targets={}
    for item in relationship:
        rid=item.get("Id");target=item.get("Target")
        if rid and target:
            target=target.lstrip("/") if target.startswith("/") else posixpath.normpath(posixpath.join("xl",target))
            if not target.startswith("xl/") or ".." in target.split("/"):
                raise ValueError("Unexpected XLSX workbook relation")
            targets[rid]=target
    out=[]
    for sheet in workbook.findall("m:sheets/m:sheet",NS):
        rid=sheet.get("{"+NS["r"]+"}id")
        if rid not in targets:raise ValueError("Missing original XLSX sheet relationship")
        name=sheet.get("name")
        out.append((name,targets[rid]))
    if not 1<=len(out)<=c["safe_XLSX_member_scope"]["max_worksheets"]:
        raise ValueError("Original worksheet count out of predeclared range")
    return out

def early_sheet_shared_strings(z:zipfile.ZipFile,path:str,max_rows:int,max_cols:int):
    early=[]
    with z.open(path) as stream:
        for event,el in ET.iterparse(stream,events=("end",)):
            if el.tag!=tag("row"):continue
            try:rownum=int(el.get("r","999999"))
            except ValueError:raise ValueError("Invalid worksheet row tag")
            if rownum>max_rows:break
            cells=[]
            for cell in el.findall("m:c",NS):
                rawref=cell.get("r","")
                match=re.match(r"^([A-Z]+)(\d+)$",rawref)
                if not match:continue
                col=0
                for letter in match.group(1):col=col*26+ord(letter)-64
                if col>max_cols:continue
                typ=cell.get("t","")
                if typ=="s":
                    value=cell.find("m:v",NS)
                    if value is not None and value.text and value.text.isdigit():
                        cells.append((col,"shared",int(value.text)))
                elif typ=="inlineStr":
                    rawtext="".join((x.text or "") for x in cell.findall(".//m:t",NS))
                    cells.append((col,"text",rawtext))
                # NO numeric, formula or date cell values are read.
            if cells:
                early.append((rownum,cells))
            el.clear()
    # Pick earliest plausible header row without inspecting later data values.
    acceptable=[x for x in early if len(x[1])>=2]
    return acceptable[0] if acceptable else (early[0] if early else None)

def early_string_values(z:zipfile.ZipFile,indices:set[int])->dict[int,str]:
    if not indices:return {}
    path="xl/sharedStrings.xml"
    if path not in z.namelist():raise ValueError("Header uses absent original SST")
    if max(indices)>10000:raise ValueError("Original source header refers to very deep string index")
    found={};index=0
    with z.open(path) as stream:
        for _,el in ET.iterparse(stream,events=("end",)):
            if el.tag!=tag("si"):continue
            if index in indices:
                found[index]="".join(x.text or "" for x in el.findall(".//m:t",NS))
            index+=1
            el.clear()
            if index>max(indices):break
    if set(found)!=indices:raise ValueError("Shared header string indices unresolvable")
    return found

def schema_only(xlsx:bytes,c:dict)->dict:
    lim=c["safe_XLSX_member_scope"]
    with zipfile.ZipFile(io.BytesIO(xlsx)) as archive:
        entries=archive.infolist()
        if len(entries)>lim["max_zip_members"] or sum(e.file_size for e in entries)>lim["max_uncompressed_total_bytes"]:
            raise ValueError("Original archive ZIP members exceed source-only limits")
        sheets=sheet_links(archive,c)
        early=[]
        wanted=set()
        for sheet,path in sheets:
            if not re.fullmatch(r"xl/worksheets/sheet\d+\.xml",path):
                raise ValueError("Unexpected nonworksheet XLSX sheet relation")
            if archive.getinfo(path).file_size>lim["max_worksheet_xml_bytes"]:
                raise ValueError("Original worksheet XML size exceeds cap")
            selected=early_sheet_shared_strings(archive,path,lim["first_text_rows_per_sheet"],
                                                lim["max_text_columns_scanned"])
            early.append((sheet,selected))
            if selected:
                wanted.update(v for col,mode,v in selected[1] if mode=="shared")
        strings=early_string_values(archive,wanted)
        out=[]
        for sheet,header in early:
            vals=[]
            rownum=None
            if header:
                rownum=header[0]
                for col,mode,value in header[1]:
                    text=strings[value] if mode=="shared" else value
                    if not text.strip():continue
                    if re.fullmatch(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?",text.strip()):
                        continue
                    # Do not expose original source data rows; report only plausible
                    # header-like strings. Never interpret numeric biological values.
                    vals.append({"column":col,"early_text_label":text[:120]})
            tokens={re.sub(r"[^a-z0-9]+","_",v["early_text_label"].lower()).strip("_") for v in vals}
            tag=sorted(tokens & set(c["candidate_2007_field_tag_keys"]))
            out.append({"sheet":sheet,
                "first_potential_text_header_row":rownum,
                "potential_header_text_labels":vals,
                "explicit_field_plant_tag_heading_candidates":tag,
                "numeric_fitness_effect_values_opened":False})
    return {"sheets":out,
            "any_explicit_plant_mark_candidate_header":any(x["explicit_field_plant_tag_heading_candidates"] for x in out)}

def run():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if c["status"]!="FROZEN_BEFORE_SECONDARY_ORIGINAL_XLSX_XML_SHEET_HEADER_EXTRACTION" or not c["safe_XLSX_member_scope"]["numeric_cells_or_formulas_NOT_SELECTED"]:
        raise ValueError("Data-bearing workbook access rule not frozen")
    raw=publisher_original_bytes(c)
    preview=schema_only(raw,c)
    return {
      "schema":"ttf_YatesFraser2014_original_Dryad_81k4v_Excel_early_text_headings_only_v01",
      "status":"DERIVATIVE_METAANALYSIS_ORIGINAL_SCHEMA_ONLY_NOT_RAABOVA_2007_TAG_CERTIFICATION",
      "publisher_source_DOI":c["original"]["original_data_doi"],
      "verified_official_Zenodo_record_id":c["original"]["official_zenodo_record_id"],
      "original_XLSX_bytes":len(raw),
      "original_XLSX_MD5_verified":c["original"]["original_md5"],
      "original_XLSX_SHA256":hashlib.sha256(raw).hexdigest(),
      "schema":preview,
      "2007_original_individual_field_tag_table_certified":False,
      "2007_vs_2014_original_individual_shared_count":None,
      "all_sheet_numeric_biological_response_values_unread":True,
      "original_secondary_meta_analysis_is_not_independent_primary_experiment":True,
      "previous_shared_2004_adult_cohort_and_44of50_stopped_gate_unchanged":True
    }

if __name__=="__main__":
    print("YATESFRASER2014_ORIGINAL_XLSX_HEADER_ONLY="+json.dumps(run(),sort_keys=True))
