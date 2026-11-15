"""Exact-SHA original Dryad CSV HEADER and ZIP member names, not fitness rows."""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from pathlib import Path
from urllib.request import Request,urlopen

ROOT="https://datadryad.org"

def sha256_hex(data):
    return hashlib.sha256(data).hexdigest()

def classify_header(headers):
    words=[str(x).lower().replace(".","_").replace("-","_") for x in headers]
    def coords(prefixes):
        return bool(
           any("lat" in x and any(p in x for p in prefixes) for x in words)
           and any(("lon" in x or "long" in x) and any(p in x for p in prefixes) for x in words))
    return {
       "possible_origin_latlon_fields":coords(("origin","source","pop","provenance","seed")),
       "possible_garden_latlon_fields":coords(("garden","site","transplant","destination","common")),
       "has_species_label":any("species" in s or "taxon" in s or "spname" in s for s in words),
       "has_study_label":any("study" in s or "paper" in s or "ref" in s or "citation" in s for s in words),
       "interpretation":"HEURISTIC_HEADER_NAMES_ONLY_NOT_CONFIRMED_ORIGIN_GARDEN_COORDINATES"
    }

def inspect_verified_source(raw,kind):
    if kind=="csv":
        first=raw.split(b"\n",1)[0]
        if len(first)>50000:raise ValueError("Header longer than safety limit")
        fields=next(csv.reader([first.decode("utf-8-sig").rstrip("\r")],strict=True))
        if not fields or len(fields)>300 or len(set(fields))!=len(fields):
            raise ValueError("Invalid or duplicate CSV source headers")
        return {"source_format":"csv","headers":fields,"source_column_count":len(fields),
                "header_only_content_inspection":True,"row_2_or_later_fields_parsed":False,
                "fitted_response_or_location_values_inspected":False}
    if kind=="zip":
        if not zipfile.is_zipfile(io.BytesIO(raw)):
            raise ValueError("Published ZIP invalid")
        with zipfile.ZipFile(io.BytesIO(raw)) as f:names=f.namelist()
        if len(names)>5000:raise ValueError("ZIP directory listing larger than limit")
        return {"source_format":"zip","member_names":names,"member_count":len(names),
                "archive_member_contents_read":False,"header_only_content_inspection":True,
                "row_2_or_later_fields_parsed":False,
                "fitted_response_or_location_values_inspected":False}
    raise ValueError("Only locked CSV or ZIP source schema types supported")

def original_download(file_id,max_bytes):
    if not isinstance(file_id,int) or file_id<=0:
        raise ValueError("Invalid frozen file ID")
    req=Request(f"{ROOT}/api/v2/files/{file_id}/download",headers={
        "User-Agent":"TTF-original-transplant-source-headers-only-v01"})
    h=hashlib.sha256()
    data=[];size=0
    with urlopen(req,timeout=80) as response:
        while True:
            piece=response.read(1048576)
            if not piece:break
            size+=len(piece)
            if size>max_bytes:raise ValueError("Original file larger than frozen cap")
            h.update(piece)
            data.append(piece)
    return b"".join(data),h.hexdigest()

def run(contract):
    if contract.get("status")!="FROZEN_AFTER_OFFICIAL_FILE_CATALOG_BEFORE_DATA_CONTENT":
        raise ValueError("Header-only source contract not prospectively frozen")
    reports=[]
    for source in contract["selected_source_files"]:
        result={"candidate":source["candidate"],"file_id":source["file_id"],
                "file_name":source["file_name"],"original_version_id":source["version_id"],
                "source_metadata_rows_only":True,"biological_effect_values_accessed":False}
        try:
            raw,digest=original_download(source["file_id"],source["max_download_bytes"])
            if len(raw)!=source["size_bytes"] or digest!=source["sha256"]:
                raise ValueError("Original Dryad file size or SHA256 differs from catalog")
            result.update(inspect_verified_source(raw,source["kind"]))
            if source["kind"]=="csv":
                result["column_name_inference"]=classify_header(result["headers"])
            result["sha256_original_verified"]=digest
            result["original_source_byte_count"]=len(raw)
            result["status"]=("SHA_VERIFIED_SOURCE_HEADER_ONLY" if source["kind"]=="csv"
                              else "SHA_VERIFIED_ARCHIVE_MEMBER_DIRECTORY_ONLY")
        except Exception as exc:
            result.update({"status":"SOURCE_SCHEMA_UNAVAILABLE_NOT_BIOLOGICAL_NULL",
                           "error_type":type(exc).__name__,"error_message":str(exc)[:250]})
        reports.append(result)
        print("TRANSPLANT_SOURCE_SCHEMA",source["candidate"],result["status"],flush=True)
    return {"schema":"ttf_Dryad_new_transplant_original_file_headers_only_v01",
            "status":"SOURCE_HEADER_AUDIT_DONE_WITH_INDIVIDUAL_STOP_FLAGS",
            "source_file_reports":reports,
            "original_TTF_Haldane_outcome_still_closed":True,
            "no_fitness_or_plasticity_response_values_analyzed":True,
            "source_study_overlap_not_resolved":True,
            "new_phenotype_or_fitness_prediction_not_authorized":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    report=run(json.loads(args.contract.read_text()))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
