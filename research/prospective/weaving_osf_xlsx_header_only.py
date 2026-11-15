"""Check exact original 2022 OSF XLSX file identity and sheet HEADERS ONLY.

Never inspect data row 2+, acclimation response ratios, experimental thermal
limits, species-specific response values, or source coordinates. All bytes
downloaded in an ephemeral runner; public output contains only file/sheet
structure and source SHA.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlparse
import xml.etree.ElementTree as ET
import zipfile

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"weaving_2022_fixed_source_schema_contract_v01.json"
NS_SS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
NS_REL="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
NS_PACKAGE="{http://schemas.openxmlformats.org/package/2006/relationships}"


def get_public_osf_file(file_id):
    if not re.fullmatch("[a-f0-9]{24}",file_id):
        raise ValueError("Invalid OSF source file ID")
    api=f"https://api.osf.io/v2/files/{file_id}/"
    request=Request(api,headers={"Accept":"application/vnd.api+json",
                                 "User-Agent":"TTF-OSF-exact-source-header-only/0.1"})
    with urlopen(request,timeout=40) as response:
        obj=json.load(response)
    data=obj.get("data") or {}
    if data.get("id")!=file_id:
        raise ValueError("OSF source ID differs from frozen metadata")
    url=(data.get("links") or {}).get("download")
    if not url or urlparse(url).scheme!="https":
        raise ValueError("OSF did not return an HTTPS source download link")
    req=Request(url,headers={"User-Agent":"TTF-OSF-exact-source-header-only/0.1"})
    with urlopen(req,timeout=100) as response:
        # Max 2MiB; original source is 351661 B.
        raw=response.read(2*1024*1024)
    return raw


def shared_strings(z):
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    strings=[]
    with z.open("xl/sharedStrings.xml") as f:
        for _,el in ET.iterparse(f,events=("end",)):
            if el.tag==NS_SS+"si":
                strings.append("".join(x.text or "" for x in el.iter(NS_SS+"t")))
                el.clear()
    return strings


def sheet_descriptors(z):
    workbook=ET.fromstring(z.read("xl/workbook.xml"))
    rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    lookup={a.get("Id"):a.get("Target") for a in rel.findall(NS_PACKAGE+"Relationship")}
    out=[]
    for element in workbook.iter(NS_SS+"sheet"):
        rid=element.get(NS_REL+"id")
        target=lookup.get(rid) or ""
        if target.startswith("/xl/"):
            sheetpath=target.lstrip("/")
        elif target.startswith("worksheets/"):
            sheetpath="xl/"+target
        else:
            raise ValueError("Unsupported sheet relationship target")
        if not sheetpath.startswith("xl/worksheets/") or sheetpath not in z.namelist():
            raise ValueError("Original worksheet ZIP part unavailable")
        out.append((element.get("name") or "",sheetpath))
    return out


def read_only_first_row(z,path,shared):
    with z.open(path) as f:
        dim=None
        header=None
        for event,el in ET.iterparse(f,events=("start","end")):
            if event=="start" and el.tag==NS_SS+"dimension":
                dim=el.get("ref")
            if event=="end" and el.tag==NS_SS+"row":
                row_index=el.get("r")
                if row_index!="1":
                    # Never access phenotype data rows.
                    break
                cells=[]
                for cell in el.findall(NS_SS+"c"):
                    ref=cell.get("r") or ""
                    col=re.match(r"[A-Z]+",ref)
                    if not col:
                        continue
                    value=cell.find(NS_SS+"v")
                    if cell.get("t")=="s" and value is not None:
                        ix=int(value.text)
                        literal=shared[ix] if 0<=ix<len(shared) else ""
                    elif cell.get("t")=="inlineStr":
                        literal="".join(x.text or "" for x in cell.iter(NS_SS+"t"))
                    else:
                        literal=value.text if value is not None else ""
                    cells.append({"column":col.group(0),"header":literal})
                header=cells
                break
            if event=="end" and el.tag not in (NS_SS+"c",NS_SS+"v",NS_SS+"is",NS_SS+"t"):
                el.clear()
    names=[c["header"] for c in (header or [])]
    normalized=[re.sub(r"[^a-z]","",str(a).lower()) for a in names]
    contains_lat=any(x in {"latitude","lat","originlatitude","sourcelatitude","collectionlatitude"} for x in normalized)
    contains_lon=any(x in {"longitude","lon","long","originlongitude","sourcelongitude","collectionlongitude"} for x in normalized)
    contains_taxon=any(x in {"species","specie","taxon","speciesname","scientificname","speciesid","fullsp"} for x in normalized)
    contains_standardized_latitude = "standlat" in normalized
    contains_latitude_bin = "latbin" in normalized
    return {"sheet_dimension_from_source":dim,
            "first_row_header_cells":header or [],
            "contains_latitude_header":contains_lat,
            "contains_longitude_header":contains_lon,
            "contains_taxon_header":contains_taxon,
            "contains_standardized_latitude_header":contains_standardized_latitude,
            "contains_latitude_bin_header":contains_latitude_bin,
            "row_level_biological_response_values_read":False}


def file_schema(raw):
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("Original OSF file is not XLSX ZIP")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names=z.namelist()
        if "xl/workbook.xml" not in names:
            raise ValueError("Not an XLSX workbook")
        shared=shared_strings(z)
        sheets=[{"sheet_name":name,**read_only_first_row(z,path,shared)} for name,path in sheet_descriptors(z)]
    return sheets


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    contract=json.loads(CONTRACT.read_text())
    if contract.get("status")!="PRE_BINARY_SPREADSHEET_CONTENT_HEADER_ONLY_AUDIT":
        raise ValueError("Unfrozen original source header-audit plan")
    src=contract["source_file"]
    report={"schema":"ttf_weaving_2022_original_supplementary_data_1_header_only_v01",
            "source_file_id":src["id"],"source_file_name":src["name"],
            "source_sha256_expected":src["sha256"],"source_size_expected":src["size_bytes"],
            "biological_response_rows_inspected":False,
            "source_origin_coordinates_inspected":False,
            "real_genetic_or_plasticity_effect_test_performed":False}
    try:
        raw=get_public_osf_file(src["id"])
        report["source_download_bytes"]=len(raw)
        report["source_sha256_observed"]=hashlib.sha256(raw).hexdigest()
        if len(raw)!=src["size_bytes"] or report["source_sha256_observed"]!=src["sha256"]:
            raise ValueError("Source original OSF byte size/SHA256 differs")
        sheets=file_schema(raw)
        report["sheets"]=sheets
        report["sheet_count"]=len(sheets)
        report["origin_geolocation_status"]=("EXACT_LATITUDE_LONGITUDE_COLUMNS_PRESENT" if any(z["contains_latitude_header"] and z["contains_longitude_header"] for z in sheets) else "NO_EXACT_ORIGIN_LATITUDE_LONGITUDE_PAIR_SOURCE_ONLY")
        report["status"]="EXACT_SOURCE_HEADER_ONLY_AUDIT_PASS"
    except Exception as e:
        report["status"]="SOURCE_HEADER_AUDIT_STOP_NOT_BIOLOGICAL_NULL"
        report["error_type"]=type(e).__name__
        report["error_message"]=str(e)[:400]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps(report,sort_keys=True))


if __name__=="__main__":
    main()
