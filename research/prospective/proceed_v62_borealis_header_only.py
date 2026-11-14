"""PROCEED v6.2 official datafile header-only schema audit.

Uses the public Dataverse dataset catalog to verify exact registered file
ID, name, size and MD5, then streams JUST one header line per selected
source tabular file. Does not read effect-size rows, download a full file,
extract climate, fit a model, inspect original genotypes or adapt a frozen TTF
test after an outcome.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed_v62_source_header_contract_v01.json"
MAX_HEADER=131072


def validated_get_catalog(url):
    p=urlparse(url)
    if p.scheme!="https" or p.netloc!="borealisdata.ca" or not p.path.startswith("/api/datasets/"):
        raise ValueError("Only official Dataverse dataset metadata is allowed")
    req=Request(url,headers={"Accept":"application/json",
                             "User-Agent":"TTF-PROCEED-v62-headers-only-v01"})
    with urlopen(req,timeout=40) as resp:
        b=resp.read(3_000_001)
    if len(b)>3_000_000:
        raise ValueError("Catalog metadata too large")
    return json.loads(b)


def verify_frozen_catalog(payload,expected):
    if payload.get("status")!="OK":
        raise ValueError("Borealis source metadata not OK")
    d=payload.get("data") or {}
    if d.get("id")!=expected["borealis_dataset_id"]:
        raise ValueError("Dataset ID drift")
    ver=d.get("latestVersion") or {}
    label=f"{ver.get('versionNumber')}.{ver.get('versionMinorNumber')}"
    if label!=expected["repository_dataset_version"]:
        raise ValueError("Borealis dataset version changed after freeze")
    files={f.get("dataFile",{}).get("id"):f for f in ver.get("files",[])}
    reports=[]
    for expected_file in expected["files"]:
        f=files.get(expected_file["id"]) or {}
        df=f.get("dataFile") or {}
        if (df.get("filename")!=expected_file["name"] or
            df.get("filesize")!=expected_file["bytes"] or
            df.get("md5")!=expected_file["md5"] or
            df.get("contentType")!="text/tab-separated-values"):
            raise ValueError("Source file ID/name/type/size/md5 changed: "+expected_file["name"])
        if f.get("restricted"):
            raise ValueError("No access to restricted data in this screen")
        reports.append({"file_id":expected_file["id"],"filename":expected_file["name"],
                        "size_bytes":expected_file["bytes"],"source_md5":expected_file["md5"]})
    return reports


def stream_one_header(file_id, maxsize=MAX_HEADER):
    if not isinstance(file_id,int) or file_id<=0:
        raise ValueError("Invalid frozen source file ID")
    url=f"https://borealisdata.ca/api/access/datafile/{file_id}"
    request=Request(url,headers={"Accept":"text/tab-separated-values",
                                "User-Agent":"TTF-PROCEED-v62-header-only-v01"})
    with urlopen(request,timeout=60) as response:
        typ=response.headers.get("Content-Type","")
        # Dataverse ingested tabular files may stream as text/tab-separated-values
        line=response.readline(maxsize+1)
    if len(line)>maxsize:
        raise ValueError("TSV header exceeds header-only maximum")
    header=line.decode("utf-8-sig").strip("\r\n")
    if "\t" not in header:
        raise ValueError("Datafile header missing TSV delimiters; likely API error/redirect")
    names=header.split("\t")
    if any(not a.strip() for a in names):
        raise ValueError("Blank original source data column name")
    if len(set(names))!=len(names):
        raise ValueError("Duplicate column names in original source")
    return {"file_id":file_id,"source_content_type":typ,
            "header_column_count":len(names),"header_columns":names,
            "rows_two_plus_read":False}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    c=json.loads(CONTRACT.read_text())
    if c["status"]!="FROZEN_BOREALIS_PROCEED_V62_NEW_SOURCE_HEADER_ONLY":
        raise ValueError("No authorized source-only v6.2 header plan")
    e=c["current_source"]
    output={"schema":"ttf_proceed_v62_borealis_verified_source_headers_v01",
            "dataset_doi":e["dataset_doi"],"dataset_release":e["repository_dataset_version"],
            "rate_or_trait_rows_read":False,"genetic_values_read":False,
            "LGM_climate_sampled":False,"adaptive_evolution_test_run":False}
    try:
        catalog=validated_get_catalog("https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        checked=verify_frozen_catalog(catalog,e)
        headers=[stream_one_header(item["file_id"]) for item in checked]
        for h,item in zip(headers,checked):
            h.update(item)
        output["source_files"]=headers
        output["source_file_count"]=len(headers)
        output["status"]="SOURCE_HEADERS_VERIFIED_NO_RATE_VALUES_READ"
    except Exception as e:
        output["status"]="SOURCE_SCHEMA_ABSTAIN_NO_BIOLOGICAL_RESULT"
        output["error_type"]=type(e).__name__
        output["error_message"]=str(e)[:400]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    print(json.dumps(output,sort_keys=True))


if __name__=="__main__":
    main()
