"""Public PROCEED v6.1 source file-catalog inventory: metadata ONLY.

Never fetch any phenotypic rate table or effect-size row. This is source
discovery for a *different possible* adaptation study, not a reopened TTF
genetic response or confirmatory test.
"""
from __future__ import annotations

import argparse
import json
import urllib.request
import urllib.parse
from pathlib import Path

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed_v61_source_metadata_contract_v01.json"


def api_metadata_get(url,maximum):
    parsed=urllib.parse.urlparse(url)
    if parsed.scheme!="https" or parsed.hostname!="borealisdata.ca":
        raise ValueError("Read only official Borealis dataset-metadata host")
    if "/api/datasets/" not in parsed.path:
        raise ValueError("Only dataset-level Dataverse metadata requests")
    request=urllib.request.Request(
        url,headers={"Accept":"application/json","User-Agent":"TTF-PROCEED-public-dataset-metadata-only-v01"})
    with urllib.request.urlopen(request,timeout=40) as r:
        raw=r.read(maximum+1)
    if len(raw)>maximum:raise ValueError("Metadata exceeded source-only maximum")
    return json.loads(raw)


def source_catalog(payload):
    if payload.get("status")!="OK":
        raise ValueError("Borealis dataset API source metadata unavailable")
    data=payload.get("data") or {}
    version=data.get("latestVersion") or {}
    files=version.get("files") or []
    entries=[]
    for raw in files:
        f=raw.get("dataFile") or {}
        entries.append({
            "filename":raw.get("label") or f.get("filename") or "",
            "file_id":f.get("id"),
            "mime_type":f.get("contentType"),
            "size_bytes":f.get("filesize"),
            "md5_metadata":f.get("md5"),
            "restricted":raw.get("restricted"),
            "version":f.get("version"),
            "description":(raw.get("description") or "")[:400],
        })
    entries.sort(key=lambda x:(x["filename"],str(x["file_id"])))
    blocks=(version.get("metadataBlocks") or {})
    citation=blocks.get("citation") or {}
    fields=citation.get("fields") or []
    simple_fields={}
    for r in fields:
        if r.get("typeName") in ("title","subtitle","publicationDate","dateOfDeposit","datasetContact","author"):
            if r.get("typeName") in ("title","publicationDate","dateOfDeposit","subtitle"):
                simple_fields[r["typeName"]]=r.get("value")
    return {"schema":"ttf_proceed_v61_borealis_dataset_public_metadata_v01",
            "status":"DATAVERSE_METADATA_LIST_SUCCESS_RESPONSE_NOT_OPENED",
            "dataset_id":data.get("id"),
            "persistent_id":data.get("persistentId"),
            "release_number":str(version.get("versionNumber"))+"."+str(version.get("versionMinorNumber")),
            "version_state":version.get("versionState"),
            "source_date":version.get("releaseTime") or version.get("createTime"),
            "file_count":len(entries),
            "files":entries,
            "citation_metadata":simple_fields,
            "phenotypic_effect_size_rows_opened":False,
            "genetic_study_results_opened":False,
            "historical_climate_values_sampled":False,
            "history_adaptation_statistical_effect_fitted":False,
            "not_a_climate_or_speciation_result":True}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    rules=json.loads(CONTRACT.read_text())
    if rules.get("status")!="FROZEN_PRE_PHENOTYPIC_RESPONSE_READ_DATAVERSE_METADATA_ONLY":
        raise ValueError("Source discovery authorization absent")
    result={"schema":"ttf_proceed_v61_borealis_dataset_public_metadata_v01",
            "phenotypic_effect_size_rows_opened":False,
            "genetic_study_results_opened":False,
            "historical_climate_values_sampled":False,
            "history_adaptation_statistical_effect_fitted":False,
            "not_a_climate_or_speciation_result":True}
    try:
        p=api_metadata_get(rules["source"]["dataset_api"],rules["limits"]["max_bytes"])
        result=source_catalog(p)
        if result["persistent_id"]!="doi:10.5683/SP3/NXSL3Q":
            raise ValueError("Returned dataset persistent DOI does not match source")
    except Exception as e:
        result["status"]="SOURCE_METADATA_UNAVAILABLE_NOT_ECOLOGICAL_NULL"
        result["error_type"]=type(e).__name__
        result["error_message"]=str(e)[:350]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
