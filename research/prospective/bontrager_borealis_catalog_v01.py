"""Public UBC Borealis mirror metadata only, not transplant fitness responses."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from urllib.request import Request,urlopen

SOURCE_DOI="doi:10.5683/SP2/DU3ZTQ"
URL="https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP2%2FDU3ZTQ"

def parse_catalog(body,expected_doi=SOURCE_DOI):
    if body.get("status")!="OK":raise ValueError("Official Borealis metadata not OK")
    d=body.get("data") or {}
    if d.get("persistentId")!=expected_doi:
        raise ValueError("DOI mismatch: do not replace source based on response data")
    version=d.get("latestVersion") or {}
    entries=version.get("files")
    if not isinstance(entries,list):raise ValueError("Missing official Borealis file metadata")
    files=[]
    for item in entries:
        f=item.get("dataFile") or {}
        files.append({"file_id":f.get("id"),"label":item.get("label"),
                      "file_name":f.get("filename"),"size_bytes":f.get("filesize"),
                      "md5":f.get("md5"),"restricted":item.get("restricted"),
                      "mime_type":f.get("contentType")})
    return {
       "schema":"ttf_Bontrager_UBC_Borealis_public_file_catalog_v01",
       "status":"BOREALIS_PUBLIC_SOURCE_METADATA_RETRIEVED",
       "persistentId":d.get("persistentId"),"dataset_id":d.get("id"),
       "version":str(version.get("versionNumber"))+"."+str(version.get("versionMinorNumber")),
       "version_release":version.get("releaseTime"),"file_count":len(files),"files":files,
       "fitness_records_read":False,"origin_destination_coordinates_read":False,
       "independent_from_Dryad_or_Hargreaves_confirmed":False,
       "no_adaptation_test_run":True
    }

def run():
    ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    report={"schema":"ttf_Bontrager_UBC_Borealis_public_file_catalog_v01",
            "fitness_records_read":False,"origin_destination_coordinates_read":False,
            "independent_from_Dryad_or_Hargreaves_confirmed":False,
            "no_adaptation_test_run":True}
    try:
        req=Request(URL,headers={"Accept":"application/json",
                                 "User-Agent":"TTF-original-transplant-alternate-repo-catalog-v01"})
        with urlopen(req,timeout=45) as response:raw=response.read(3000000)
        report.update(parse_catalog(json.loads(raw)))
    except Exception as e:
        report.update({"status":"BOREALIS_SOURCE_METADATA_UNAVAILABLE_NOT_BIOLOGICAL_NULL",
                       "error_type":type(e).__name__,"error_message":str(e)[:240]})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(report["status"],report.get("file_count"))
if __name__=="__main__":run()
