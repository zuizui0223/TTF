"""Official Dryad metadata-only verification for fresh transplant source cohorts.

Anonymous public GET of DOI->versions->file manifest, WITHOUT fetching any
source/fitness/phenotype/genotype content. Pinned selection is fixed before
response/model inspection. Errors are recorded as access gates, not absence.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import quote,urlparse
from urllib.request import Request,urlopen

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"source_only_original_dataset_contract_v01.json"
API="https://datadryad.org/api/v2/"
MAX_BODY=500000

def public_get(url):
    u=urlparse(url)
    if u.scheme!="https" or u.hostname!="datadryad.org" or not u.path.startswith("/api/v2/"):
        raise ValueError("Not an approved official Dryad HTTPS metadata API URL")
    req=Request(url,headers={"Accept":"application/json",
                             "User-Agent":"TTF-fresh-source-DOI-metadata-only/2026-11-19"})
    with urlopen(req,timeout=14) as resp:
        if urlparse(resp.geturl()).hostname!="datadryad.org":
            raise ValueError("Unexpected metadata host redirect")
        data=resp.read(MAX_BODY+1)
        if len(data)>MAX_BODY:raise ValueError("Metadata response exceeds size cap")
        return json.loads(data)

def flatten_api_records(payload,hint):
    if isinstance(payload,list):return payload
    if not isinstance(payload,dict):raise ValueError("Not a JSON API collection")
    embedded=payload.get("_embedded",{})
    if isinstance(embedded,dict):
        for key,value in embedded.items():
            if hint in key.lower() and isinstance(value,list):
                return value
        for value in embedded.values():
            if isinstance(value,list):return value
    for key in ("records","data",hint,"items"):
        if isinstance(payload.get(key),list):return payload[key]
    raise ValueError("Unknown Dryad version/file response schema")

def get_access(url,fetcher):
    try:
        return {"status":"PUBLIC_METADATA_AVAILABLE","body":fetcher(url)}
    except (HTTPError,URLError,TimeoutError,ValueError,OSError,json.JSONDecodeError) as e:
        return {"status":"PUBLIC_METADATA_ACCESS_NOT_CONFIRMED",
                "error_type":type(e).__name__,
                "http_code":getattr(e,"code",None)}

def check_study(study,fetcher=public_get):
    doi=study["doi"]
    if not doi.startswith("10.5061/dryad.") or "/" in doi:
        raise ValueError("DOI outside frozen Dryad scope")
    version_url=API+"datasets/"+quote("doi:"+doi,safe="")+"/versions?per_page=100"
    try:
        wanted=( [(study["file"],study["public_file_stream_id_from_Dryad_page"])]
                 if "public_file_stream_id_from_Dryad_page" in study else
                 [(study["source_location_file"],study["source_location_file_stream_id_from_Dryad_page"]),
                  (study["transplant_location_file"],study["transplant_location_file_stream_id_from_Dryad_page"])]
                 if "source_location_file_stream_id_from_Dryad_page" in study else
                 [(study["archive"],None)] )
    except KeyError as e:
        raise ValueError("Frozen candidate file identity incomplete") from e
    versions=get_access(version_url,fetcher)
    status="PUBLIC_VERSION_MANIFEST_NOT_CONFIRMED"
    results=[]
    version_id=None
    if versions["status"]=="PUBLIC_METADATA_AVAILABLE":
        try:
            records=flatten_api_records(versions["body"],"version")
            # Prefer exact latest available number (no response-dependent tuning).
            ranked=sorted((v for v in records if isinstance(v,dict) and
                           isinstance(v.get("id"),int)),key=lambda v:
                           (int(v.get("versionNumber",0)),int(v["id"])),reverse=True)
            if not ranked:raise ValueError("No public version ID")
            version_id=ranked[0]["id"]
            files_url=API+f"versions/{version_id}/files?per_page=100"
            files=get_access(files_url,fetcher)
            if files["status"]=="PUBLIC_METADATA_AVAILABLE":
                file_records=flatten_api_records(files["body"],"file")
                for path,pinned_id in wanted:
                    matches=[v for v in file_records if isinstance(v,dict) and
                             (v.get("path")==path or v.get("filename")==path)]
                    if len(matches)!=1:
                        results.append({"expected_original_path":path,
                                        "status":"ORIGINAL_FILE_NAME_NOT_UNIQUELY_REGISTERED"})
                        continue
                    m=matches[0]
                    actual_id=m.get("id")
                    match=pinned_id is None or actual_id==pinned_id
                    results.append({
                      "expected_original_path":path,
                      "public_file_id_from_landing_page":pinned_id,
                      "official_file_record_id":actual_id,
                      "official_file_size_bytes":m.get("size"),
                      "official_file_digest":m.get("digest"),
                      "official_file_digest_type":m.get("digestType"),
                      "file_id_matches_frozen_page":match,
                      "status":"ORIGINAL_FILE_METADATA_MATCHES_FROZEN_SOURCE" if match
                               else "FILE_IDENTIFIER_CHANGED__REQUALIFY_BEFORE_RAW_ACCESS",
                    })
                status=("ALL_PINNED_FILES_REGISTERED_IN_OFFICIAL_METADATA" if all(
                    x["status"]=="ORIGINAL_FILE_METADATA_MATCHES_FROZEN_SOURCE"
                    for x in results) else "ORIGINAL_FILE_METADATA_PARTIALLY_UNQUALIFIED")
            else:
                status="OFFICIAL_VERSION_FOUND_BUT_FILE_MANIFEST_ACCESS_FAILED"
                results.append({"file_list_http_error_code":files.get("http_code"),
                                "file_list_error_type":files.get("error_type")})
        except (ValueError,TypeError,KeyError) as e:
            status="OFFICIAL_VERSION_SCHEMA_NOT_VALIDATED"
            results.append({"schema_error_type":type(e).__name__})
    return {
      "study_id":study["study_id"],"dataset_doi":doi,
      "latest_public_version_id":version_id,
      "original_pinned_files":len(wanted),
      "files_metadata":results,"status":status,
      "version_metadata_request_status":versions["status"],
      "version_request_http_error_code":versions.get("http_code"),
      "version_request_error_type":versions.get("error_type"),
      "file_or_response_bytes_downloaded":False,
      "original_source_GPS_verified":False,
      "new_original_fitness_panel_qualified":False,
    }

def evaluate(contract,fetcher=public_get):
    if contract.get("schema")!="ttf_fresh_source_only_original_transplant_discovery_v01":
        raise ValueError("Frozen study selection is absent")
    studies=contract["candidates"]
    if [s["study_id"] for s in studies]!=[
        "Albano2025_Trifolium_repens","CrossEckert2026_Camissoniopsis_cheiranthifolia",
        "Carley2025_Clarkia_demographic"]:
        raise ValueError("Frozen candidate set/order changed")
    return {
      "schema":"ttf_fresh_official_Dryad_DOI_file_manifest_v01",
      "status":"PUBLIC_ORIGINAL_FILE_PROVENANCE_ONLY_NOT_SOURCE_GPS_OR_LGM_EFFECT",
      "studies":[check_study(s,fetcher=fetcher) for s in studies],
      "frozen_selection_before_biological_response":True,
      "source_GPS_independently_certified":0,
      "fitness_or_genetics_response_opened":False,
      "LGM21ka_contrast_qualified":False,
      "old_Bontrager_44of50_stop_unchanged":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    d=json.loads(CONTRACT.read_text("utf-8"))
    result=evaluate(d)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))

if __name__=="__main__":
    main()
