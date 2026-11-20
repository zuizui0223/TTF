"""First-shot official Dryad source metadata and CSV identity preflight.

No genotype, fitness, growth, herbivory, survival or LGM response values are
parsed or returned. Operates only on original file identity, CSV headers,
and (for Cross/Eckert's independent geographic-only files) location rows.
An original repository file is not automatically a source-population GPS.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, re
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlparse
from urllib.request import Request, urlopen

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"source_only_original_dataset_contract_v01.json"
ALLOWED_HOSTS=("datadryad.org","www.datadryad.org","s3.amazonaws.com","storage.googleapis.com")
SOURCE_COLUMNS={"Garden","Position","Population","Rep","ContOrigin","Latitude","Longitude"}
MAX_HEADER=8192
MAX_LOCATION_CSV=16384
MAX_API=300000

def trusted(url):
    parsed=urlparse(url)
    host=(parsed.hostname or "").lower()
    if parsed.scheme!="https" or not (
        host in ALLOWED_HOSTS or host.endswith(".s3.amazonaws.com")
        or host.endswith(".amazonaws.com")
    ):
        raise ValueError("URL did not resolve to an official Dryad or object-store HTTPS host")
    return host

def fetch(url,cap):
    trusted(url)
    req=Request(url,headers={"Accept":"text/csv,application/json,*/*",
                             "User-Agent":"TTF-independent-transport-source-only/2026-11-19",
                             "Range":f"bytes=0-{cap-1}"})
    with urlopen(req,timeout=16) as resp:
        trusted(resp.geturl())
        return resp.read(cap+1),resp.geturl(),resp.status,resp.headers.get("Content-Type","")

def attempt(url,limit):
    out={"endpoint_family":urlparse(url).path.split("/")[1] if len(urlparse(url).path.split("/"))>1 else "root"}
    try:
        content,final,status,mime=fetch(url,limit)
        if len(content)>limit:raise ValueError("Returned byte limit exceeded")
        out.update({"status":"BYTES_RECEIVED",
                    "http_status":status,
                    "final_host":trusted(final),
                    "content_type":mime,
                    "sample_bytes":len(content)})
        return out,content
    except (HTTPError,URLError,TimeoutError,ValueError,OSError) as e:
        out.update({"status":"ACCESS_NOT_CONFIRMED","error_type":type(e).__name__,
                    "http_status":getattr(e,"code",None)})
        return out,None

def dryad_file_urls(file_id):
    if not isinstance(file_id,int) or file_id<1:
        raise ValueError("File ID not a positive pinned integer")
    return (
        f"https://datadryad.org/api/v2/files/{file_id}/download",
        f"https://datadryad.org/downloads/file_stream/{file_id}",
    )

def inspect_csv_header(data):
    if data is None:return None
    sample=data.decode("utf-8-sig",errors="replace")
    if sample.lstrip().startswith(("<!DOCTYPE","<!doctype","<html","{","[")):
        raise ValueError("Returned body looks like HTML or JSON, not CSV")
    first=sample.splitlines()
    if not first or len(first[0])>MAX_HEADER-1:
        raise ValueError("CSV header missing or truncated")
    columns=next(csv.reader([first[0]]))
    if len(columns)<2 or len(set(columns))!=len(columns):
        raise ValueError("CSV header invalid or repeated")
    return columns

def source_geography_only_rows(data,headers):
    # This path is called ONLY on geographically separate CSV files whose
    # Dryad README says they contain locations, not on fitness/outcome CSVs.
    sample=data.decode("utf-8-sig",errors="strict")
    rows=list(csv.DictReader(io.StringIO(sample)))
    if len(rows)>100:raise ValueError("Geography-only file unexpectedly large")
    simplified={c.lower().strip().replace("_","").replace(" ",""):c for c in headers}
    lat=next((simplified[k] for k in ("lat","latitude") if k in simplified),None)
    lon=next((simplified[k] for k in ("lon","long","longitude") if k in simplified),None)
    if not lat or not lon:return {"location_row_count":len(rows),
                                 "coordinates_parsed":False,
                                 "reason":"No clearly named lat/lon fields"}
    positions=[]
    for row in rows:
        latitude,longitude=float(row[lat]),float(row[lon])
        if not (-90<=latitude<=90 and -180<=longitude<=180):
            raise ValueError("Geography-only location coordinate out of range")
        positions.append((latitude,longitude))
    return {"location_row_count":len(rows),
            "unique_location_pairs":len(set(positions)),
            "coordinates_parsed":True,
            "independent_original_source_GPS_certified":False}

def preflight_file(file_id,name,geography_only=False):
    attempts=[]
    for url in dryad_file_urls(file_id):
        attempt_metadata,data=attempt(url,MAX_LOCATION_CSV if geography_only else MAX_HEADER)
        attempts.append(attempt_metadata)
        if data is None:continue
        try:
            header=inspect_csv_header(data)
            if geography_only:
                geo=source_geography_only_rows(data,header)
                if not geo["location_row_count"]:
                    raise ValueError("No locations in source-only file")
            else:
                missing=SOURCE_COLUMNS-set(header)
                if missing:raise ValueError("Required original population/garden/coordinate headers missing")
                geo=None
            return {
                "file":name,"public_Dryad_file_id":file_id,
                "status":"GEOGRAPHIC_ONLY_ORIGINAL_FILE_READ" if geography_only
                         else "ORIGINAL_EXPERIMENT_FILE_HEADER_CONFIRMED_NO_OUTCOME_ROWS_PARSED",
                "header_columns":header,
                "official_access_attempts":attempts,
                "download_sample_sha256":hashlib.sha256(data).hexdigest(),
                "download_sample_bytes":len(data),
                "full_file_cryptographic_integrity_verified":False,
                "location_file_summary":geo,
                "response_values_parsed_or_exported":False,
                "original_source_GPS_independently_certified":False,
            }
        except (UnicodeError,ValueError,csv.Error) as e:
            attempts[-1]["source_only_validation_error"]=str(e)
    return {"file":name,"public_Dryad_file_id":file_id,
            "status":"ORIGINAL_BYTES_OR_SOURCE_SCHEMA_NOT_CONFIRMED",
            "official_access_attempts":attempts,
            "original_source_GPS_independently_certified":False,
            "response_values_parsed_or_exported":False}

def evaluate(contract,runner=preflight_file):
    if contract.get("schema")!="ttf_fresh_source_only_original_transplant_discovery_v01":
        raise ValueError("Frozen discovery scope changed")
    candidates=contract["candidates"]
    if [r["study_id"] for r in candidates] != [
      "Albano2025_Trifolium_repens",
      "CrossEckert2026_Camissoniopsis_cheiranthifolia",
      "Carley2025_Clarkia_demographic",
    ]:raise ValueError("Candidate order or eligibility set changed")
    if (contract["certified_independent_original_source_GPS_populations"]!=0
        or contract["certified_original_experiment_cohorts"]!=0
        or contract["LGM_fitness_or_genetic_response_opened"]
        or not contract["old_bontrager_44_of_50_terminal_stop_unchanged"]):
        raise ValueError("Previously frozen science/source-only precondition changed")
    albano=candidates[0];cross=candidates[1]
    files=[
       {"study_id":albano["study_id"],
        "result":runner(albano["public_file_stream_id_from_Dryad_page"],albano["file"],False)},
       {"study_id":cross["study_id"],
        "result":runner(cross["source_location_file_stream_id_from_Dryad_page"],cross["source_location_file"],True)},
       {"study_id":cross["study_id"],
        "result":runner(cross["transplant_location_file_stream_id_from_Dryad_page"],cross["transplant_location_file"],True)},
    ]
    return {
       "schema":"ttf_fresh_Dryad_primary_source_access_preflight_v01",
       "status":"SOURCE_ONLY_TRANSPORT_CHECK_NOT_LGM_PANEL_CERTIFICATION",
       "study_candidates_before_response":len(candidates),
       "official_files_probed":len(files),
       "files":files,
       "Carley_archive_7MB_deliberately_not_downloaded":True,
       "biological_response_rows_parsed_or_modeled":False,
       "new_independent_source_population_GPS_certified":0,
       "new_independent_experiment_cohorts_certified":0,
       "LGM21ka_contrast_qualified":False,
       "prior_bontrager_44of50_stop_unchanged":True,
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    d=json.loads(CONTRACT.read_text("utf-8"))
    result=evaluate(d)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))

if __name__=="__main__":
    main()
