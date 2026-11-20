"""Response-blind WorldClim LGM ZIP public-source accessibility preflight.

Four preregistered official URLs only: WorldClim 1.4 paired baseline
and CCSM4, MIROC-ESM, MPI-ESM-P LGM (~22ka), all 2.5 arcmin.
Reads at most 512 bytes per endpoint; NEVER extracts ecological values.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from urllib.parse import urlparse
HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"lgm_paired_climate_source_only_contract_v01.json"
HOST="geodata.ucdavis.edu"
READ_MAX=512
EXPECTED={
   "WC1.4_MODERN":"https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/bio_2-5m_bil.zip",
   "CCSM4":"https://geodata.ucdavis.edu/climate/cmip5/lgm/cclgmbi_2-5m.zip",
   "MIROC-ESM":"https://geodata.ucdavis.edu/climate/cmip5/lgm/mrlgmbi_2-5m.zip",
   "MPI-ESM-P":"https://geodata.ucdavis.edu/climate/cmip5/lgm/melgmbi_2-5m.zip",
}

def classify_head(raw):
    if raw.startswith(b"PK\x03\x04"):return "ZIP_LOCAL_HEADER_OBSERVED_NOT_FULL_FILE_VERIFIED"
    if raw.lstrip().lower().startswith((b"<html",b"<!doctype html")):
        return "HTML_NOT_ZIP"
    if not raw:return "EMPTY_RESPONSE"
    return "UNRECOGNIZED_CONTENT"

def request_head(url):
    u=urlparse(url)
    if u.scheme!="https" or u.hostname!=HOST:
        raise ValueError("Official LGM archive URL not permitted")
    req=Request(url,headers={"Range":"bytes=0-511","Accept":"application/zip,*/*",
                             "User-Agent":"TTF-LGM-source-only-availability/2026-11-20"})
    with urlopen(req,timeout=18) as response:
        if urlparse(response.geturl()).hostname!=HOST:
            raise ValueError("Unexpected climate source host redirect")
        raw=response.read(READ_MAX)
        return {"http_status":response.status,
                "content_type":response.headers.get("Content-Type",""),
                "content_length_header":response.headers.get("Content-Length"),
                "content_range_header":response.headers.get("Content-Range"),
                "first_bytes_read":len(raw),
                "prefix_sha256":hashlib.sha256(raw).hexdigest(),
                "source_head_classification":classify_head(raw)}

def probe(contract,runner=request_head):
    if contract.get("schema")!="ttf_Albano2025_pre_response_LGM_climate_exposure_contract_v01":
        raise ValueError("Frozen LGM source contract missing or changed")
    urls={"WC1.4_MODERN":contract["baseline_url"]}
    for m in contract["models"]:
        if m["code"] in urls:raise ValueError("GCM duplicate")
        urls[m["code"]]=m["official_url"]
    if urls!=EXPECTED:raise ValueError("Official paired source URL or frozen GCM family drift")
    results={}
    for code,url in urls.items():
        try:
            res=runner(url)
            results[code]={"status":"OFFICIAL_ZIP_HEADER_PRESENT" if res["source_head_classification"].startswith("ZIP_LOCAL_HEADER") else "UNVERIFIED_OR_NONZIP_RESPONSE",**res}
        except (HTTPError,URLError,TimeoutError,ValueError,OSError) as e:
            results[code]={"status":"OFFICIAL_ARCHIVE_ACCESS_NOT_CONFIRMED","exception_type":type(e).__name__,
                           "http_error":getattr(e,"code",None)}
    return {
      "schema":"ttf_Albano2025_WorldClim_1_4_paired_LGM_source_transport_v01",
      "matched_baseline_period":"1960-1990",
      "LGM_snapshot_period":"about 22000 years BP, not exactly 21000 BP",
      "GCM_models":["CCSM4","MIROC-ESM","MPI-ESM-P"],
      "resolution_arcmin":2.5,
      "responses":results,
      "any_original_climate_zip_fully_downloaded":False,
      "any_LGM_or_current_bioclim_raster_values_extracted":False,
      "actual_source_GPS_records_submitted_to_CI":False,
      "biological_response_opened_or_tested":False,
      "status":"SOURCE_TRANSPORT_DIAGNOSTIC_ONLY_NOT_HISTORICAL_PREDICTOR_GAIN"
    }

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",required=True,type=Path)
    a=parser.parse_args()
    output=probe(json.loads(CONTRACT.read_text("utf-8")))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(output,sort_keys=True))
if __name__=="__main__":main()
