"""Official primary CCSM4 ~22ka LGM 2.5m BIO GeoTIFF ZIP index retry.

The original model role (CCSM4 primary) was frozen before all climate
values were opened. This is archive metadata only, NOT a fitness effect.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
from worldclim_lgm_remote_zip_directory_v01 import inspect_one,source_urls

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"lgm_paired_climate_source_only_contract_v01.json"

def evaluate(contract,inspector=inspect_one):
    if contract["models"][0]["role"]!="primary" or contract["models"][0]["code"]!="CCSM4":
        raise ValueError("Original precommitted primary GCM role changed")
    url=source_urls(contract)["CCSM4"]
    attempts=[];receipt=None
    for i in range(3):
        try:
            receipt=inspector(url)
            attempts.append({"attempt":i+1,"status":receipt["status"]})
            break
        except (TimeoutError,OSError,ValueError) as e:
            attempts.append({"attempt":i+1,"status":"OFFICIAL_ARCHIVE_ACCESS_ERROR",
                             "exception_type":type(e).__name__})
            if i<2:time.sleep(2)
    return {"schema":"ttf_Albano2025_CC_22ka_primary_original_ZIP_index_retry_v01",
            "status":"PRIMARY_CCSM4_ZIP_CENTRAL_INDEX_VERIFIED" if receipt else
                     "PRIMARY_CCSM4_ZIP_INDEX_NOT_YET_ACCESSIBLE",
            "source_url":url,"attempts":attempts,
            "zip_index":receipt,
            "original_CCSM4_raster_values_extracted":False,
            "fitness_or_genetic_responses_opened":False,
            "prior_bontrager_negative_unchanged":True}

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    r=evaluate(json.loads(CONTRACT.read_text("utf-8")))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(r,sort_keys=True))
if __name__=="__main__":main()
