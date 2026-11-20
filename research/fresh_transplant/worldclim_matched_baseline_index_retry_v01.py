"""Retry official WorldClim 1.4 current 1960–1990 BIL ZIP central directory.

Source-only, one archive, no raster/fitness values, fail-closed attempts.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
from worldclim_lgm_remote_zip_directory_v01 import inspect_one,source_urls
HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"lgm_paired_climate_source_only_contract_v01.json"

def evaluate(contract,worker=inspect_one):
    source=source_urls(contract)["WC1.4_MODERN"]
    attempts=[]
    result=None
    for i in range(3):
        try:
            result=worker(source)
            attempts.append({"attempt":i+1,"status":result["status"]})
            break
        except (OSError,ValueError,TimeoutError) as exc:
            attempts.append({"attempt":i+1,"status":"ACCESS_ERROR",
                             "error_type":type(exc).__name__})
            if i<2:time.sleep(2)
    return {
      "schema":"ttf_Albano2025_WC1_4_1960_1990_bio_BIL_zip_source_index_v01",
      "source_url":source,
      "status":"OFFICIAL_CURRENT_BASELINE_CENTRAL_INDEX_VERIFIED" if result else
               "CURRENT_BASELINE_INDEX_ACCESS_HOLD",
      "attempts":attempts,"index":result,
      "LGM_or_modern_raster_values_read":False,
      "fitness_or_genetic_response_read":False,
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    out=evaluate(json.loads(CONTRACT.read_text("utf-8")))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(out,sort_keys=True))
if __name__=="__main__":main()
