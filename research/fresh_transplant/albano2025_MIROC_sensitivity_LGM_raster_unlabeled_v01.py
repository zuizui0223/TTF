"""Response-blind MIROC-ESM LGM sensitivity model original archive and 95-cell samples.

GCM role frozen as sensitivity before observing any of three model climate
fields. At run time retrieves official ZIP index (retry), authenticates
19-member archive and CRC-pins original BIO1/BIO12 members before
retrieving only those two GeoTIFFs. No plant outcomes or original GPS.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
from worldclim_lgm_remote_zip_directory_v01 import inspect_one,source_urls
import albano2025_MPI_LGM_raster_values_unlabeled_v01 as engine

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"lgm_paired_climate_source_only_contract_v01.json"
TARGET="MIROC-ESM"
EXPECTED="https://geodata.ucdavis.edu/climate/cmip5/lgm/mrlgmbi_2-5m.zip"
NAMES={"BIO1":"mrlgmbi1.tif","BIO12":"mrlgmbi12.tif"}

def evaluate(contract,inspector=inspect_one,climate_engine=engine):
    if (contract["models"][1]["code"]!=TARGET or
        contract["models"][1]["role"]!="sensitivity" or
        contract["models"][1]["official_url"]!=EXPECTED):
        raise ValueError("Original frozen MIROC sensitivity GCM role/source changed")
    url=source_urls(contract)[TARGET]
    events=[];index=None
    for attempt in range(3):
        try:
            index=inspector(url)
            events.append({"attempt":attempt+1,"status":"OFFICIAL_ZIP_INDEX_ACQUIRED"})
            break
        except (ValueError,OSError,TimeoutError) as err:
            events.append({"attempt":attempt+1,"status":"NETWORK_OR_FORMAT_UNVERIFIED",
                           "error_type":type(err).__name__})
            if attempt<2:time.sleep(2)
    if index is None:
        return {"schema":"ttf_Albano2025_MIROC_sensitivity_source_only_LGM_v01",
                "status":"OFFICIAL_MIROC_ARCHIVE_SOURCE_ACCESS_HOLD",
                "source_events":events,
                "biological_response_read":False,"LGM_values_extracted":False}
    if index["archive_file_count"]!=19 or index["official_archive_total_bytes"]<200000000:
        raise ValueError("MIROC official 2.5m 19-member ZIP index changed")
    members={}
    for code,name in NAMES.items():
        found=[x for x in index["files"] if x["filename"]==name]
        if len(found)!=1 or found[0]["method"]!=8:
            raise ValueError("Official original MIROC BIO1/BIO12 ZIP members unqualified")
        z=found[0]
        members[code]={"name":name,"local_header_offset":z["local_header_offset"],
                       "compressed_bytes":z["compressed_bytes"],
                       "uncompressed_bytes":z["uncompressed_bytes"],
                       "crc32":z["CRC32_hex"]}
    original=(climate_engine.URL,climate_engine.ARCHIVE_BYTES,climate_engine.ENTRIES)
    try:
        climate_engine.URL=url
        climate_engine.ARCHIVE_BYTES=index["official_archive_total_bytes"]
        climate_engine.ENTRIES=members
        grid=json.loads(climate_engine.CELLS.read_text("utf-8"))
        sampled=climate_engine.evaluate(grid)
    finally:
        climate_engine.URL,climate_engine.ARCHIVE_BYTES,climate_engine.ENTRIES=original
    sampled.update({
      "schema":"ttf_Albano2025_MIROC_ESM_sensitivity_LGM_2_5m_cells_v01",
      "status":"FROZEN_MIROC_SENSITIVITY_LGM_RASTER_SAMPLED_NOT_BIOLOGICAL_EFFECT",
      "GCM":TARGET,
      "GCM_role":"sensitivity_precommitted",
      "source_events":events,
      "official_original_zip_archive_bytes":index["official_archive_total_bytes"],
      "original_PINNED_zip_member_metadata":{k:{x:v[x] for x in (
          "local_header_offset","compressed_bytes","uncompressed_bytes","crc32","name")}
          for k,v in members.items()},
    })
    return sampled

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();d=json.loads(CONTRACT.read_text("utf-8"))
    x=evaluate(d)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(x,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(x,sort_keys=True))
if __name__=="__main__":main()
