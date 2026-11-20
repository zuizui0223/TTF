"""Sample precommitted primary CCSM4 LGM BIO1/BIO12 on 95 unlabeled Albano cells.

Uses the exact MPI-tested ZIP Range, deflate, CRC32 and GeoTIFF workflow.
Only the GCM source and official central-index member identities differ.
No original plant ID/latitude/longitude, phenotype or fitness rows loaded.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import albano2025_MPI_LGM_raster_values_unlabeled_v01 as engine

HERE=Path(__file__).resolve().parent
MANIFEST=HERE/"source_registry/albano2025_CCSM4_primary_LGM_tiff_members_v01.json"
CONTRACT=HERE/"lgm_paired_climate_source_only_contract_v01.json"

def evaluate(manifest,contract,grid=None,engine_module=engine):
    if (manifest.get("schema")!="ttf_Albano2025_CCSM4_primary_original_ZIP_two_GeoTIFF_member_index_v01"
        or manifest["original_archive_total_bytes"]!=219610901
        or manifest["original_archive_member_count"]!=19):
        raise ValueError("Frozen primary CCSM4 source index identity drift")
    if (contract["models"][0]["code"],contract["models"][0]["role"])!=("CCSM4","primary"):
        raise ValueError("Original precommitted CCSM4 GCM role changed")
    url=manifest["original_archive_url"]
    if url!=contract["models"][0]["official_url"] or not url.endswith("/cclgmbi_2-5m.zip"):
        raise ValueError("Original primary CCSM4 official source URL changed")
    source={x["filename"]:x for x in manifest["members"]}
    if set(source)!={"cclgmbi1.tif","cclgmbi12.tif"}:
        raise ValueError("Wrong original BIO variables in primary GCM member family")
    entries={}
    for bio,name in [("BIO1","cclgmbi1.tif"),("BIO12","cclgmbi12.tif")]:
        z=source[name]
        if z["method"]!=8:raise ValueError("Primary ZIP member must be deflated")
        entries[bio]={
          "name":z["filename"],"local_header_offset":z["local_header_offset"],
          "compressed_bytes":z["compressed_bytes"],
          "uncompressed_bytes":z["uncompressed_bytes"],"crc32":z["CRC32_hex"]}
    if grid is None:grid=json.loads(engine_module.CELLS.read_text("utf-8"))
    old=(engine_module.URL,engine_module.ARCHIVE_BYTES,engine_module.ENTRIES)
    try:
        engine_module.URL=url
        engine_module.ARCHIVE_BYTES=manifest["original_archive_total_bytes"]
        engine_module.ENTRIES=entries
        out=engine_module.evaluate(grid)
    finally:
        engine_module.URL,engine_module.ARCHIVE_BYTES,engine_module.ENTRIES=old
    out.update({
       "schema":"ttf_Albano2025_CCSM4_PRIMARY_source_only_2_5m_LGM_cells_v01",
       "status":"FROZEN_PRIMARY_CCSM4_LGM_RASTER_SAMPLED_NOT_BIOLOGICAL_EFFECT",
       "GCM":"CCSM4","GCM_role":"primary_precommitted",
       "source_manifest_reference":"source_registry/albano2025_CCSM4_primary_LGM_tiff_members_v01.json",
    })
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    result=evaluate(json.loads(MANIFEST.read_text("utf-8")),
                    json.loads(CONTRACT.read_text("utf-8")))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
