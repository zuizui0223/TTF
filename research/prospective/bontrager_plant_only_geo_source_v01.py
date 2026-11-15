"""New plant-only transplant source metadata stratum; no fitness values.

This does not replace the frozen all-taxon source census. Restricts only
official source functional_group labels before looking at observed fitness
or other response columns. Original Dryad/UBC source ZIP identity required.
"""
from __future__ import annotations
import argparse,io,json,zipfile
from pathlib import Path
from urllib.request import urlopen,Request
from bontrager_mirror_zip_v01 import verify_original_zip,ORIGINAL_URL
from bontrager_transplant_geo_census_v01 import metadata_census

HERE=Path(__file__).resolve().parent
AMENDMENT=HERE/"bontrager_plant_functional_group_amendment_v01.json"
SOURCE=HERE/"bontrager_mirror_zip_contract_v01.json"
GEO_CONTRACT=HERE/"bontrager_transplant_geo_contract_v01.json"


def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    amended=json.loads(AMENDMENT.read_text())
    base=json.loads(GEO_CONTRACT.read_text())
    source=json.loads(SOURCE.read_text())
    if amended.get("status")!="FROZEN_AFTER_SOURCE_LABEL_DISTRIBUTION_BEFORE_ANY_FITNESS_RESPONSE":
        raise ValueError("Plant-only original taxonomy qualification not frozen")
    if amended["source_file_sha256"]!=source["dryad_sha256"]:
        raise ValueError("Plant-only original publisher archive differs from first all-taxon source")
    out={"schema":"ttf_Bontrager_plant_only_original_transplant_source_geography_v01",
         "fitness_measurement_sign_or_magnitude_opened":False,
         "adaptive_fitness_or_speciation_prediction_performed":False,
         "old_original_frozen_Haldane_status_modified":False,
         "new_fitness_response_opening_authorized":False,
         "original_unselected_all_taxa_source_receipt_preserved":True}
    try:
        req=Request(ORIGINAL_URL,headers={"User-Agent":"TTF-plant-only-Bontrager-source-qualification-v01"})
        with urlopen(req,timeout=90) as resp:raw=resp.read(source["max_download_bytes"]+1)
        verify_original_zip(raw,source)
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            with z.open(base["source_authority"]["source_member"]) as file:
                rows=file.read(25_000_000)
        if len(rows)>=25_000_000:raise ValueError("Source CSV exceeds source-only cap")
        groups=set(amended["plant_original_functional_group_categories_exact"])
        out.update(metadata_census(rows,base,allowed_functional_groups=groups))
        if out["source_metadata_exclusion_counts"].get("excluded_nonplant_original_functional_group")!=92:
            raise ValueError("Original 92 nonplant rows changed; source-only STOP")
        if sum(out["functional_group_metadata_count"].values())!=7452:
            raise ValueError("Exactly 7452 source-functional-group plant records expected")
        out["source_population_stratum"]="ORIGINAL_FUNCTIONAL_GROUP_VASCULAR_PLANTS_ONLY"
        out["plant_functional_group_labels"]=sorted(groups)
        out["plant_stratum_qualification_status"]=out["status"]
        out["fitness_measurement_sign_or_magnitude_opened"]=False
        out["new_fitness_response_opening_authorized"]=False
        out["source_zip_sha256"]=source["dryad_sha256"]
    except Exception as exc:
        out.update(status="PLANT_SOURCE_METADATA_STOP_NOT_ADAPTATION_NULL",
                   error_type=type(exc).__name__,error_message=str(exc)[:330])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("PLANT_ONLY_TRANSPLANT_CELLS",out["status"],
          out.get("qualified_binomial_taxa"),out.get("qualified_original_source_studies"),
          out.get("shared_garden_fitness_type_cells_with_multiple_source_IDs_and_coordinates"))
if __name__=="__main__":main()
