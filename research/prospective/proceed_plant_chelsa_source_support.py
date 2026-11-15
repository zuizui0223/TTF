"""Response-blind 70-plant 21ka/current CHELSA source support; no trait rate fitted.

Reproduces the original study–species disjoint split from the source-only
census. Samples four exactly SHA-pinned CHELSA-TraCE21k bioclim layers at
21ka and present on the ORIGINAL PROCEED collection localities. No GBIF
centroids, imputation, genetic-effect data or trait-change rates are used.
All origin coordinates and pointwise climate values stay ephemeral in CI.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,urllib.request
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
from proceed_v62_metadata_row_census import (
    fetch_original_datafile_bytes,valid_binomial,valid_coordinates)
from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT,NATIVE_LOCK,ROW_CONTRACT
from proceed_v62_plant_genetic_source_census import (
    component_member_sets,split_connected_components)
from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
from site_feature_core import raster_values

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed_plant_chelsa_source_support_contract_v01.json"
NAMES=["bio01_-190","bio07_-190","bio12_-190","bio15_-190",
       "bio01_20","bio07_20","bio12_20","bio15_20"]


def component_assignment(study_species):
    comps=component_member_sets(study_species)
    def key(comp):
        raw="|".join(sorted(comp[0]))
        h=hashlib.sha256(("ttf-proceed-v62-plant-component-v01|"+raw).encode()).hexdigest()
        return (-len(comp[0]),h)
    comps.sort(key=key)
    mapping={}
    members={"train":set(),"heldout":set()}
    for spp,studies in comps:
        arm="train" if len(members["train"])<=len(members["heldout"]) else "heldout"
        for sp in spp:
            if sp in mapping:raise ValueError("Duplicate species role")
            mapping[sp]=arm
        members[arm].update(spp)
    digest=hashlib.sha256(json.dumps({k:sorted(v) for k,v in members.items()},
                                  sort_keys=True,separators=(",",":")).encode()).hexdigest()
    summary=split_connected_components(study_species)
    if digest!=summary["source_only_component_membership_sha256"]:
        raise ValueError("Component role rule inconsistent with original")
    return mapping,summary,digest


def source_geometry(raw,contract):
    allow=json.loads(ROW_CONTRACT.read_text())["source_metadata_columns_authorized"]
    prohibited=set(json.loads(ROW_CONTRACT.read_text())["outcome_fields_must_not_be_decoded"])
    if set(allow)&prohibited:raise ValueError("Phenotype response field in source metadata schema")
    reader=csv.reader(io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",newline=""),
                      delimiter=",",strict=True)
    header=next(reader)
    if len(header)!=69 or len(set(header))!=69 or not set(allow).issubset(header):
        raise ValueError("Original PROCEED native source header changed")
    ix={k:header.index(k) for k in allow}
    rows=[];study_species=defaultdict(set)
    for values in reader:
        if not values:continue
        if len(values)!=69:raise ValueError("Original source record width changed")
        # Do not access outcome columns in any row.
        m={k:values[pos].strip() for k,pos in ix.items()}
        if (m["phylum"]!="Streptophyta" or m["genphen"]!="Genetic" or m["taxa"]!="Plant"):
            continue
        a=valid_coordinates(m["sample1_latitude"],m["sample1_longitude"])
        b=valid_coordinates(m["sample2_latitude"],m["sample2_longitude"])
        sp=valid_binomial(m["sp_pub"]) or valid_binomial(m["sp_ncbi"])
        ref=m["ref_ID"]
        if not (a and b and sp and ref):continue
        study_species[ref].add(sp)
        # Coordinate representation is lon/lat for rasterio sampling.
        rows.append((sp,ref,(a[1],a[0]),(b[1],b[0])))
    return rows,study_species


def verified_raster(a,path):
    url=a["url"]
    parsed=urlparse(url)
    if parsed.scheme!="https" or parsed.hostname!="os.zhdk.cloud.switch.ch":
        raise ValueError("Source raster URL not original pinned CHELSA host")
    expected_bytes=int(a["size_bytes"]);expected_sha=a["sha256"]
    if path.exists() and path.stat().st_size==expected_bytes:
        h=hashlib.sha256()
        with path.open("rb") as f:
            for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
        if h.hexdigest()==expected_sha:return
    request=urllib.request.Request(url,headers={"User-Agent":"TTF-PROCEED70-response-blind-CHELSA-source-v01"})
    h=hashlib.sha256();n=0
    with urllib.request.urlopen(request,timeout=180) as r,path.open("wb") as f:
        while True:
            b=r.read(1024*1024)
            if not b:break
            n+=len(b)
            if n>expected_bytes:
                path.unlink(missing_ok=True)
                raise ValueError("CHELSA source exceeded SHA-pinned exact byte size")
            f.write(b);h.update(b)
    if n!=expected_bytes or h.hexdigest()!=expected_sha:
        path.unlink(missing_ok=True)
        raise ValueError("CHELSA original exact byte and SHA256 mismatch")


def source_climate_coverage(cohort,studies,climates,output):
    roles,split,digest=component_assignment(studies)
    if len(cohort)!=1561 or len(roles)!=70 or split["species_per_role"]!={"train":35,"heldout":35}:
        raise ValueError("Original 70-plant response-blind panel geometry differs")
    coords=sorted({xy for _,_,first,second in cohort for xy in (first,second)})
    indices={site:i for i,site in enumerate(coords)}
    xy=np.array(coords,float)
    if xy.shape[1]!=2 or not np.isfinite(xy).all():
        raise ValueError("Original georeferenced coordinates invalid")
    x=np.column_stack(climates)
    if x.shape!=(len(xy),8):
        raise ValueError("Original climate source matrix not 8 columns")
    complete=np.isfinite(x).all(axis=1)
    species_full=set();studies_full=set();pairs_complete=0;genes=set();roles_full={"train":set(),"heldout":set()}
    pairs_climate_full=[]
    for sp,study,c1,c2 in cohort:
        i=indices[c1];j=indices[c2]
        if complete[i] and complete[j]:
            species_full.add(sp);studies_full.add(study);genes.add(sp.split()[0])
            roles_full[roles[sp]].add(sp)
            pairs_complete+=1;pairs_climate_full.append((i,j))
    return {
        "status":None,
        "original_source_qualified_species":len(roles),
        "original_source_qualified_studies":len(studies),
        "original_source_study_species_split_sha256":digest,
        "original_geo_pairs_with_valid_species_study":len(cohort),
        "unique_original_georeferenced_lonlat_sites":len(coords),
        "all_eight_CHELSA_site_complete":int(complete.sum()),
        "all_eight_CHELSA_site_incomplete":int((~complete).sum()),
        "all_eight_climate_complete_pair_records":pairs_complete,
        "climate_complete_species":len(species_full),
        "climate_complete_source_studies":len(studies_full),
        "climate_complete_genera":len(genes),
        "climate_complete_species_per_original_role":{r:len(v) for r,v in roles_full.items()},
        "original_lonlat_byte_sha256":hashlib.sha256(xy.astype("<f8").tobytes()).hexdigest(),
        "climate_raw_numeric_matrix_sha256":hashlib.sha256(x.astype("<f8").tobytes()).hexdigest(),
        "raw_climate_or_coordinates_publicly_exported":False,
        "phenotype_rate_values_opened":False
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--work",required=True,type=Path)
    args=parser.parse_args()
    contract=json.loads(CONTRACT.read_text())
    if contract.get("status")!="FROZEN_BEFORE_SOURCE_ONLY_ORIGINAL_CHELSA_CLIMATE_VALUES":
        raise ValueError("Source-only climate sampling was not frozen")
    src=json.loads(SOURCE_CONTRACT.read_text())["current_source"]
    native=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    asset_path=HERE/"chelsa_trace21k_fixed8_assets_v01.json"
    asset=json.loads(asset_path.read_text())
    if len(asset["assets"])!=8 or asset["schema"]!="ttf_jensen80_original_C1_CHELSA8_v01":
        raise ValueError("Focal historical climate input family changed")
    labels=[f"{a['variable']}_{a['time_index']}" for a in asset["assets"]]
    if labels!=NAMES:raise ValueError("Original historical climate predictor band order changed")
    args.work.mkdir(parents=True,exist_ok=True);args.output.parent.mkdir(parents=True,exist_ok=True)
    out={"schema":"ttf_proceed70_plant_genetic_original_sites_CHELSA8_support_v01",
         "phenotype_rate_values_opened":False,"genetic_fitness_effect_fitted":False,
         "prior_TTF_heldout_response_reopened":False,
         "raw_climate_or_coordinates_publicly_exported":False}
    try:
        original_catalog=validated_get_catalog(
            "https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(original_catalog,src)
        raw,md5,sha256=fetch_original_datafile_bytes(1024206)
        if md5!=native["observed_md5"] or sha256!=native["observed_sha256"] or len(raw)!=native["observed_size_bytes"]:
            raise ValueError("Original plant data MD5/SHA256/bytes changed")
        cohort,studies=source_geometry(raw,contract)
        roles,split,digest=component_assignment(studies)
        expected=contract["authority"]
        if (len(roles)!=expected["plant_source_qualified_taxa"] or
            len(studies)!=expected["plant_source_studies"] or
            digest!=expected["plant_role_membership_sha256"]):
            raise ValueError("Frozen original plant study×species panel changed")
        coords=sorted({xy for _,_,a,b in cohort for xy in (a,b)})
        xy=np.asarray(coords,float)
        climate=[];receipts=[]
        for i,a in enumerate(asset["assets"]):
            path=args.work/f"climate_{i}.tif"
            verified_raster(a,path)
            value=raster_values(path,xy)
            climate.append(value)
            receipts.append({"name":labels[i],"sha256":a["sha256"],
                             "source_bytes":a["size_bytes"],
                             "finite_unique_original_sites":int(np.isfinite(value).sum())})
            path.unlink()
            print("CHELSA8_PLANT_SOURCE_VERIFIED",labels[i],int(np.isfinite(value).sum()),flush=True)
        out.update(source_climate_coverage(cohort,studies,climate,args.output))
        out["source_asset_receipts"]=receipts
        out["source_native_proceed_md5"]=md5
        out["source_native_proceed_sha256"]=sha256
        gates=contract["frozen_minimums"]
        passed=(out["climate_complete_species"]>=gates["complete_species_total"]
                and out["climate_complete_source_studies"]>=gates["complete_studies_total"]
                and out["climate_complete_genera"]>=gates["complete_genera_total"]
                and min(out["climate_complete_species_per_original_role"].values())>=
                    min(gates["complete_train_species"],gates["complete_heldout_species"]))
        out["status"]="PASS_CHELSA8_PLANT_SOURCE_SUPPORT_NO_ADAPTATION_TEST" if passed else gates["fail_status"]
    except Exception as e:
        out["status"]="NOT_EVALUABLE_ORIGINAL_CHELSA_SOURCE_OR_INPUT_REPLAY"
        out["error_type"]=type(e).__name__
        out["error_message"]=str(e)[:400]
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:out.get(k) for k in ("status","unique_original_georeferenced_lonlat_sites",
          "all_eight_climate_complete_pair_records","climate_complete_species",
          "climate_complete_source_studies","climate_complete_species_per_original_role")},sort_keys=True))


if __name__=="__main__":main()
