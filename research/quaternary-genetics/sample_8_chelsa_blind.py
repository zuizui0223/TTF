"""Sample eight original C1 CHELSA rasters at He-blind terrestrial population sites."""
from __future__ import annotations
import argparse, collections, csv, hashlib, json, math, os
from pathlib import Path
import numpy as np
import rasterio
import requests

def download_one(asset, path):
    h=hashlib.sha256();size=0
    with requests.get(asset["url"],stream=True,timeout=(30,300)) as response:
        response.raise_for_status()
        with path.open("wb") as f:
            for chunk in response.iter_content(2**20):
                if chunk:
                    size+=len(chunk);h.update(chunk);f.write(chunk)
    if size!=int(asset["size_bytes"]) or h.hexdigest()!=asset["sha256"]:
        path.unlink(missing_ok=True)
        raise ValueError("original-C1-asset SHA/size failure")
    return size

def sample(exact_geo,manifest,out):
    assets=manifest["assets"]
    if len(assets)!=8 or {a["variable"] for a in assets}!={"bio01","bio07","bio12","bio15"}:
        raise ValueError("wrong CHELSA assets")
    sites=exact_geo["adult_terrestrial_climate_sites"]
    if not sites or not all(all(k in s for k in ("species_id","study_id","latitude","longitude")) for s in sites):
        raise ValueError("missing or changed coordinate manifest")
    groups={(m["species_id"],m["study_id"]) for m in exact_geo["adult_terrestrial_frozen_manifest"]}
    if any((s["species_id"],s["study_id"]) not in groups for s in sites):
        raise ValueError("outside frozen prospective eligible groups")
    points=[(float(r["longitude"]),float(r["latitude"])) for r in sites]
    target=out/"downloads";target.mkdir(parents=True,exist_ok=True)
    v=np.empty((len(points),8),float)
    receipts=[]
    for i,asset in enumerate(assets):
        path=target/f"original_{i}.tif"
        size=download_one(asset,path)
        with rasterio.open(path) as ds:
            if ds.crs is None or not ds.crs.is_geographic:raise ValueError("non-geographic raster")
            vals=np.asarray([float(x[0]) for x in ds.sample(points)],dtype=float)
            if ds.nodata is not None:vals[vals==float(ds.nodata)]=np.nan
            if np.isinf(vals).any():raise ValueError("infinite raster climate")
            v[:,i]=vals
        path.unlink(missing_ok=True)
        receipts.append({"variable":asset["variable"],"time_index":asset["time_index"],
              "sha256":asset["sha256"],"size_bytes":size,
              "finite_site_values":int(np.isfinite(v[:,i]).sum())})
        print("ORIGINAL RASTER VERIFIED AND SAMPLED",asset["variable"],asset["time_index"],size,flush=True)
    valid=np.isfinite(v).all(axis=1)
    with (out/"site_keys.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.writer(f);w.writerow(["species_id","study_id","latitude","longitude"])
        for x in sites:w.writerow([x["species_id"],x["study_id"],x["latitude"],x["longitude"]])
    np.save(out/"climate_values.npy",v,allow_pickle=False)
    bykey=collections.Counter((s["species_id"],s["study_id"]) for s,b in zip(sites,valid) if b)
    retained={k for k,n in bykey.items() if n>=5}
    receipt={"schema":"gendivrange_original_C1_eight_rasters_terrestrial_response_blind_v01",
       "status":"REAL_SAMPLED_CLIMATE_NO_GENETIC_RESPONSE",
       "sites":len(sites),"sites_all_eight_climates_finite":int(valid.sum()),
       "qualified_groups_with_five_finite_sites":len(retained),
       "input_exact_geographic_qualified_groups":len(groups),
       "climate_input_sha256":hashlib.sha256(json.dumps(sites,sort_keys=True).encode()).hexdigest(),
       "climate_matrix_sha256":hashlib.sha256(v.tobytes()).hexdigest(),
       "climate_site_keys_sha256":hashlib.sha256((out/"site_keys.csv").read_bytes()).hexdigest(),
       "raster_source_receipts":receipts,"genetic_He_values_opened":False,
       "genetic_response_authorized":False}
    (out/"climate_receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:x for k,x in receipt.items() if k!="raster_source_receipts"},sort_keys=True))
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geo",type=Path,required=True)
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    sample(json.loads(a.geo.read_text()),json.loads(a.manifest.read_text()),a.out)
