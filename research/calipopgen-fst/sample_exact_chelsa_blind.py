"""Sample frozen CHELSA 21ka/0BP at CaliPopGen Fst-qualified sites, response-blind."""
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np
import rasterio,requests

def acquire(asset,path):
    h=hashlib.sha256();size=0
    with requests.get(asset["url"],stream=True,timeout=(30,360)) as response:
        response.raise_for_status()
        with path.open("wb") as out:
            for chunk in response.iter_content(1024*1024):
                if chunk:
                    size+=len(chunk);h.update(chunk);out.write(chunk)
    if size!=int(asset["size_bytes"]) or h.hexdigest()!=asset["sha256"]:
        path.unlink(missing_ok=True)
        raise RuntimeError("Original C1 climate raster hash/size mismatch")
    return size

def run(q,assets,out):
    if q["status"]!="PASS_RESPONSE_BLIND_GEOGRAPHIC_SUPPORT" or q["Fst_numeric_values_opened"]:
        raise ValueError("No response-blind qualifying cohort")
    if len(assets["assets"])!=8 or assets["Fst_values_opened"]:
        raise ValueError("Incorrect exact CHELSA asset contract")
    groups=q["eligible_groups"]
    coords=sorted({tuple(map(float,c)) for g in groups for c in g["site_coordinates"]})
    if not coords:raise ValueError("no sites")
    loc={p:i for i,p in enumerate(coords)}
    values=np.full((len(coords),8),np.nan,dtype=np.float64)
    downloads=out/"_temp_COGs";downloads.mkdir(parents=True,exist_ok=True)
    receipts=[]
    for i,asset in enumerate(assets["assets"]):
        path=downloads/f"original_{i}.tif"
        byte_count=acquire(asset,path)
        with rasterio.open(path) as ds:
            if ds.crs is None or not ds.crs.is_geographic:
                raise ValueError("CHELSA is not geographical")
            samples=np.array([float(x[0]) for x in ds.sample([(lon,lat) for lat,lon in coords])])
            if ds.nodata is not None:samples[samples==float(ds.nodata)]=np.nan
            samples[~np.isfinite(samples)]=np.nan
            values[:,i]=samples
        path.unlink(missing_ok=True)
        receipts.append({k:asset[k] for k in ("variable","time_index","sha256","size_bytes")})
        print("REAL CHELSA FILE HASH VERIFIED",i+1,asset["variable"],asset["time_index"],byte_count,flush=True)
    finite=np.isfinite(values).all(axis=1)
    accepted=[]
    for g in groups:
        valid=[p for p in g["pair_coordinates"] if finite[loc[(p[0],p[1])]] and finite[loc[(p[2],p[3])]]]
        sites={tuple(p[:2]) for p in valid}|{tuple(p[2:]) for p in valid}
        if len(valid)<5 or len(sites)<4:continue
        accepted.append({k:v for k,v in g.items() if k not in ("pair_coordinates","site_coordinates") }|
                        {"climate_valid_pair_coordinates":valid,"climate_valid_pair_count":len(valid)})
    species_by_role=[{g["species"] for g in accepted if g["role"]==i} for i in (0,1)]
    groups_by_role=[sum(g["role"]==i for g in accepted) for i in (0,1)]
    nsp=len(species_by_role[0]|species_by_role[1])
    passed=(nsp>=30 and len(accepted)>=35 and min(len(v) for v in species_by_role)>=10 and min(groups_by_role)>=12)
    status="PASS_CLIMATE_RESPONSE_BLIND" if passed else "NOT_EVALUABLE_AFTER_CLIMATE_MASK"
    with (out/"site_keys.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["latitude","longitude"])
        w.writerows(coords)
    np.save(out/"climate_values.npy",values,allow_pickle=False)
    receipt={"schema":"ttf_calipopgen_fst_climate_sample_v01",
       "status":status,"unique_georeferenced_sites":len(coords),
       "sites_with_eight_climate_values":int(finite.sum()),
       "eligible_groups_before_climate":len(groups),
       "eligible_groups_after_climate":len(accepted),
       "eligible_species_after_climate":nsp,
       "eligible_species_per_role":[len(v) for v in species_by_role],
       "eligible_groups_per_role":groups_by_role,
       "climate_value_matrix_sha256":hashlib.sha256(values.tobytes()).hexdigest(),
       "site_order_sha256":hashlib.sha256((out/"site_keys.csv").read_bytes()).hexdigest(),
       "original_CHELSA_assets":receipts,"Fst_values_opened":False}
    (out/"climate_receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    (out/"eligible_pair_coordinates.json").write_text(json.dumps(
      {"schema":"ttf_calipopgen_fst_qualified_pair_coordinates_with_climate_v01",
       "status":status,"groups":accepted,"Fst_values_opened":False},indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k!="original_CHELSA_assets"},sort_keys=True))
    if not passed:print("STOP BEFORE GENETIC OUTCOME OPENING: climate cohort failed")
if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--qualified",type=Path,required=True)
    a.add_argument("--assets",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    p=a.parse_args();p.out.mkdir(parents=True,exist_ok=True)
    run(json.loads(p.qualified.read_text()),json.loads(p.assets.read_text()),p.out)
