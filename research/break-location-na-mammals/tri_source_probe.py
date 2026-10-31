"""Check real EarthEnv TRI raster metadata and SHA, without accessing genetic breaks."""
import argparse,hashlib,json,urllib.error,urllib.request
from pathlib import Path

def inspect(candidate,work):
    url=candidate["candidate_url"]
    if url!="https://data.earthenv.org/topography/tri_1KMmd_GMTEDmd.tif":
        raise ValueError("unapproved TRI URI")
    maxsize=int(candidate["max_probe_bytes"])
    r={"source":url,"status":"UNKNOWN","SHA256_verified":False,
       "break_geometry_or_genetic_predictor_read":False}
    req=urllib.request.Request(url,method="HEAD",headers={"User-Agent":"TTF-mammal-break-TRI-source-audit-v01"})
    try:
        with urllib.request.urlopen(req,timeout=45) as h:
            r["head_status"]=h.status
            r["content_type"]=h.headers.get("Content-Type")
            length=h.headers.get("Content-Length")
            r["content_length"]=int(length) if length and length.isdigit() else None
            r["range_support"]=h.headers.get("Accept-Ranges")
    except Exception as e:
        r["head_error"]=f"{type(e).__name__}: {str(e)[:180]}"
    if r.get("content_length",maxsize+1)>maxsize:
        r["status"]="LARGE_OR_UNKNOWN_FILE_NOT_DOWNLOADED"
        return r
    target=work/"tri_source.tif"
    n=0;digest=hashlib.sha256()
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"TTF-TRI-source-audit-v01"}),timeout=100) as resp,target.open("wb") as file:
            while 1:
                block=resp.read(1024*1024)
                if not block:break
                n+=len(block)
                if n>maxsize:raise ValueError("Exceeded frozen size cap")
                file.write(block);digest.update(block)
        r.update({"downloaded_bytes":n,"sha256":digest.hexdigest(),"SHA256_verified":True})
        import rasterio
        with rasterio.open(target) as ds:
            r["raster"]={"width":ds.width,"height":ds.height,"bands":ds.count,
                "crs":str(ds.crs),"bounds":[float(x) for x in ds.bounds],
                "nodata":str(ds.nodata),"dtype":ds.dtypes[0],
                "pixel_size":[float(x) for x in ds.res]}
        r["status"]="TRI_SOURCE_HASHED_AND_GEOGRAPHIC_SCHEMA_VALID"
    except Exception as e:
        r["status"]="TRI_SOURCE_DOWNLOAD_OR_RASTER_SCHEMA_FAILED"
        r["downloaded_bytes"]=n
        r["error"]=f"{type(e).__name__}: {str(e)[:220]}"
    finally:
        target.unlink(missing_ok=True)
    return r

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--candidate",type=Path,required=True)
    parser.add_argument("--work",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    a.work.mkdir(parents=True,exist_ok=True);a.out.parent.mkdir(parents=True,exist_ok=True)
    cfg=json.loads(a.candidate.read_text())
    result={"schema":"ttf_jensen2024_tri1km_exact_source_receipt_v01",**inspect(cfg,a.work),
      "source_vintage_2018_aggregation_algorithm_equivalence":False,
      "empirical_break_location_model_tested":False}
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
