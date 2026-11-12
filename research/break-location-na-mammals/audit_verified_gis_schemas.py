"""Audit provider GIS *schema* under frozen SHA, without genetic response access."""
import argparse,hashlib,io,json,urllib.request,zipfile
from pathlib import Path
from urllib.parse import urlparse
import shapefile
import rasterio
from pyproj import CRS

ALLOWED={"dmap-prod-oms-edc.s3.us-east-1.amazonaws.com","naturalearth.s3.amazonaws.com","datadownload-production.s3.amazonaws.com"}

def get_verified(item,sha_entry,work):
    u=item["url"]
    if urlparse(u).hostname not in ALLOWED:raise ValueError("source host changed")
    path=work/(item["term"]+".zip");h=hashlib.sha256();n=0
    with urllib.request.urlopen(urllib.request.Request(u,headers={"User-Agent":"TTF-GIS-source-schema-v01"}),timeout=120) as resp,path.open("wb") as f:
        while 1:
            block=resp.read(1024*1024)
            if not block:break
            n+=len(block);h.update(block);f.write(block)
            if n>350_000_000:raise ValueError("source oversized")
    if h.hexdigest()!=sha_entry["sha256"] or n!=sha_entry["bytes"]:
        path.unlink(missing_ok=True)
        raise ValueError("publisher GIS candidate SHA/size changed")
    return path

def shapefile_schema(z,term):
    shp=next(n for n in z.namelist() if n.lower().endswith(".shp"))
    stem=shp[:-4]
    opts={}
    for ext in (".shp",".shx",".dbf"):
        opts[ext.lstrip(".")]=io.BytesIO(z.read(stem+ext))
    rdr=shapefile.Reader(**opts)
    names=[f[0] for f in rdr.fields if f[0]!="DeletionFlag"]
    prj=stem+".prj"
    crs=None
    if prj in z.namelist():
        wkt=z.read(prj).decode("utf8",errors="replace")
        try:
            CRS.from_wkt(wkt)
            crs="VALID_WKT:"+str(CRS.from_wkt(wkt).to_epsg())
        except Exception:crs="CRS_NOT_PARSED"
    return {"source_term":term,"file_type":"ESRI_SHAPEFILE",
            "feature_count":len(rdr),"shape_type":rdr.shapeTypeName,
            "dbf_fields":names,
            "source_crs":crs,"schema_only":True}

def unep_schema(z,work):
    names=z.namelist()
    tif=next(n for n in names if n.endswith("/Mountain_UNEP-WCMC_2002/mountain_2002.tif"))
    info=z.getinfo(tif)
    out={"source_term":"unep_mountains","file_type":"ARCGIS_GRID_PLUS_TIFF_NOT_SHP",
        "original_tiff_uncompressed_bytes":info.file_size,
        "tiff_path":tif,"arcinfo_grid_files":sum(n.lower().endswith(".adf") for n in names)}
    target=work/"unep_mountain_mask.tif"
    with z.open(tif) as r,target.open("wb") as f:
        while 1:
            block=r.read(4*1024*1024)
            if not block:break
            f.write(block)
    try:
        with rasterio.open(target) as ds:
            out.update({"raster_crs":str(ds.crs),"raster_width":ds.width,
                "raster_height":ds.height,"raster_dtype":list(ds.dtypes),
                "raster_nodata":ds.nodata,
                "raster_resolution":[float(v) for v in ds.res],
                "raster_bounds":[float(v) for v in ds.bounds],
                "raster_band_count":ds.count})
            if ds.count>0 and ds.crs is not None and ds.crs.is_geographic:
                sample={"Montana":(-110.0,46.0),"Florida":(-80.5,26.0),
                        "Himalaya":(86.0,28.0)}
                out["synthetic_reference_pixel_values"]={name:[float(v) for v in next(ds.sample([p]))]
                      for name,p in sample.items()}
    finally:target.unlink(missing_ok=True)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sources",type=Path,required=True)
    ap.add_argument("--lock",type=Path,required=True)
    ap.add_argument("--work",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args();s=json.loads(a.sources.read_text());lock=json.loads(a.lock.read_text())
    expected={v["term"]:v for v in lock["asset_hashes"]}
    a.work.mkdir(parents=True,exist_ok=True);a.out.parent.mkdir(parents=True,exist_ok=True)
    reports=[]
    for item in s["source_candidates"]:
        zippath=get_verified(item,expected[item["term"]],a.work)
        with zipfile.ZipFile(zippath) as z:
            detail=unep_schema(z,a.work) if item["term"]=="unep_mountains" else shapefile_schema(z,item["term"])
        zippath.unlink()
        detail["verified_publisher_sha256"]=expected[item["term"]]["sha256"]
        detail["licensing_notes"]=item["rights"]
        reports.append(detail)
        print("SOURCE GEO SCHEMA VERIFIED",detail["source_term"],detail["file_type"],flush=True)
    result={"schema":"ttf_jensen_publisher_barrier_GIS_schema_v01",
            "status":"FIVE_VERIFIED_GIS_FORMATS_IDENTIFIED_NO_GENETIC_OUTCOME",
            "sources":reports,"tri_source":"EarthEnv TRI receipt verified separately",
            "break_or_background_coordinates_opened":False,"genetic_model_run":False}
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"sources":[{"name":v["source_term"],
        "fields":len(v.get("dbf_fields",[])),"features":v.get("feature_count"),
        "crs":v.get("source_crs") or v.get("raster_crs")} for v in reports]}))
if __name__=="__main__":main()
