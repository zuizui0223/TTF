"""Inspect georeferenced ArcInfo GRID layers inside UNEP-WCMC 2002 mountain archive.
The display TIFF in this archive is RGB and is NOT a mountain distance source.
No Jensen genetic breaks, case points, or backgrounds are read.
"""
import argparse,hashlib,json,urllib.request,zipfile
from pathlib import Path
import rasterio
from pyproj import Transformer

URL="https://datadownload-production.s3.amazonaws.com/Mountains_and_Tree_Cover_in_Mountain_Regions_2002.zip"
SHA="9a2c8910cc9b694ef967bac15994312b9518ecb071c262414d199ad4a3e66307"
BYTES=50857843
def verified_fetch(path):
    size=0;digest=hashlib.sha256()
    with urllib.request.urlopen(urllib.request.Request(URL,headers={"User-Agent":"TTF-WCMC-ArcGRID-geometry-audit"}),timeout=120) as response,path.open("wb") as file:
        while True:
            chunk=response.read(1048576)
            if not chunk:break
            size+=len(chunk);digest.update(chunk);file.write(chunk)
            if size>BYTES:raise ValueError("source exceeds frozen bytes")
    if size!=BYTES or digest.hexdigest()!=SHA:raise ValueError("UNEP original bytes mismatched")
    return size
def inspect_grid(path):
    try:
        with rasterio.open(path) as ds:
            report={"status":"OPEN_AS_NATIVE_GIS_GRID","driver":ds.driver,
                "crs":str(ds.crs),"bands":ds.count,"width":ds.width,"height":ds.height,
                "resolution":[float(v) for v in ds.res],
                "bounds":[float(x) for x in ds.bounds],
                "dtypes":list(ds.dtypes),"nodata":ds.nodata}
            if ds.crs:
                tr=Transformer.from_crs("EPSG:4326",ds.crs,always_xy=True)
                refs={"Montana USA":(-110.0,46.0),"Florida USA":(-80.5,26.0),
                      "Himalaya":(86.0,28.0)}
                report["non_genetic_reference_samples"]={}
                for name,(lon,lat) in refs.items():
                    x,y=tr.transform(lon,lat)
                    value=next(ds.sample([(x,y)]))
                    report["non_genetic_reference_samples"][name]=[float(v) for v in value]
            return report
    except Exception as e:
        return {"status":"GRID_OPEN_FAILED","error_type":type(e).__name__,"error":str(e)[:220]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args();args.work.mkdir(parents=True,exist_ok=True);args.out.parent.mkdir(parents=True,exist_ok=True)
    zippath=args.work/"WCMC2002.zip";verified_fetch(zippath)
    with zipfile.ZipFile(zippath) as z:
        files=z.namelist()
        if any(x.startswith("/") or ".." in Path(x).parts for x in files):
            raise ValueError("unsafe WCMC ZIP paths")
        grids=sorted({x.rsplit("/",1)[0] for x in files if x.endswith("/hdr.adf")})
        if not grids:raise ValueError("no ArcInfo grids")
        z.extractall(args.work/"extracted")
    entries=[]
    for g in grids:
        info=inspect_grid(args.work/"extracted"/g)
        entries.append({"grid_directory":g,**info})
        print("MOUNTAIN ARCGIS GRID",g,info["status"],info.get("driver"),info.get("crs"),flush=True)
    report={"schema":"ttf_wcmc2002_original_ArcInfo_mountain_grid_schema_v01",
        "status":"NATIVE_ARCINFO_MOUNTAIN_LAYERS_AUDITED_NO_BREAK_RESPONSE",
        "source_sha256":SHA,"grids":entries,
        "any_valid_georeferenced_grid":any(x.get("status")=="OPEN_AS_NATIVE_GIS_GRID" and x.get("crs")!="None" for x in entries),
        "do_not_use_ungeoreferenced_RGB_TIF_as_mountain_indicator":True,
        "unep_nonredistribution_complied":True,
        "empirical_break_AUC_tested":False}
    args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"grid_count":len(entries),"any_georeferenced":report["any_valid_georeferenced_grid"],"empirical_model":False}))
if __name__=="__main__":main()
