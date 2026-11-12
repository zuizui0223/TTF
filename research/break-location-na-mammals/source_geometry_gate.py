"""Response-blind geometry schema and source taxon key qualification, Jensen2024.
Reads Shapefile headers/CRS and counts, but never evaluates genetic break positions.
"""
import argparse,collections,io,json,re,zipfile
from pathlib import Path
import shapefile
from pyproj import CRS
from source_integrity import fetch

def norm(s):return re.sub(r"[^a-z0-9]","",s.lower())
def shapebundles(z,kind):
    out={}
    for name in z.namelist():
        if not name.lower().endswith(".shp"):continue
        parts=Path(name).parts
        species=parts[2] if kind=="break" else parts[1]
        out[norm(species)]={"common_name":species,"shp":name}
    return out
def describe(z,item,expect):
    stem=item["shp"][:-4]
    for ext in (".shp",".shx",".dbf",".prj"):
        if stem+ext not in z.namelist():
            return {"valid":False,"issue":"missing "+ext}
    try:
        reader=shapefile.Reader(shp=io.BytesIO(z.read(stem+".shp")),
                    shx=io.BytesIO(z.read(stem+".shx")),
                    dbf=io.BytesIO(z.read(stem+".dbf")))
        geom=reader.shapeTypeName
        n=len(reader)
        crs=CRS.from_wkt(z.read(stem+".prj").decode("utf-8-sig","replace"))
        epsg=crs.to_epsg()
        expected=({"POLYLINE","POLYLINEZ","POLYLINEM"} if expect=="break"
                   else {"POLYGON","POLYGONZ","POLYGONM"})
        if not any(t in geom.upper() for t in expected):
            return {"valid":False,"issue":"unexpected "+geom,"epsg":epsg,"features":n}
        if n<1:return {"valid":False,"issue":"no geometry features","epsg":epsg}
        return {"valid":True,"shape_type":geom,"features":n,
            "crs_epsg":epsg,"crs_is_geographic":crs.is_geographic}
    except Exception as e:
        return {"valid":False,"issue":type(e).__name__+":"+str(e)[:100]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True)
    m=json.loads(a.manifest.read_text())
    files={x["role"]:x for x in m["files"]}
    paths={}
    for name in ("pleistocene_break_lines","surveyed_or_clipped_ranges"):
        e=files[name];p=a.work/(str(e["id"])+".zip")
        fetch(e,p);paths[name]=p
    with zipfile.ZipFile(paths["pleistocene_break_lines"]) as bz,zipfile.ZipFile(paths["surveyed_or_clipped_ranges"]) as rz:
        breaks=shapebundles(bz,"break");ranges=shapebundles(rz,"range")
        paired=sorted(set(breaks)&set(ranges))
        results=[]
        for key in paired:
            lb=describe(bz,breaks[key],"break")
            rr=describe(rz,ranges[key],"range")
            results.append({"species_common_name":breaks[key]["common_name"],
              "breaks":lb,"ranges":rr,
              "schema_qualified":lb["valid"] and rr["valid"]})
    good=[v for v in results if v["schema_qualified"]]
    report={"schema":"ttf_jensen2024_break_range_schema_qualification_v01",
      "status":"GEOMETRIC_SCHEMA_PREQUALIFICATION_ONLY",
      "break_archive_species":len(breaks),"range_archive_species":len(ranges),
      "exact_common_name_joined_species":len(paired),
      "qualified_source_schema_species":len(good),
      "unmatched_break_names":[breaks[k]["common_name"] for k in sorted(set(breaks)-set(ranges))],
      "shape_type_break":dict(collections.Counter(x["breaks"].get("shape_type","invalid") for x in results)),
      "shape_type_ranges":dict(collections.Counter(x["ranges"].get("shape_type","invalid") for x in results)),
      "crs_break":dict(collections.Counter(str(x["breaks"].get("crs_epsg")) for x in results)),
      "crs_ranges":dict(collections.Counter(str(x["ranges"].get("crs_epsg")) for x in results)),
      "unqualified_records":[x for x in results if not x["schema_qualified"]],
      "source_geometry_attribute_or_line_position_evaluated":False,
      "empirical_history_model_evaluated":False,
      "meets_frozen_minimum_50_species":len(good)>=50}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k not in ("unqualified_records",)},sort_keys=True))
if __name__=="__main__":main()
