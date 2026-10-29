"""Jensen2024 break-line geometry vs clipped surveyed range pre-qualification.
Pure spatial QC only; does not read any paleoclimate or fit genetic models.
"""
import argparse,io,json,math,re,zipfile
from pathlib import Path
import shapefile
from pyproj import CRS,Transformer
from shapely.geometry import shape
from shapely.ops import transform,unary_union
from shapely.validation import make_valid
from source_integrity import fetch
def norm(s):return re.sub("[^a-z0-9]","",s.lower())
def members(z,role):
    out={}
    for name in z.namelist():
        if not name.endswith(".shp"):continue
        parts=Path(name).parts
        common=parts[2] if role=="break" else parts[1]
        out[norm(common)]=(common,name[:-4])
    return out
def load_group(z,stem):
    reader=shapefile.Reader(shp=io.BytesIO(z.read(stem+".shp")),
           shx=io.BytesIO(z.read(stem+".shx")),
           dbf=io.BytesIO(z.read(stem+".dbf")))
    crs=CRS.from_wkt(z.read(stem+".prj").decode("utf-8-sig","replace"))
    geoms=[]
    for f in reader.iterShapes():
        try:
            g=shape(f.__geo_interface__)
            if not g.is_empty:geoms.append(g)
        except Exception:continue
    if not geoms:raise ValueError("No readable source geometry")
    geom=unary_union(geoms)
    to_wgs=Transformer.from_crs(crs,CRS.from_epsg(4326),always_xy=True)
    return transform(to_wgs.transform,geom),len(geoms)
def qualify(bzip,rzip,contract):
    if contract.get("status")!="PRE_SPATIAL_COORDINATE_OPENING_SOURCE_QA":
        raise ValueError("unfrozen geographical rule")
    with zipfile.ZipFile(bzip) as bz,zipfile.ZipFile(rzip) as rz:
        br=members(bz,"break");ra=members(rz,"range")
        joined=sorted(set(br)&set(ra))
        result=[]
        for key in joined:
            common,stem=br[key]
            rec={"common_name":common}
            try:
                line,nline=load_group(bz,stem)
                poly,npoly=load_group(rz,ra[key][1])
                repaired=False
                if not poly.is_valid:
                    poly=make_valid(poly)
                    repaired=True
                if poly.is_empty or line.is_empty or poly.geom_type not in ("Polygon","MultiPolygon"):
                    raise ValueError("empty or non-polygon range geometry")
                ref=poly.representative_point()
                if not (-90<ref.y<90 and -180<ref.x<180):
                    raise ValueError("invalid range coordinate")
                local=CRS.from_proj4("+proj=aeqd +lat_0="+str(ref.y)+
                    " +lon_0="+str(ref.x)+" +datum=WGS84 +units=m +no_defs")
                transformer=Transformer.from_crs(4326,local,always_xy=True)
                l=transform(transformer.transform,line)
                p=transform(transformer.transform,poly)
                length=float(l.length/1000)
                distance=float(l.distance(p)/1000)
                if not (math.isfinite(length) and math.isfinite(distance)):
                    raise ValueError("nonfinite geometry metrics")
                valid=(length>0 and p.area>0 and distance<=50)
                rec.update({"line_source_features":nline,"range_source_features":npoly,
                            "line_length_km":length,
                            "break_to_surveyed_range_nearest_distance_km":distance,
                            "topology_repaired":repaired,
                            "within_predeclared_50km_buffer":valid})
            except Exception as e:
                rec.update({"within_predeclared_50km_buffer":False,
                            "error":type(e).__name__+": "+str(e)[:130]})
            result.append(rec)
    ok=sum(v["within_predeclared_50km_buffer"] for v in result)
    return {"schema":"ttf_jensen2024_spatial_break_range_overlap_gate_v01",
        "status":"PASS_SPATIAL_SOURCE_SUPPORT" if ok>=contract["minimum_species"] else "NOT_EVALUABLE_SPATIAL_SOURCE_SUPPORT",
        "file_name_joined_species":len(joined),"spatially_qualified_species":ok,
        "predeclared_required_species":contract["minimum_species"],
        "spatially_disqualified":[v for v in result if not v["within_predeclared_50km_buffer"]],
        "topology_repair_count":sum(v.get("topology_repaired",False) for v in result),
        "spatial_distances_km_summary":{
          "min":min((v["break_to_surveyed_range_nearest_distance_km"] for v in result if "break_to_surveyed_range_nearest_distance_km" in v),default=None),
          "max":max((v["break_to_surveyed_range_nearest_distance_km"] for v in result if "break_to_surveyed_range_nearest_distance_km" in v),default=None)},
        "jensen_break_locations_opened_for_spatial_QA_only":True,
        "historic_climate_predictors_tested":False,
        "genetic_break_response_model_fit":False}
def main():
    a=argparse.ArgumentParser()
    for n in ("manifest","rule","out","work"):a.add_argument("--"+n,required=True,type=Path)
    q=a.parse_args()
    m=json.loads(q.manifest.read_text());rule=json.loads(q.rule.read_text())
    if not rule["frozen_before_break_geometry_coordinates"]:raise ValueError("unfrozen")
    f={x["role"]:x for x in m["files"]}
    q.work.mkdir(parents=True,exist_ok=True)
    b=q.work/"breaks.zip";r=q.work/"surveyed-ranges.zip"
    fetch(f["pleistocene_break_lines"],b)
    fetch(f["surveyed_or_clipped_ranges"],r)
    report=qualify(b,r,rule)
    q.out.parent.mkdir(parents=True,exist_ok=True)
    q.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True))
if __name__=="__main__":main()
