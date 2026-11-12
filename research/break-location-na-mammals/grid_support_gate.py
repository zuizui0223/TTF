"""Sample Pleistocene break points and within-survey background; report support counts only.

No genetic-climate response value or covariate is read, and no raw IUCN-derived point coordinates are exported.
"""
import argparse,json,math,zipfile
from pathlib import Path
import numpy as np
from pyproj import CRS,Transformer
from shapely.geometry import Point
from shapely.prepared import prep
from shapely.ops import transform
from shapely.validation import make_valid
from source_integrity import fetch
from source_overlap_gate import members,load_group
def count_sampled(line,poly,config):
    if not poly.is_valid:poly=make_valid(poly)
    anchor=poly.representative_point()
    local=CRS.from_proj4("+proj=aeqd +lat_0="+str(anchor.y)+
              " +lon_0="+str(anchor.x)+" +datum=WGS84 +units=m +no_defs")
    tr=Transformer.from_crs(4326,local,always_xy=True)
    line=transform(tr.transform,line);poly=transform(tr.transform,poly)
    source_line_length_m=line.length
    positive_step=float(config["positive_spacing_km"])*1000
    background_step=float(config["background_spacing_km"])*1000
    buffer_dist=float(config["exclusion_buffer_km"])*1000
    if source_line_length_m<=0 or poly.area<=0:return (0,0,0)
    contained=prep(poly)
    pos=0
    for d in np.arange(positive_step/2,source_line_length_m,positive_step):
        if contained.covers(line.interpolate(float(d))):pos+=1
    # An inference about negative sites is not justified close to an uncertain genetic break
    forbidden=prep(line.buffer(buffer_dist))
    minx,miny,maxx,maxy=poly.bounds
    xs=np.arange(math.ceil(minx/background_step)*background_step,maxx,background_step)
    ys=np.arange(math.ceil(miny/background_step)*background_step,maxy,background_step)
    if len(xs)*len(ys)>400000:raise ValueError("grid exceeds safety bound")
    bg=0
    for y in ys:
        for x in xs:
            p=Point(float(x),float(y))
            if contained.contains(p) and not forbidden.intersects(p):bg+=1
    return pos,bg,int(len(xs)*len(ys))
def main():
    p=argparse.ArgumentParser()
    for key in ("manifest","rule","sampling","out","work"):
        p.add_argument("--"+key,required=True,type=Path)
    a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True)
    manifest=json.loads(a.manifest.read_text());rule=json.loads(a.rule.read_text())
    conf=json.loads(a.sampling.read_text())
    if conf["status"]!="FROZEN_BEFORE_POSITIVE_BACKGROUND_SAMPLING":
        raise ValueError("sampling preregistration missing")
    if rule["status"]!="PRE_SPATIAL_COORDINATE_OPENING_SOURCE_QA":
        raise ValueError("spatial gate unverified")
    assets={v["role"]:v for v in manifest["files"]}
    b=a.work/"break.zip";r=a.work/"range.zip"
    fetch(assets["pleistocene_break_lines"],b)
    fetch(assets["surveyed_or_clipped_ranges"],r)
    out=[]
    with zipfile.ZipFile(b) as bz,zipfile.ZipFile(r) as rz:
        bb=members(bz,"break");rr=members(rz,"range")
        for k in sorted(set(bb)&set(rr)):
            name=bb[k][0]
            obj={"species_common_name":name}
            try:
                line,_=load_group(bz,bb[k][1])
                poly,_=load_group(rz,rr[k][1])
                if not poly.is_valid:poly=make_valid(poly)
                positive,background,candidates=count_sampled(line,poly,conf)
                obj.update({"positive_break_grid_points":positive,
                    "background_surveyed_grid_points":background,
                    "background_bbox_candidate_grid_points":candidates})
                obj["eligible"]=(positive>=conf["min_positive_points"]
                    and background>=conf["min_background_points"]
                    and background>=conf["min_background_per_positive"]*positive)
            except Exception as e:
                obj["eligible"]=False
                obj["error"]=type(e).__name__+":"+str(e)[:150]
            out.append(obj)
    accepted=[x for x in out if x["eligible"]]
    report={"schema":"ttf_jensen2024_within_survey_sampling_support_v01",
      "status":"PASS_SAMPLED_BREAK_BACKGROUND_SUPPORT" if len(accepted)>=conf["minimum_species_after_sampling"] else "NOT_EVALUABLE_SAMPLED_BREAK_BACKGROUND_SUPPORT",
      "paired_species_input":len(out),"qualified_species":len(accepted),
      "source_baseline_sample_policy":"10km breakpoint spacing; 20km within-survey grid; background excludes uncertain 50km break buffer",
      "total_break_points":sum(x["positive_break_grid_points"] for x in accepted),
      "total_surveyed_background_points":sum(x["background_surveyed_grid_points"] for x in accepted),
      "ineligible":[x for x in out if not x["eligible"]],
      "full_sampling_counts":out,
      "no_genetic_break_locations_or_IUCN_derived_background_coordinates_exported":True,
      "historic_climate_or_empirical_model_opened":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k not in ("full_sampling_counts","ineligible")},sort_keys=True))
if __name__=="__main__":main()
