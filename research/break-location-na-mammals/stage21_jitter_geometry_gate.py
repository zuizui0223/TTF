"""Audit exact Jensen 80-species 50km positional jitter, no climate or break model.

All author/IUCN geometries and sampled positive-point coordinates stay on the
ephemeral Actions runner. Output is counts and displacement magnitudes only.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,shapely
from pyproj import CRS,Transformer
from shapely.ops import transform
from shapely.validation import make_valid
from source_integrity import fetch
from source_overlap_gate import members,load_group,norm
from jitter_case_geometry import jitter_cases,repair_survey_polygon
import zipfile

def main():
    ap=argparse.ArgumentParser()
    for k in ("panel","sources","contract","work","output"):
        ap.add_argument("--"+k,required=True,type=Path)
    a=ap.parse_args()
    panel=json.loads(a.panel.read_text())
    contract=json.loads(a.contract.read_text())
    if contract["status"]!="LOCKED_PRE_EMPIRICAL_BREAK_AUC":
        raise ValueError("Jitter protocol not frozen")
    rules=contract["case_jitter"]
    a.work.mkdir(parents=True,exist_ok=True);a.output.parent.mkdir(parents=True,exist_ok=True)
    source={x["role"]:x for x in json.loads(a.sources.read_text())["files"]}
    b=a.work/"original_breaks.zip";r=a.work/"original_ranges.zip"
    fetch(source["pleistocene_break_lines"],b)
    fetch(source["surveyed_or_clipped_ranges"],r)
    records=[x for x in panel["qualified_species_records"] if x["synthetic_estimator_min10_cases"]]
    if len(records)!=80:raise ValueError("wrong 80-species panel")
    summaries=[];all_distances=[];n_cases=0;n_moved=0;roles=[0,0]
    with zipfile.ZipFile(b) as bz,zipfile.ZipFile(r) as rz:
        bb=members(bz,"break");rr=members(rz,"range")
        for i,s in enumerate(records):
            k=norm(s["common_name"])
            line,_=load_group(bz,bb[k][1]);poly,_=load_group(rz,rr[k][1])
            if not poly.is_valid:poly=make_valid(poly)
            anchor=poly.representative_point()
            crs=CRS.from_proj4("+proj=aeqd +lat_0="+str(anchor.y)+
                " +lon_0="+str(anchor.x)+" +datum=WGS84 +units=m +no_defs")
            tr=Transformer.from_crs(4326,crs,always_xy=True)
            lin=transform(tr.transform,line)
            polygon=transform(tr.transform,poly)
            distances=np.arange(5000.,lin.length,10000.)
            points=shapely.line_interpolate_point(lin,distances)
            points=points[shapely.covers(polygon,points)]
            xy=shapely.get_coordinates(points)
            if len(xy)!=s["case_points"]:
                raise ValueError("Original case-grid sites differ from frozen source, species index "+str(i))
            seed_key=rules["seed_namespace"]+s["scientific_name"].lower()
            seed=int.from_bytes(hashlib.sha256(seed_key.encode()).digest()[:8],"big")
            jitter_polygon=repair_survey_polygon(polygon)
            moved,diagnostic=jitter_cases(xy,jitter_polygon,np.random.default_rng(seed),
                radius_m=rules["source_uncertainty_radius_km"]*1000,
                max_attempts=rules["max_attempts"])
            displacement=np.linalg.norm(moved-xy,axis=1)/1000.
            if (displacement>rules["source_uncertainty_radius_km"]+1e-8).any() or not shapely.covers(jitter_polygon,shapely.points(moved)).all():
                raise RuntimeError("Jitter invalid: source range or radius violated")
            roles[int(s["split_role"])]+=1
            n_cases+=len(xy);n_moved+=diagnostic["successful"];all_distances.append(displacement)
            summaries.append({"scientific_name":s["scientific_name"],"role":int(s["split_role"]),
              "case_grid_points":len(xy),**diagnostic})
            if (i+1)%10==0:print("SOURCE_JITTER_GEO_VERIFIED",i+1,flush=True)
    displacements=np.concatenate(all_distances)
    median_nonzero=float(np.median(displacements[displacements>0])) if (displacements>0).any() else 0.
    passed=(n_cases==contract["panel"]["positive_case_site_records"] and roles==[44,36]
      and n_moved/n_cases>=rules["nonzero_required_fraction"]
      and median_nonzero>=rules["global_nonzero_median_km_at_least"])
    out={"schema":"ttf_jensen80_break_50km_position_jitter_geometry_gate_v01",
      "status":"PASS_GEOMETRIC_50KM_JITTER_FEASIBILITY" if passed else "NOT_EVALUABLE_50KM_POSITIONAL_JITTER",
      "species":len(records),"species_roles":roles,"original_positive_site_records":n_cases,
      "jittered_case_sites":n_moved,"original_site_fallbacks":n_cases-n_moved,
      "jitter_acceptance_fraction":n_moved/n_cases,
      "nonzero_shift_median_km":median_nonzero,
      "displacement_p95_km":float(np.quantile(displacements,.95)),
      "rules_hash_sha256":hashlib.sha256(a.contract.read_bytes()).hexdigest(),
      "per_species_aggregate":summaries,
      "case_or_range_coordinates_exported":False,
      "jittered_climate_or_barrier_values_extracted":False,
      "genetic_history_auc_evaluated":False,
      "prior_Fst_He_COI_outcomes_reopened":False}
    a.output.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k!="per_species_aggregate"},sort_keys=True))
if __name__=="__main__":main()
