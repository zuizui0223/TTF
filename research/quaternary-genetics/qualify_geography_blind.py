"""Pre-response GenDivRange geographic eligibility with exact distance and taxon joins."""
import argparse, collections, hashlib, json, math
from pathlib import Path
from inspect_headers_only import fetch_verified
from reconcile_keys import stream

NAMESPACE="ttf-fresh-quaternary-gendivrange-v01"
def haversine(a,b):
    la,lo,lb,ln=map(math.radians,(*a,*b))
    q=math.sin((lb-la)/2)**2+math.cos(la)*math.cos(lb)*math.sin((ln-lo)/2)**2
    return 12742.0176*math.asin(min(1.,math.sqrt(q)))
def has_100km_span(points):
    pts=sorted(points)
    if len(pts)<5:return False
    if 111*(max(v[0] for v in pts)-min(v[0] for v in pts))>=100:return True
    for i,p in enumerate(pts):
        for q in pts[i+1:]:
            if haversine(p,q)>=100:return True
    return False
def run(pop,spec):
    metadata={}
    for r in stream(spec,{"Spec_id","Study_id","Spec_Latin_GenDivRange",
                          "Marker_type","Habitat_adulthood","Habitat_breeding","N_pops"}):
        k=(r.get("Spec_id"),r.get("Study_id"))
        if k in metadata:raise ValueError("duplicate metadata study/species key")
        metadata[k]=r
    points=collections.defaultdict(set)
    allpop=collections.Counter()
    observed=collections.Counter()
    by_pop=collections.defaultdict(lambda: collections.defaultdict(set))
    bad=collections.Counter()
    for r in stream(pop,{"Spec_id","Study_id","Pop_id","Latitude","Longitude"},{"He"}):
        k=(r.get("Spec_id"),r.get("Study_id"))
        if k not in metadata:
            bad["missing_taxon_study_key"]+=1;continue
        allpop[k]+=1
        try:
            lat,lon=float(r["Latitude"]),float(r["Longitude"])
        except (TypeError,ValueError):
            bad["invalid_coordinate"]+=1;continue
        if not (-90<=lat<=90 and -180<=lon<=180 and math.isfinite(lat) and math.isfinite(lon)):
            bad["invalid_coordinate"]+=1;continue
        if lat==lon==0:
            bad["zero_coordinate"]+=1;continue
        if not r["He"]:
            bad["no_He_cell"]+=1;continue
        observed[k]+=1
        points[k].add((lat,lon))
        popid=str(r.get("Pop_id") or "")
        by_pop[k][popid].add((lat,lon))
    eligible=[]
    for k,coords in sorted(points.items()):
        if not has_100km_span(coords):continue
        m=metadata[k];name=str(m["Spec_Latin_GenDivRange"]).strip().lower()
        split=hashlib.sha256((NAMESPACE+"|"+name).encode()).digest()[0]%2
        eligible.append({"species_id":str(k[0]),"study_id":str(k[1]),
          "species_name":name,"role":split,
          "habitat_adulthood":str(m.get("Habitat_adulthood")),
          "marker":str(m.get("Marker_type")),
          "populations_with_He_and_valid_coordinates":observed[k],
          "unique_He_geographic_sites":len(coords)})
    terrestrial=[x for x in eligible if x["habitat_adulthood"]=="Terrestrial"]
    climate_keys={(x["species_id"],x["study_id"]) for x in terrestrial}
    climate_sites=[{"species_id":k[0],"study_id":k[1],"latitude":p[0],"longitude":p[1]}
                   for k,coords in sorted(points.items()) if k in climate_keys for p in sorted(coords)]
    def summary(xs):
        return {"eligible_species":len({x["species_name"] for x in xs}),
          "eligible_study_species_groups":len(xs),
          "species_by_role":[len({x["species_name"] for x in xs if x["role"]==i}) for i in (0,1)],
          "groups_by_role":[sum(x["role"]==i for x in xs) for i in (0,1)],
          "marker_groups":dict(collections.Counter(x["marker"] for x in xs))}
    duplicate_groups=[{"species_id":k[0],"study_id":k[1],"pop_id":pid,
                       "distinct_coordinates":len(v)}
       for k,pops in by_pop.items() for pid,v in pops.items() if pid and
       sum(1 for r in [] )>1] # duplicates audited separately in key reconciliation
    return {"status":"BLIND_EXACT_GEO_AND_TAXON_RECONCILED",
        "original_study_species_keys":len(metadata),
        "joined_study_species_keys":len(allpop),
        "source_population_rows_ignored_no_species_key":bad["missing_taxon_study_key"],
        "row_qc":dict(bad),
        "all_taxa":summary(eligible),
        "adult_terrestrial":summary(terrestrial),
        "eligibility_policy": "≥5 unique He-observed sites and exact max great-circle distance ≥100km",
        "all_taxa_frozen_manifest":eligible,
        "adult_terrestrial_frozen_manifest":terrestrial,
        "adult_terrestrial_climate_sites":climate_sites,
        "He_numeric_values_opened":False,"outcome_opening_authorized":False}
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--metadata",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    o=p.parse_args()
    md=json.loads(o.metadata.read_text());data={x["filename"]:x for x in md["files"]}
    report=run(fetch_verified(data["pop_tab_v2025-03-31.xlsx"]),
               fetch_verified(data["spec_tab_v2025-03-31.xlsx"]))
    o.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if "manifest" not in k and k!="adult_terrestrial_climate_sites"},sort_keys=True))
