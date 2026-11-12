"""Read only published GIS codebook attributes to predeclare barriers.
Does not sample any genetic break, any surveyed background, or their coordinate.
"""
from __future__ import annotations
import argparse,collections,hashlib,io,json,re,urllib.request,zipfile
from pathlib import Path
import shapefile

TARGETS={"ecoregion_level1","ecoregion_level2","naturalearth_rivers","naturalearth_lakes"}
NAMEWORDS=re.compile(r"superior|huron|michigan|ontario|erie|great\s+bear|great\s+slave|winnipeg|st\.?\s*lawrence",re.I)
def verified_zip(item,sha):
    total=0;h=hashlib.sha256()
    r=urllib.request.Request(item["url"],headers={"User-Agent":"TTF-geography-source-codebook-v01"})
    data=io.BytesIO()
    with urllib.request.urlopen(r,timeout=120) as stream:
        while True:
            block=stream.read(1024*1024)
            if not block:break
            total+=len(block)
            if total>sha["bytes"]:raise ValueError("original source size mismatch")
            data.write(block);h.update(block)
    if total!=sha["bytes"] or h.hexdigest()!=sha["sha256"]:
        raise ValueError("original source SHA mismatch")
    data.seek(0);return zipfile.ZipFile(data)

def attributes(z):
    members=z.namelist()
    stem=next(n[:-4] for n in members if n.lower().endswith(".shp"))
    args={s:io.BytesIO(z.read(stem+"."+s)) for s in ("shp","dbf","shx")}
    reader=shapefile.Reader(**args)
    return [x.as_dict() for x in reader.iterRecords()]

def summarize(key,rows):
    if key.startswith("ecoregion"):
        field="NA_L1CODE" if key=="ecoregion_level1" else "NA_L2CODE"
        names="NA_L1NAME" if key=="ecoregion_level1" else "NA_L2NAME"
        combos=sorted({(str(x.get(field)),str(x.get(names))) for x in rows})
        return {"role":key,"source_polygon_records":len(rows),
                "distinct_region_codes":len({x[0] for x in combos}),
                "geographical_region_names":combos[:80],
                "no_genetic_response_opened":True}
    matches=[]
    for x in rows:
        name=str(x.get("name") or x.get("name_en") or "")
        if NAMEWORDS.search(name):
            matches.append({"name":name,"scalerank":x.get("scalerank"),"featurecla":x.get("featurecla")})
    return {"role":key,"source_feature_records":len(rows),
            "rank_distribution":dict(collections.Counter(str(x.get("scalerank")) for x in rows)),
            "potential_major_water_name_matches":matches[:80],
            "no_species_coordinates_read":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--sources",type=Path,required=True)
    p.add_argument("--sha-lock",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    src=json.loads(a.sources.read_text())
    verified={x["term"]:x for x in json.loads(a.sha_lock.read_text())["asset_hashes"]}
    results=[]
    for item in src["source_candidates"]:
        term=item["term"]
        if term not in TARGETS:continue
        with verified_zip(item,verified[term]) as z:
            summary=summarize(term,attributes(z))
        results.append(summary)
        print("ORIGINAL GIS CODEBOOK",term,summary.get("distinct_region_codes"),flush=True)
    result={"schema":"ttf_jensen_barrier_GIS_codebook_v01",
        "status":"PUBLIC_GIS_ATTRIBUTES_VERIFIED_PRE_BREAK_MODEL",
        "sources":results,"actual_break_outcomes_or_backgrounds_read":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,default=str)+"\n")
    print(json.dumps({"status":result["status"],"terms":[x["role"] for x in results]}))
if __name__=="__main__":main()
