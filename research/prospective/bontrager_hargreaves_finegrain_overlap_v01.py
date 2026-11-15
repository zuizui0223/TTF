"""Response-blind fine-grain cross-original metadata overlap: Bontrager x Hargreaves.

Two authentic published source datasets must match exact pinned SHA256/MD5
BEFORE any row is decoded. No fitness values, signs, sample response, climate
association or genetic response is indexed, calculated or modeled.

This measures overlap of source-study, species, source-ID, transplant-garden
and fitness-METRIC TYPE labels; it does not equate IDs to independent experiments.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import zipfile
from collections import defaultdict, Counter
from pathlib import Path
from urllib.parse import urlparse, unquote
from urllib.request import urlopen, Request

HERE = Path(__file__).resolve().parent
CONTRACT = HERE / "bontrager_hargreaves_finegrain_overlap_contract_v01.json"
MAX_CSV = 25_000_000

def norm_taxon(value: str) -> str | None:
    candidate = " ".join(value.strip().replace("_", " ").split())
    if re.fullmatch(r"[A-Za-z-]+ [A-Za-z-]+", candidate):
        return candidate.casefold()
    return None

def valid_coord(lon: str, lat: str):
    try:
        x, y = float(lon), float(lat)
    except (ValueError, TypeError, OverflowError):
        return None
    return (x, y) if math.isfinite(x) and math.isfinite(y) and -180 <= x <= 180 and -90 <= y <= 90 else None

def haversine(a,b):
    lo1,la1=a;lo2,la2=b
    dlat=math.radians(la2-la1)
    dlon=math.radians((lo2-lo1+180)%360-180)
    z=math.sin(dlat/2)**2+math.cos(math.radians(la1))*math.cos(math.radians(la2))*math.sin(dlon/2)**2
    return 6371.0088*2*math.asin(math.sqrt(max(0.,min(z,1.))))

def zenodo_file_metadata(entry:dict) -> dict:
    rid=entry["zenodo_record"]
    url=f"https://zenodo.org/api/records/{rid}"
    request=Request(url,headers={"Accept":"application/json",
                    "User-Agent":"TTF-cross-original-source-label-metadata-20261115"})
    with urlopen(request,timeout=35) as response:
        raw=response.read(500_001)
    if len(raw)>500_000:
        raise ValueError("Public preservation record metadata exceeds cap")
    info=json.loads(raw)
    if str(info.get("id"))!=str(rid):
        raise ValueError("Preservation record identifier differs")
    meta=info.get("metadata") or {}
    identifiers=[v.get("identifier") for k in ("related_identifiers","identifiers","alternate_identifiers")
                 for v in (meta.get(k) or []) if isinstance(v,dict)]
    def ndoi(s):
        if not isinstance(s,str):return ""
        return re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:)", "", s.strip().lower()).rstrip("/")
    if entry["doi"] not in {ndoi(v) for v in identifiers}:
        raise ValueError("No exact structured original Dryad DOI for Zenodo source")
    group=info.get("files") or {}
    if isinstance(group,list):
        files=group
    elif isinstance(group,dict) and isinstance(group.get("entries"),dict):
        files=list(group["entries"].values())
    else:
        raise ValueError("Unsupported official Zenodo file listing")
    selected=[f for f in files if isinstance(f,dict) and f.get("key")==entry["file"]]
    if len(selected)!=1:
        raise ValueError("Exact preserved file not uniquely available")
    item=selected[0]
    if item.get("size")!=entry["size"] or str(item.get("checksum","")).lower()!="md5:"+entry["md5"]:
        raise ValueError("Published original source checksum/length mismatch")
    links=item.get("links") or {}
    link=links.get("content") or links.get("download") or links.get("self")
    if not isinstance(link,str):
        raise ValueError("Official original file lacks a usable public link")
    parsed=urlparse(link)
    path=unquote(parsed.path)
    allowed=(path.startswith(f"/api/records/{rid}/files/") and path.endswith("/content")
        or path.startswith("/api/files/") and path.endswith("/"+entry["file"])
        or path.startswith(f"/records/{rid}/files/") and path.endswith("/"+entry["file"]))
    if not(parsed.scheme=="https" and parsed.hostname=="zenodo.org" and parsed.port in (None,443)
           and parsed.username is None and parsed.password is None and allowed):
        raise ValueError("Non-official or unexpected original file URL")
    return link

def original_bytes(entry):
    link=zenodo_file_metadata(entry)
    req=Request(link,headers={"Accept":"*/*","User-Agent":"TTF-cross-original-source-only-20261115"})
    with urlopen(req,timeout=75) as response:
        raw=response.read(entry["size"]+1)
    if (len(raw)!=entry["size"] or hashlib.md5(raw).hexdigest()!=entry["md5"]
        or hashlib.sha256(raw).hexdigest()!=entry["sha256"]):
        raise ValueError("Original source bytes failed published and frozen digests")
    return raw

def read_bontrager_csv(archive_bytes, member):
    with zipfile.ZipFile(io.BytesIO(archive_bytes)) as z:
        if member not in z.namelist():
            raise ValueError("Original Bontrager ZIP member absent")
        m=z.getinfo(member)
        if m.file_size>MAX_CSV:
            raise ValueError("Bontrager source CSV member exceeds cap")
        raw=z.read(member)
    if len(raw)>MAX_CSV:
        raise ValueError("Original source member exceeds allowed CSV size")
    return raw

def metadata_census(raw:bytes, *, source:str, c:dict):
    if source not in ("bontrager","hargreaves"):
        raise ValueError("Unknown original source cohort")
    if source=="hargreaves":
        first=raw.split(b"\n",1)[0]+b"\n"
        if hashlib.sha256(first).hexdigest()!=c["original"]["hargreaves"]["expected_header_sha256"]:
            raise ValueError("Hargreaves first original header SHA drift")
        try:
            decoded=raw.decode("utf-8-sig","strict")
            encoding="utf-8-sig"
        except UnicodeDecodeError:
            decoded=raw.decode("cp1252","strict")
            encoding="cp1252"
    else:
        decoded=raw.decode("utf-8-sig","strict")
        encoding="utf-8-sig"
    reader=csv.reader(io.StringIO(decoded,newline=""),strict=True)
    header=next(reader)
    if len(header)!=c["original"][source]["expected_csv_columns"] or len(set(header))!=47:
        raise ValueError("Original 47-field source schema changed")
    allowed=c["strict_metadata_only_fields"]
    if tuple(allowed)!=("study","functional_group","taxon","source","site","fitness_type","lat_site","long_site"):
        raise ValueError("Metadata allowlist includes a forbidden outcome column")
    if not set(allowed).issubset(header):
        raise ValueError("Required source/garden original metadata fields absent")
    ix={field:header.index(field) for field in allowed}
    kept=set();sites=defaultdict(set);study_all=set()
    fields=c["exact_join_hierarchy"]
    levels={tuple(key):set() for key in fields}
    nrows=nplant=n_binomial=0
    blank_counts=Counter()
    for line in reader:
        if not line:
            continue
        nrows+=1
        if nrows>100_000 or len(line)!=47:
            raise ValueError("Original source row limit or width mismatch")
        m={field:line[pos].strip() for field,pos in ix.items()}
        if m["functional_group"] not in c["original_plant_groups"]:
            continue
        nplant+=1
        if m["study"]:
            study_all.add(m["study"])
        taxon=norm_taxon(m["taxon"])
        if taxon:
            n_binomial+=1
        vals={"study":m["study"],"taxon":taxon,
              "source":m["source"],"site":m["site"],"fitness_type":m["fitness_type"]}
        for key in levels:
            if all(vals.get(f) for f in key):
                levels[key].add(tuple(vals[f] for f in key))
            else:
                blank_counts["/".join(key)]+=1
        if vals["study"] and taxon and vals["site"]:
            sites[(vals["study"],taxon,vals["site"])].add(valid_coord(m["long_site"],m["lat_site"]))
    return {
        "original_rows":nrows,"plant_functional_rows":nplant,
        "binomial_plant_rows":n_binomial, "plant_study_labels":study_all,
        "levels":levels,"garden_positions":sites, "encoding":encoding,
        "missingness_rows_by_level":dict(sorted(blank_counts.items())),
    }

def count_matched_garden_coords(a:dict,b:dict):
    keys=set(a)&set(b)
    counts=Counter()
    for key in keys:
        x=a[key];y=b[key]
        if len(x)!=1 or len(y)!=1 or None in x or None in y:
            counts["missing_or_ambiguous"]+=1
            continue
        distance=haversine(next(iter(x)),next(iter(y)))
        if distance<=1:
            counts["distance_le_1km"]+=1
        elif distance<=10:
            counts["distance_gt_1km_le_10km"]+=1
        else:
            counts["distance_gt_10km"]+=1
    return {"shared_source_study_species_garden_keys":len(keys), **dict(sorted(counts.items()))}

def compare(a:dict,b:dict,reference:dict,c:dict):
    labels_ref={item["original_study_label"] for item in reference["catalog"]}
    ref_sha=hashlib.sha256(json.dumps(reference,ensure_ascii=False,
                     sort_keys=True,indent=2).encode("utf-8")).hexdigest()
    # Original GitHub file has a terminating newline, use byte SHA in run().
    if len(labels_ref)!=c["prior_reference_counts"]["bontrager_qualified_studies"]:
        raise ValueError("Reference study-label cohort size drift")
    if len(b["plant_study_labels"])!=c["prior_reference_counts"]["hargreaves_plant_study_labels"]:
        raise ValueError("Hargreaves source study-label cohort size drift")
    exact_qualified=len(labels_ref & b["plant_study_labels"])
    if exact_qualified!=c["prior_reference_counts"]["exact_qualified_study_label_overlap"]:
        raise ValueError("Prior frozen 102 source label comparison not reproduced")
    both_full=a["plant_study_labels"]&b["plant_study_labels"]
    h_only=b["plant_study_labels"]-a["plant_study_labels"]
    b_only=a["plant_study_labels"]-b["plant_study_labels"]
    report={}
    for level in c["exact_join_hierarchy"]:
        key=tuple(level)
        left=a["levels"][key];right=b["levels"][key]
        intersection=left&right
        report["__".join(level)]={
            "bontrager_all_plant_unique_keys":len(left),
            "hargreaves_plant_unique_keys":len(right),
            "exact_shared_keys":len(intersection),
            "hargreaves_keys_not_shared_with_bontrager":len(right-intersection),
            "bontrager_keys_not_shared_with_hargreaves":len(left-intersection),
        }
    return {
        "schema":"ttf_Bontrager_Hargreaves_two_original_source_finegrain_overlap_v01",
        "status":"SOURCE_ONLY_FINEGRAIN_EXPERIMENT_LABEL_OVERLAP_NOT_REPLICATE_CERTIFICATION",
        "original_source_record_counts":{
            "bontrager_rows":a["original_rows"],
            "hargreaves_rows":b["original_rows"],
            "bontrager_plant_rows":a["plant_functional_rows"],
            "hargreaves_plant_rows":b["plant_functional_rows"],
        },
        "original_text_encodings":{
            "bontrager":a["encoding"],"hargreaves":b["encoding"]
        },
        "qualified_Bontrager_reference_labels":len(labels_ref),
        "qualified_source_study_overlap_reproduced":exact_qualified,
        "all_Bontrager_plant_original_study_labels":len(a["plant_study_labels"]),
        "all_Hargreaves_plant_original_study_labels":len(b["plant_study_labels"]),
        "all_Bontrager_Hargreaves_shared_original_study_labels":len(both_full),
        "Hargreaves_study_labels_not_in_any_Bontrager_original_plant_row":len(h_only),
        "Hargreaves_study_labels_only_unmatched_to_qualified_Bontrager_but_found_in_full_original":len((b["plant_study_labels"]-labels_ref)&a["plant_study_labels"]),
        "Bontrager_plant_study_labels_not_in_Hargreaves":len(b_only),
        "Hargreaves_unlinked_full_study_labels_for_followup":sorted(h_only),
        "nested_source_key_overlaps":report,
        "shared_transplant_garden_coordinate_agreement":count_matched_garden_coords(a["garden_positions"],b["garden_positions"]),
        "missingness_rows_by_source":{"bontrager":a["missingness_rows_by_level"],"hargreaves":b["missingness_rows_by_level"]},
        "all_exact_source_keys_shared_do_NOT_certify_identical_experiments":True,
        "unlinked_study_labels_do_NOT_certify_independent_experiments":True,
        "response_columns_indexed_or_fitted":[],
        "fitness_effect_or_direction_opened":False,
        "raw_population_or_garden_coordinates_exported":False,
        "Bontrager_44_of_50_frozen_gate_unchanged":True,
    }

def run():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    c=json.loads(CONTRACT.read_text())
    if c.get("status")!="FROZEN_BEFORE_TWO_ORIGINAL_SOURCE_DATA_ROW_KEY_AUDIT":
        raise ValueError("No new source-only pre-data authorization")
    refpath=HERE.parents[1]/c["qualified_reference"]
    refbytes=refpath.read_bytes()
    if hashlib.sha256(refbytes).hexdigest()!=c["reference_SHA256"]:
        raise ValueError("Original qualified study reference catalog SHA drift")
    reference=json.loads(refbytes)
    # All eight original response-independent metadata field choices frozen above.
    raw_b=original_bytes(c["original"]["bontrager"])
    raw_h=original_bytes(c["original"]["hargreaves"])
    raw_bcsv=read_bontrager_csv(raw_b,c["original"]["bontrager"]["member"])
    a=metadata_census(raw_bcsv,source="bontrager",c=c)
    b=metadata_census(raw_h,source="hargreaves",c=c)
    out=compare(a,b,reference,c)
    out["original_input_SHA256_verified"]={
        "bontrager":c["original"]["bontrager"]["sha256"],
        "hargreaves":c["original"]["hargreaves"]["sha256"]}
    out["reference_catalog_SHA256_verified"]=c["reference_SHA256"]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    summary={k:out[k] for k in (
        "qualified_source_study_overlap_reproduced",
        "all_Bontrager_plant_original_study_labels",
        "all_Bontrager_Hargreaves_shared_original_study_labels",
        "Hargreaves_study_labels_not_in_any_Bontrager_original_plant_row",
        "Hargreaves_study_labels_only_unmatched_to_qualified_Bontrager_but_found_in_full_original")}
    summary["nested_source_key_overlaps"]={k:v["exact_shared_keys"] for k,v in out["nested_source_key_overlaps"].items()}
    summary["shared_transplant_garden_coordinate_agreement"]=out["shared_transplant_garden_coordinate_agreement"]
    print("TTF_TWO_ORIGINAL_SOURCE_FINEGRAIN_METADATA=",json.dumps(summary,sort_keys=True))

if __name__=="__main__":
    run()
