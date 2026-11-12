"""PROCEED v6.2 response-blind metadata/coordinate census (no rate inspection).

Only the listed source *metadata* columns are decoded to analysis variables.
No response column is referenced for any species, location, design or
eligibility decision. The full original source bytes remain private to an
ephemeral CI runner and are not an output artifact.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,urllib.request
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import urlparse
from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed_v62_metadata_row_census_contract_v01.json"
MAX_BYTES=20_000_000
CATEGORY_FIELDS=("system","phylum","experimental","exp_type","genphen",
                 "design","trait_type","data_type","climate_change")


def fetch_original_datafile_bytes(fileid, maximum=MAX_BYTES):
    if not isinstance(fileid,int) or fileid <= 0: raise ValueError("Invalid frozen file ID")
    url=f"https://borealisdata.ca/api/access/datafile/{fileid}"
    req=urllib.request.Request(url,headers={
        "Accept":"text/tab-separated-values",
        "User-Agent":"TTF-PROCEED-v62-response-blind-metadata-census"})
    h1=hashlib.sha256();h2=hashlib.md5(); chunks=[];count=0
    with urllib.request.urlopen(req,timeout=100) as resp:
        while True:
            piece=resp.read(1024*1024)
            if not piece:break
            count+=len(piece)
            if count>maximum:raise ValueError("Source file larger than source-only byte limit")
            h1.update(piece);h2.update(piece);chunks.append(piece)
    return b"".join(chunks),h2.hexdigest(),h1.hexdigest()


def valid_binomial(raw):
    value=" ".join(str(raw or "").replace("_"," ").split())
    if not re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+(?: [a-z][A-Za-z-]+)?",value):
        return None
    if any(token in {"sp","spp","cf","aff"} for token in value.lower().split()):
        return None
    return value


def valid_coordinates(lat,lon):
    try:
        y=float(lat);x=float(lon)
        return (y,x) if math.isfinite(y) and math.isfinite(x) and -90<=y<=90 and -180<=x<=180 else None
    except (ValueError,TypeError,OverflowError):
        return None


def safe_counts(counter,max_labels=60):
    if len(counter)>max_labels:
        return {"distinct_labels":len(counter),"status":"MANY_DISTINCT_LABELS_NOT_ASSESSED"}
    return {str(k or "<missing>"):int(v) for k,v in sorted(counter.items())}


def connected_components(study_to_species):
    """Bipartite response-blind graph: repeated studies cannot leak across species split."""
    sp_to_studies=defaultdict(set)
    for study,species in study_to_species.items():
        for s in species:sp_to_studies[s].add(study)
    unseen=set(sp_to_studies)
    components=[]
    while unseen:
        seed=min(unseen); todo=[seed]; seen_sp=set();seen_st=set()
        while todo:
            sp=todo.pop()
            if sp in seen_sp:continue
            seen_sp.add(sp);unseen.discard(sp)
            for study in sp_to_studies[sp]:
                if study in seen_st:continue
                seen_st.add(study)
                todo.extend(study_to_species[study]-seen_sp)
        components.append({"species_count":len(seen_sp),"study_count":len(seen_st)})
    return sorted(components,key=lambda x:(-x["species_count"],-x["study_count"]))


def metadata_only_summary(raw,contract):
    import inspect
    allowed=contract["source_metadata_columns_authorized"]
    prohibited=set(contract["outcome_fields_must_not_be_decoded"])
    if prohibited.intersection(allowed):
        raise ValueError("Response column leaked into allowlist")
    source=io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",newline="")
    reader=csv.reader(source,delimiter="\t",strict=True)
    header=next(reader,[])
    if len(header)<50 or len(header)!=len(set(header)):
        raise ValueError("Unexpected missing or duplicate original file headers")
    if not set(allowed).issubset(header):
        raise ValueError("Required metadata columns missing from exact original source")
    names={k:header.index(k) for k in allowed}
    all_count=0;nometa=0;noid=0;original_binomial=set();withbinomial=set()
    locations={"both":0,"only_first":0,"only_second":0,"neither":0}
    category={n:Counter() for n in CATEGORY_FIELDS}
    terrestrial_sys=Counter()
    species_studies=defaultdict(set)
    species_both_studies=defaultdict(set)
    study_species_both=defaultdict(set)
    studies_any=set()
    exact_location_pairs=set()
    sample1_coords=set()
    sample2_coords=set()
    metadata_presence=Counter()
    duplicate_es=Counter()
    original_genera=set()
    bad_rows=0
    phylum_map={}
    for values in reader:
        if not values or all(not v.strip() for v in values):continue
        all_count+=1
        if len(values)!=len(header):
            bad_rows+=1
            continue
        # Index ONLY metadata names. No response effect / n / mean / sd ever
        # appears as an index in this function.
        row={key:values[index].strip() for key,index in names.items()}
        study=row["ref_ID"];effect=row["es_ID"]
        if not study or not effect:
            noid+=1
        else:
            duplicate_es[effect]+=1
        for term in CATEGORY_FIELDS:
            category[term][row[term]]+=1
        if row["years"]:metadata_presence["years"]+=1
        if row["generations"]:metadata_presence["generations"]+=1
        if row["gtime"]:metadata_presence["gtime"]+=1
        if row["ref_ID"]:metadata_presence["ref_ID"]+=1
        if row["genphen"]:metadata_presence["genphen"]+=1
        if row["experimental"]:metadata_presence["experimental"]+=1
        if row["exp_type"]:metadata_presence["exp_type"]+=1
        if row["design"]:metadata_presence["design"]+=1
        if row["climate_change"]:metadata_presence["climate_change"]+=1
        species=valid_binomial(row["sp_pub"]) or valid_binomial(row["sp_ncbi"])
        if species:
            withbinomial.add(species)
            original_genera.add(species.split()[0])
            if study: species_studies[species].add(study)
            phylum_map[species]=row["phylum"]
        loc1=valid_coordinates(row["sample1_latitude"],row["sample1_longitude"])
        loc2=valid_coordinates(row["sample2_latitude"],row["sample2_longitude"])
        locations["both" if loc1 and loc2 else "only_first" if loc1 else "only_second" if loc2 else "neither"]+=1
        if loc1:sample1_coords.add(loc1)
        if loc2:sample2_coords.add(loc2)
        if loc1 and loc2:
            exact_location_pairs.add((loc1,loc2))
            if species and study:
                species_both_studies[species].add(study)
                study_species_both[study].add(species)
        if study:studies_any.add(study)

    qualified=sorted(species_both_studies)
    roles={"train":0,"heldout":0}
    for species in qualified:
        d=hashlib.sha256(("ttf-proceed-v62-geography-source-only-v01|"+species.lower()).encode()).digest()
        roles["train" if d[0]&1 else "heldout"]+=1
    cross_studies=Counter()
    for study,sp_list in study_species_both.items():
        rs=set()
        for species in sp_list:
            d=hashlib.sha256(("ttf-proceed-v62-geography-source-only-v01|"+species.lower()).encode()).digest()
            rs.add("train" if d[0]&1 else "heldout")
        cross_studies[len(rs)]+=1
    components=connected_components(study_species_both)
    gate=contract["prospective_gate"]
    # NO biological-response opening authorized here, even if geometry passes.
    geometry_ok=(len(qualified)>=gate["minimum_terrestrial_georeferenced_species"]
                 and len(study_species_both)>=gate["minimum_source_studies"]
                 and len({s.split()[0] for s in qualified})>=gate["minimum_independent_genera"]
                 and min(roles.values())>=gate["minimum_species_per_hash_heldout_arm"])
    # Never mark terrestrial when source system labels have not been interpreted.
    return {
        "schema":"ttf_proceed_v62_source_metadata_row_geometry_census_v01",
        "status":"UNCLASSIFIED_SYSTEM_GEOGRAPHY_POSSIBLE_CODEBOOK_REQUIRED" if geometry_ok else "NOT_EVALUABLE_GEOGRAPHIC_SUPPORT_WITH_CURRENT_QUALIFICATION",
        "all_source_rows":all_count,"malformed_rows":bad_rows,
        "effect_IDs_missing_or_study_ID_missing":noid,
        "duplicate_effect_id_keys":sum(n>1 for n in duplicate_es.values()),
        "total_published_binomial_labels":len(withbinomial),
        "species_with_both_georeferenced_original_endpoints":len(qualified),
        "distinct_genus_with_both_endpoints":len({x.split()[0] for x in qualified}),
        "study_ids_with_both_georeferenced_endpoints":len(study_species_both),
        "source_study_ids_any_record":len(studies_any),
        "species_hash_roles_both_coords":roles,
        "study_label_cross_species_arm":cross_studies.get(2,0),
        "source_study_species_graph_components":components,
        "record_geolocation_status":locations,
        "distinct_original_coordinate_pairs":len(exact_location_pairs),
        "distinct_first_origin_coordinates":len(sample1_coords),
        "distinct_second_origin_coordinates":len(sample2_coords),
        "metadata_presence_counts":dict(metadata_presence),
        "original_metadata_categories":{term:safe_counts(ct) for term,ct in category.items()},
        "prequalified_original_terrestrial_status":"PENDING_SOURCE_SYSTEM_CODEBOOK",
        "genphen_or_experiment_category_heritable_interpretation":"NOT_ASSIGNED_UNTIL_OFFICIAL_README",
        "geographic_qualification_not_inference_authorization":True,
        "rate_or_phenotype_values_summarized":False,
        "genetic_values_opened":False,
        "old_TTF_frozen_outcomes_reopened":False,
        "source_sites_or_effects_publicly_exported":False
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    c=json.loads(CONTRACT.read_text())
    if c.get("status")!="PRE_RESPONSE_SOURCE_CENSUS_FROZEN_AFTER_HEADER_5_OF_5_PASS":
        raise ValueError("Source-only protocol not frozen")
    f=c["header_authority"]["rate_source"]
    e={
        "schema":"ttf_proceed_v62_source_metadata_row_geometry_census_v01",
        "source_file":f["name"],
        "rate_or_phenotype_values_summarized":False,
        "genetic_values_opened":False,
        "old_TTF_frozen_outcomes_reopened":False
    }
    try:
        src=json.loads((HERE/"proceed_v62_source_header_contract_v01.json").read_text())
        official=src["current_source"]
        catalog=validated_get_catalog("https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(catalog,official)
        raw,md5,sha256=fetch_original_datafile_bytes(f["id"])
        e["observed_source_size_bytes"]=len(raw)
        e["observed_source_md5"]=md5
        e["observed_source_sha256"]=sha256
        if len(raw)!=f["expected_bytes"] or md5!=f["expected_md5"]:
            raise ValueError("Original source rate-table size/MD5 differs from frozen Borealis original")
        e.update(metadata_only_summary(raw,c))
        e["source_file"]=f["name"]
        e["source_md5_verified"]=True
        e["observed_source_sha256"]=sha256
    except Exception as exc:
        e["status"]="SOURCE_RETRIEVAL_OR_METADATA_CENSUS_STOP_NOT_BIOLOGICAL_NULL"
        e["error_type"]=type(exc).__name__
        e["error_message"]=str(exc)[:350]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(e,sort_keys=True,indent=2)+"\n")
    summary={k:e.get(k) for k in ("status","all_source_rows","record_geolocation_status",
              "species_with_both_georeferenced_original_endpoints",
              "study_ids_with_both_georeferenced_endpoints","source_md5_verified")}
    print(json.dumps(summary,sort_keys=True))


if __name__=="__main__":main()
