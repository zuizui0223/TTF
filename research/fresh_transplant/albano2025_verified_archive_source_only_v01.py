"""Albano 2025: SHA-pinned original Dryad ZIP source-only audit.

No fitness, genotype, growth, seed or LGM response fields are selected.
The CSV reader necessarily scans rows but only the named provenance
columns are accessed, computed from or exported.
Keep the individual source-coordinate registry LOCAL, not in Git.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,zipfile
from pathlib import Path

ORIGINAL_CSV="Albano_et_al_2025_EM_Experiment_Data.csv"
ORIGINAL_BYTES=421805
ORIGINAL_SHA256="b3796db37145dd32ab148f5c2105ca0646923feab153d8857d4f90a5539f080c"
SOURCE_FIELDS=("Garden","Position","Population","Rep","ContOrigin",
               "ContGarden","Latitude","Longitude")
EXPECTED_GARDENS={"Lafayette","Mississauga","Montpellier","Uppsala"}
EXPECTED_ORIGINS={"NAM","EUR"}
R_EARTH_KM=6371.0088

def haversine(a,b):
    lat1,lon1=map(math.radians,a)
    lat2,lon2=map(math.radians,b)
    dlat=lat2-lat1
    dlon=lon2-lon1
    h=math.sin(dlat/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return 2*R_EARTH_KM*math.asin(min(1.0,math.sqrt(max(0.0,h))))

def audit_zip(path:Path, verify_pinned:bool=True):
    with zipfile.ZipFile(path) as archive:
        entries=[entry for entry in archive.infolist()
                 if entry.filename==ORIGINAL_CSV and not entry.is_dir()]
        if len(entries)!=1:raise ValueError("Exact original CSV missing or duplicate")
        entry=entries[0]
        if entry.file_size>2_000_000:raise ValueError("Original CSV inflated-size cap")
        with archive.open(entry) as fh:raw=fh.read(2_000_001)
    digest=hashlib.sha256(raw).hexdigest()
    if verify_pinned and (len(raw)!=ORIGINAL_BYTES or digest!=ORIGINAL_SHA256):
        raise ValueError("Original Dryad SHA256 or byte size mismatch")
    reader=csv.reader(io.StringIO(raw.decode("utf-8-sig"),newline=""))
    try:header=next(reader)
    except StopIteration:raise ValueError("Original CSV empty")
    if len(header)!=len(set(header)) or not set(SOURCE_FIELDS).issubset(header):
        raise ValueError("Original provenance field names missing/duplicated")
    index={field:header.index(field) for field in SOURCE_FIELDS}
    rows=[]
    for lineno,source in enumerate(reader,2):
        if not source:continue
        if len(source)!=len(header):raise ValueError(f"CSV width mismatch at line {lineno}")
        # Only the source/geography/identity columns are selected.
        record={field:source[index[field]].strip() for field in SOURCE_FIELDS}
        if any(not v for v in record.values()):raise ValueError("Missing source identity/coordinate field")
        try:
            lat=float(record["Latitude"]);lon=float(record["Longitude"])
        except ValueError as exc:raise ValueError("Non-numeric source coordinates") from exc
        if not(math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180):
            raise ValueError("Source coordinate out of bounds")
        record["latitude"]=lat;record["longitude"]=lon
        rows.append(record)
    if not rows:raise ValueError("No provenance rows")
    gardens={r["Garden"] for r in rows}
    origins={r["ContOrigin"] for r in rows}
    if verify_pinned and (gardens!=EXPECTED_GARDENS or origins!=EXPECTED_ORIGINS):
        raise ValueError("Original garden/origin strata changed")
    individual_keys=collections.Counter((r["Garden"],r["Position"]) for r in rows)
    if any(n!=1 for n in individual_keys.values()):
        raise ValueError("Garden+Position not a unique original individual key")
    rep_keys=collections.Counter((r["Garden"],r["Population"],r["Rep"]) for r in rows)
    by_population=collections.defaultdict(list)
    for r in rows:by_population[r["Population"]].append(r)
    site_rows=[]
    for population,observed in by_population.items():
        loc={(r["latitude"],r["longitude"]) for r in observed}
        origin={r["ContOrigin"] for r in observed}
        if len(loc)!=1 or len(origin)!=1:raise ValueError("Population source location/origin drift")
        lat,lon=next(iter(loc))
        site_rows.append({"population_id_verbatim":population,
                          "latitude":lat,"longitude":lon,
                          "origin":next(iter(origin)),
                          "gardens_observed":sorted({r["Garden"] for r in observed}),
                          "original_csv_individual_rows":len(observed)})
    # Different source IDs sharing a coordinate are not independent sites.
    unique_coordinates={(s["latitude"],s["longitude"]) for s in site_rows}
    locations=[(s["latitude"],s["longitude"]) for s in site_rows]
    near_distances=[haversine(locations[i],locations[j])
                    for i in range(len(locations)) for j in range(i+1,len(locations))
                    if site_rows[i]["origin"]==site_rows[j]["origin"]]
    site_rows.sort(key=lambda s:(s["origin"],int(s["population_id_verbatim"])
                                 if s["population_id_verbatim"].isdigit()
                                 else s["population_id_verbatim"]))
    coverage=collections.Counter(len(s["gardens_observed"]) for s in site_rows)
    observed_pairs=sum(coverage[k]*k for k in coverage)
    possible_pairs=len(site_rows)*len(gardens)
    row_by_garden=collections.Counter(r["Garden"] for r in rows)
    origin_by_pop=collections.Counter(s["origin"] for s in site_rows)
    site_pair_within_1km=sum(d<1.0 for d in near_distances)
    duplicate_rep_groups=sum(n>1 for n in rep_keys.values())
    excess_rep_rows=sum(n-1 for n in rep_keys.values() if n>1)
    receipt={
        "schema":"ttf_Albano2025_SHA_pinned_original_source_only_audit_v01",
        "repository_doi":"10.5061/dryad.bk3j9kdqk",
        "official_Dryad_version_id":394175,
        "official_original_file_id":4337649,
        "original_csv_path":ORIGINAL_CSV,
        "original_csv_bytes":len(raw),"sha256_computed_from_uploaded_original_CSV":digest,
        "official_repository_sha256_MATCHED":verify_pinned,
        "original_csv_columns_total":len(header),
        "rows_original_experiment":len(rows),
        "original_source_population_IDs":len(site_rows),
        "distinct_exact_coordinate_pairs":len(unique_coordinates),
        "source_population_coordinate_conflicts":0,
        "source_population_continent_conflicts":0,
        "source_geographic_provenance_columns_verified":list(SOURCE_FIELDS),
        "source_populations_by_continent":dict(sorted(origin_by_pop.items())),
        "gardens":sorted(gardens),
        "individual_rows_by_garden":dict(sorted(row_by_garden.items())),
        "population_garden_pairs_observed":observed_pairs,
        "population_garden_pairs_full_cross_possible":possible_pairs,
        "population_garden_pairs_absent_from_original_rows":possible_pairs-observed_pairs,
        "population_garden_coverage_frequency":{str(k):coverage[k] for k in sorted(coverage)},
        "individual_Garden_Position_unique":True,
        "Garden_Population_Rep_not_unique":duplicate_rep_groups>0,
        "Garden_Population_Rep_colliding_key_groups":duplicate_rep_groups,
        "Garden_Population_Rep_excess_rows_on_colliding_keys":excess_rep_rows,
        "unique_site_coordinates_nearest_pair_km":round(min(near_distances),6),
        "same_continent_site_pairs_under_1km":site_pair_within_1km,
        "actual_independently_surveyed_source_GPS_provenance_verified":False,
        "source_id_to_reported_coordinate_mapping_verified_in_genuine_file":verify_pinned,
        "all_96_sources_occupy_independent_21ka_climate_cells_proven":False,
        "four_gardens_full_population_factorial":possible_pairs==observed_pairs,
        "phenotype_fitness_genotype_or_LGM_values_modeled":False,
        "historical_21ka_exposure_obtained":False,
        "independent_LGM_fitness_effect_tested":False,
        "old_bontrager_44of50_stop_unchanged":True,
        "status":"ORIGINAL_SOURCE_ID_COORDINATES_AUTHENTICATED_GPS_FIELD_ORIGIN_UNCERTIFIED_LGM_NOT_TESTED",
    }
    detailed={
        "schema":"ttf_Albano2025_LOCAL_SOURCE_SITE_REGISTRY_v01",
        "source_csv_sha256":digest,
        "note":"Data rows reproduce only source IDs, original-source coordinates, origin and garden presence. Do not commit to TTF repository; do not infer field-survey GPS or 21ka histories.",
        "site_count":len(site_rows),
        "sites":site_rows,
    }
    return receipt,detailed

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--local-site-registry",type=Path)
    a=p.parse_args()
    receipt,detail=audit_zip(a.zip)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    if a.local_site_registry:
        a.local_site_registry.parent.mkdir(parents=True,exist_ok=True)
        a.local_site_registry.write_text(json.dumps(detail,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(receipt,sort_keys=True))

if __name__=="__main__":main()
