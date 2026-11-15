"""Bontrager original transplant SOURCE METADATA qualification, no fitness value read.

Original Dryad SHA256 = mirrored Borealis source ZIP MD5 required.
The CSV reader necessarily tokenizes whole lines, but the forbidden
'fitness' column is NEVER indexed, decoded numerically, filtered, summarized
or used for matching. Only allowed study, taxon, source/garden identifiers,
lat/lon and original response-TYPE categories are touched. Output aggregate.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen

from bontrager_mirror_zip_v01 import verify_original_zip,ORIGINAL_URL

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"bontrager_transplant_geo_contract_v01.json"
SOURCE=HERE/"bontrager_mirror_zip_contract_v01.json"


def binomial(name):
    n=" ".join(str(name or "").strip().replace("_"," ").split())
    return n if re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+",n) else None


def original_lonlat(lon,lat):
    try:
        x=float(lon);y=float(lat)
    except (ValueError,TypeError,OverflowError):
        return None
    if not (math.isfinite(x) and math.isfinite(y) and -180<=x<=180 and -90<=y<=90):
        return None
    return (x,y)


def connected_components(study_species):
    byspecies=defaultdict(set)
    for study,ss in study_species.items():
        for s in ss:byspecies[s].add(study)
    unseen=set(byspecies);blocks=[]
    while unseen:
        todo=[min(unseen)];taxa=set();papers=set()
        while todo:
            s=todo.pop()
            if s in taxa:continue
            taxa.add(s);unseen.discard(s)
            for k in byspecies[s]:
                if k not in papers:
                    papers.add(k);todo.extend(study_species[k]-taxa)
        blocks.append((taxa,papers))
    return blocks


def component_disjoint_split(study_species):
    blocks=connected_components(study_species)
    blocks.sort(key=lambda x:(-len(x[0]),
                hashlib.sha256("|".join(sorted(x[0])).encode()).hexdigest()))
    stax={"train":set(),"heldout":set()}
    sstudy={"train":set(),"heldout":set()}
    for taxa,papers in blocks:
        role="train" if len(stax["train"])<=len(stax["heldout"]) else "heldout"
        stax[role].update(taxa);sstudy[role].update(papers)
    if stax["train"]&stax["heldout"] or sstudy["train"]&sstudy["heldout"]:
        raise AssertionError("Study or taxon leaked into both source-only roles")
    return {
        "connected_study_species_components":len(blocks),
        "largest_component_species":max((len(c[0]) for c in blocks),default=0),
        "largest_component_studies":max((len(c[1]) for c in blocks),default=0),
        "species_per_disjoint_role":{k:len(v) for k,v in stax.items()},
        "original_studies_per_disjoint_role":{k:len(v) for k,v in sstudy.items()},
        "species_and_studies_disjoint":True,
        "source_only_membership_sha256":hashlib.sha256(json.dumps(
           {k:sorted(v) for k,v in stax.items()},sort_keys=True,
           separators=(",",":")).encode()).hexdigest()
    }


def metadata_census(original_csv,contract):
    allowed=contract["allowed_metadata_only_columns"]
    if set(allowed)&set(contract["absolute_forbidden_response_columns"]):
        raise ValueError("Biological response in source metadata allowlist")
    reader=csv.reader(io.TextIOWrapper(io.BytesIO(original_csv),
                                       encoding="utf-8-sig",newline=""),
                      delimiter=",",strict=True)
    header=next(reader,[])
    if len(header)!=47 or len(set(header))!=47 or not set(allowed).issubset(header):
        raise ValueError("Original 47-field source header changed")
    at={name:header.index(name) for name in allowed}
    rows=0;source_types=Counter();reason=Counter()
    study_species=defaultdict(set)
    source_ids=defaultdict(set);garden_ids=defaultdict(set)
    candidates=[]
    categories=Counter()
    for cells in reader:
        if not cells:continue
        rows+=1
        if len(cells)!=len(header):
            raise ValueError("Original CSV row width changed")
        # Index *only* declared source metadata. Never inspect cells[index_of_fitness].
        m={name:cells[i].strip() for name,i in at.items()}
        source_types[m["fitness_type"] or "<missing>"]+=1
        categories[m["functional_group"] or "<missing>"]+=1
        study=m["study"];taxon=binomial(m["taxon"])
        if not(study and taxon and m["source"] and m["site"] and m["fitness_type"]):
            reason["missing_study_or_binomial_or_site_or_type"]+=1
            continue
        sc=original_lonlat(m["long_source"],m["lat_source"])
        gc=original_lonlat(m["long_site"],m["lat_site"])
        if sc is None or gc is None:
            reason["missing_or_invalid_original_coordinate"]+=1
            continue
        reason["nominal_eligible_original_locations"]+=1
        # Study×species labels are not independent data by themselves.
        sourcekey=(study,taxon,m["source"])
        gardenkey=(study,taxon,m["site"])
        source_ids[sourcekey].add(sc)
        garden_ids[gardenkey].add(gc)
        candidates.append((study,taxon,m["source"],m["site"],
                           m["fitness_type"],sc,gc))
    conflict_source={k for k,z in source_ids.items() if len(z)>1}
    conflict_garden={k for k,z in garden_ids.items() if len(z)>1}
    experiment_cells=defaultdict(set)
    cell_source_coords=defaultdict(set)
    unique_fivepart_rows=set()
    unique_geo_pairs=set()
    for study,taxon,source,site,ft,sc,gc in candidates:
        if ((study,taxon,source) in conflict_source or
            (study,taxon,site) in conflict_garden):
            reason["excluded_ambiguous_source_or_garden_ID"]+=1
            continue
        cell=(study,taxon,site,ft)
        experiment_cells[cell].add(source)
        cell_source_coords[cell].add(sc)
        unique_fivepart_rows.add((study,taxon,source,site,ft))
        unique_geo_pairs.add((study,taxon,source,site,ft,sc,gc))
        study_species[study].add(taxon)
    qualified={
        k for k,v in experiment_cells.items()
        if len(v)>=2 and len(cell_source_coords[k])>=2
    }
    qualifying_studies={k[0] for k in qualified}
    qualifying_taxa={k[1] for k in qualified}
    qualifying_genera={s.split()[0] for s in qualifying_taxa}
    qualifying_graph=defaultdict(set)
    qualifying_sources=set()
    qualifying_categories=Counter()
    for k in qualified:
        qualifying_graph[k[0]].add(k[1])
        qualifying_sources.update((k[0],k[1],s) for s in experiment_cells[k])
        qualifying_categories[k[3]]+=1
    roles=component_disjoint_split(qualifying_graph)
    pre=contract["future_macro_gate"]
    possible=(len(qualifying_taxa)>=pre["minimum_binomial_species"]
        and len(qualifying_studies)>=pre["minimum_original_studies"]
        and len(qualifying_genera)>=pre["minimum_genera"]
        and min(roles["species_per_disjoint_role"].values())>=min(
              pre["minimum_species_train"],pre["minimum_species_heldout"]))
    return {
      "schema":"ttf_Bontrager_new_transplant_original_metadata_geographic_contrast_v01",
      "status":("GEOGRAPHY_MACRO_SOURCE_POSSIBLE_REQUIRES_HISTORY_CLIMATE_AND_STUDY_OVERLAP"
                if possible else pre["if_fail"]),
      "original_source_rows":rows,
      "original_source_column_count":len(header),
      "source_metadata_exclusion_counts":dict(sorted(reason.items())),
      "raw_fitness_type_metadata_count":dict(sorted(source_types.items())),
      "functional_group_metadata_count":dict(sorted(categories.items())),
      "source_ID_with_multiple_coordinates":len(conflict_source),
      "garden_ID_with_multiple_coordinates":len(conflict_garden),
      "original_study_species_source_ID_keys_with_coordinates":len(source_ids),
      "original_study_species_garden_ID_keys_with_coordinates":len(garden_ids),
      "unique_study_species_source_site_fitness_type_rows":len(unique_fivepart_rows),
      "unique_study_species_source_site_fitness_type_geo_rows":len(unique_geo_pairs),
      "original_site_same_study_taxon_fitness_type_cells":len(experiment_cells),
      "shared_garden_fitness_type_cells_with_multiple_source_IDs_and_coordinates":len(qualified),
      "qualifying_cells_by_fitness_type":dict(sorted(qualifying_categories.items())),
      "qualified_original_source_studies":len(qualifying_studies),
      "qualified_binomial_taxa":len(qualifying_taxa),
      "qualified_genera":len(qualifying_genera),
      "qualified_distinct_study_species_source_populations":len(qualifying_sources),
      "study_species_component_split":roles,
      "origin_and_garden_WGS84_coordinates_present_not_climatically_matched":True,
      "source_studies_independent_of_Hargreaves_NOT_verified":True,
      "fitness_measurement_sign_or_magnitude_opened":False,
      "adaptive_fitness_or_speciation_prediction_performed":False,
      "raw_source_and_garden_coordinates_exported":False,
      "unresolved_published_local_adaptation_meta_selection_bias":True
    }


def run():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    contract=json.loads(CONTRACT.read_text())
    lock=json.loads(SOURCE.read_text())
    if contract["status"]!="FROZEN_BEFORE_ORIGINAL_SOURCE_ROW_METADATA_CENSUS":
        raise ValueError("No original cohort source-only contract")
    receipt={"schema":"ttf_Bontrager_new_transplant_original_metadata_geographic_contrast_v01",
             "fitness_measurement_sign_or_magnitude_opened":False,
             "adaptive_fitness_or_speciation_prediction_performed":False,
             "raw_source_and_garden_coordinates_exported":False}
    try:
        req=Request(ORIGINAL_URL,headers={"User-Agent":"TTF-Bontrager-transplant-origin-garden-source-census-v01"})
        with urlopen(req,timeout=90) as response:
            original=response.read(lock["max_download_bytes"]+1)
        verify_original_zip(original,lock)
        with zipfile.ZipFile(io.BytesIO(original)) as z:
            with z.open(contract["source_authority"]["source_member"]) as f:
                # Keep original source bytes ephemeral, parse ONLY allowlisted metadata cells.
                raw=f.read(25_000_000)
        if len(raw)>=25_000_000:
            raise ValueError("CSV source exceeds bounded source-only cohort size")
        receipt.update(metadata_census(raw,contract))
        receipt["original_archive_Dryad_sha256"]=lock["dryad_sha256"]
        receipt["original_archive_Borealis_md5"]=lock["mirror_md5"]
    except Exception as exc:
        receipt.update(status="SOURCE_OR_GEOGRAPHIC_METADATA_STOP_NOT_BIOLOGICAL_NULL",
                       error_type=type(exc).__name__,
                       error_message=str(exc)[:340])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("BONTRAGER_ORIGINAL_META_CENSUS",receipt["status"],receipt.get("qualified_binomial_taxa"),
          receipt.get("qualified_original_source_studies"),
          receipt.get("shared_garden_fitness_type_cells_with_multiple_source_IDs_and_coordinates"))
if __name__=="__main__":run()
