"""Response-blind geospatial eligibility census for Barley 2021 acclimation data.

This reads ONLY allowed source metadata columns into analysis variables.
The published phenotype ARR/thermal-limit responses are neither inspected,
summarized nor used for source qualification. It does not test a biological
history/plasticity association; no TTF frozen genetic outcome is accessed.
"""
from __future__ import annotations
import argparse
from collections import defaultdict, Counter
import csv
import hashlib
import io
import json
import math
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CONTRACT = ROOT / "quaternary_phenotype_source_screen_v01.json"
AMENDMENT = ROOT / "barley_pretreatment_taxonomy_amendment_v02.json"


def git_blob_sha1(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()


def source_census(csv_bytes: bytes, rules: dict) -> dict:
    """Reject ambiguous geometry. No ARR or thermal-limit column is accessed."""
    amendment = json.loads(AMENDMENT.read_text())
    if amendment["status"] != "FROZEN_PRE_ARR_AMENDMENT_FOR_METADATA_DISCOVERY_ERRORS":
        raise ValueError("No authorized source-only taxonomy amendment")
    animal_phyla = set(amendment["allowlisted_animal_phyla"])
    text = csv_bytes.decode("utf-8-sig")
    reader = csv.reader(io.StringIO(text, newline=""))
    cols = next(reader, [])
    permitted = rules["prephenotype_allowed_fields"]
    prohibited = rules["response_fields_prohibited_until_separate_authorization"]
    if len(set(cols)) != len(cols) or not set(permitted).issubset(cols):
        raise ValueError("Original source header missing or ambiguous")
    if any(field in permitted for field in prohibited):
        raise ValueError("Prohibited outcome field in the geometry allowlist")
    ix = {name: cols.index(name) for name in permitted}
    all_records = 0
    terrestrial_rows = 0
    terrestrial_species = set()
    nonterrestrial = Counter()
    terrestrial_nonanimal_rows = Counter()
    terrestrial_nonanimal_taxa = set()
    terrestrial_animal_phylum_rows = Counter()
    missing_coords = 0
    invalid_coords = 0
    missing_study_pop = 0
    per_population_coords = defaultdict(set)
    per_study_species_sites = defaultdict(set)
    per_study_species_populations = defaultdict(set)
    per_study_species_phylum = {}
    all_ecosystems = set()
    max_columns = len(cols)

    for values in reader:
        if not values or all(not x.strip() for x in values):
            continue
        all_records += 1
        if len(values) != max_columns:
            raise ValueError("Malformed row with unexpected number of CSV columns")
        # IMPORTANT: only allowlisted metadata fields are *used*.
        # This does not inspect or derive any genetic or plasticity outcome.
        a = {name: values[idx].strip() for name, idx in ix.items()}
        ecosystem = a["ecosystem"].lower()
        all_ecosystems.add(ecosystem)
        if ecosystem != "terrestrial":
            nonterrestrial[ecosystem or "missing"] += 1
            continue
        terrestrial_rows += 1
        phylum = a["phylum"].strip()
        if phylum not in animal_phyla:
            terrestrial_nonanimal_rows[phylum or "missing"] += 1
            terrestrial_nonanimal_taxa.add(a["taxon"])
            continue
        terrestrial_animal_phylum_rows[phylum] += 1
        species = " ".join(a["taxon"].replace("_", " ").split())
        study, pop = a["study"], a["source_population"]
        if not species or not study or not pop:
            missing_study_pop += 1
            continue
        terrestrial_species.add(species)
        try:
            lat = float(a["latitude"])
            lon = float(a["longitude"])
        except (ValueError, TypeError):
            missing_coords += 1
            continue
        if not math.isfinite(lat) or not math.isfinite(lon):
            missing_coords += 1
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            invalid_coords += 1
            continue
        key = (study, species)
        per_study_species_sites[key].add((lat, lon))
        per_study_species_populations[key].add(pop)
        per_study_species_phylum[key] = a["phylum"]
        per_population_coords[(study, species, pop)].add((lat, lon))

    geocoded_groups = {k for k, points in per_study_species_sites.items() if len(points) >= 2 and len(per_study_species_populations[k]) >= 2}
    geocoded_species = {k[1] for k in per_study_species_sites}
    qualified_species = {k[1] for k in geocoded_groups}
    qualified_studies = {k[0] for k in geocoded_groups}
    genera = {name.strip().split()[0] for name in qualified_species if len(name.strip().split()) >= 2}
    nonbinomial_taxa = sorted(x for x in terrestrial_species if len(x.strip().split()) < 2)
    terrestrial_label_pairs = sorted((str(study), str(taxon)) for study, taxon in per_study_species_sites)
    source_split = Counter()
    for name in sorted(qualified_species):
        z = hashlib.sha256(("ttf-q-evo-phenotype-v01|" + name.strip().lower()).encode()).digest()
        source_split["train" if z[0] & 1 else "heldout"] += 1
    split_counts = {"train": source_split["train"], "heldout": source_split["heldout"]}
    duplicate_conflicts = sum(len(coords) > 1 for coords in per_population_coords.values())
    conflicting_meta_keys = [list(key) for key, points in per_population_coords.items() if len(points) > 1]
    macro_ok = (len(qualified_species) >= 50 and len(qualified_studies) >= 15 and
                len(genera) >= 10 and min(split_counts.values()) >= 20 and
                duplicate_conflicts == 0)
    return {
        "schema": "ttf_barley_response_blind_terrestrial_origin_geometry_v01",
        "status": "MACRO_SOURCE_GEOMETRY_ELIGIBLE_RESPONSE_CLOSED" if macro_ok
            else "PILOT_ONLY_NO_MACRO_HISTORY_VS_PLASTICITY_CLAIM",
        "total_csv_data_rows": all_records,
        "source_header_columns": max_columns,
        "terrestrial_rows": terrestrial_rows,
        "terrestrial_nonanimal_rows_by_phylum": dict(sorted(terrestrial_nonanimal_rows.items())),
        "terrestrial_nonanimal_taxon_labels": sorted(terrestrial_nonanimal_taxa),
        "terrestrial_animal_rows_by_phylum": dict(sorted(terrestrial_animal_phylum_rows.items())),
        "other_ecosystem_rows": dict(sorted(nonterrestrial.items())),
        "ecosystem_labels": sorted(all_ecosystems),
        "terrestrial_species_any_coordinate_status": len(terrestrial_species),
        "terrestrial_species_with_at_least_one_georeferenced_record": len(geocoded_species),
        "study_species_groups_with_georeferenced_origin": len(per_study_species_sites),
        "original_study_species_groups_with_two_georeferenced_populations": len(geocoded_groups),
        "eligible_terrestrial_species_two_populations": len(qualified_species),
        "original_studies_with_two_georeferenced_populations": len(qualified_studies),
        "genera_among_eligible_species": len(genera),
        "terrestrial_taxon_metadata_labels": sorted(terrestrial_species),
        "nonbinomial_terrestrial_taxon_labels": nonbinomial_taxa,
        "metadata_only_study_taxon_groups": [list(x) for x in terrestrial_label_pairs],
        "eligible_species_hash_train_heldout": split_counts,
        "distinct_origin_population_keys": len(per_population_coords),
        "conflicting_origin_coordinates_for_same_pop_key": duplicate_conflicts,
        "metadata_only_conflicting_source_population_keys": sorted(conflicting_meta_keys),
        "metadata_only_amendment": amendment["schema"],
        "missing_study_species_population_metadata_rows": missing_study_pop,
        "terrestrial_missing_or_nonfinite_coord_rows": missing_coords,
        "terrestrial_out_of_bounds_coord_rows": invalid_coords,
        "author_source_meta_n_species": rules["source"]["publication_n_species"],
        "author_source_meta_n_studies": rules["source"]["publication_n_studies"],
        "author_source_meta_n_populations": rules["source"]["publication_n_populations"],
        "geometry_eligibility_only": True,
        "phenotype_ARR_values_evaluated": False,
        "thermal_limit_values_evaluated": False,
        "history_climate_extracted": False,
        "genetic_response_opened": False,
        "only_terrestrial_CHELSA_domain_allowed": True,
        "response_opening_authorized_by_this_census": False,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, help="Local source CSV, requires exact Git blob hash")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    c = json.loads(CONTRACT.read_text())
    if c["status"] != "PRE_ROW_LEVEL_ARR_VIEW_SOURCE_GEOGRAPHY_ONLY_PILOT":
        raise ValueError("Source-only screen protocol no longer frozen")
    src = c["source"]
    url = (f"https://raw.githubusercontent.com/{src['github_repo']}/"
           f"{src['github_commit']}/{src['source_file']}")
    if args.source:
        raw = args.source.read_bytes()
    else:
        request = urllib.request.Request(url, headers={"User-Agent": "TTF-external-pilot-source-only-census/0.1"})
        with urllib.request.urlopen(request, timeout=60) as response:
            raw = response.read(2_000_000)
    actual_sha1 = git_blob_sha1(raw)
    if actual_sha1 != src["source_git_blob_sha1"]:
        raise ValueError("Original Git blob identity mismatch BEFORE any response geometry audit")
    summary = source_census(raw, c)
    summary["source_git_blob_sha1_verified"] = actual_sha1
    summary["source_sha256"] = hashlib.sha256(raw).hexdigest()
    summary["source_commit"] = src["github_commit"]
    summary["source_url"] = url
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
