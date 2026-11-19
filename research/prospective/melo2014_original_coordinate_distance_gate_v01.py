"""Melo et al. 2014 original source-coordinate versus narrative distance audit.

Original article DOI 10.1111/nph.12779, Methods 'The system'.
Published coordinates are transcribed verbatim, without repair or proxy.
Source-only: no fitness, genetic, LGM, phenotypic or RI response is read.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHORTLIST = HERE / "hargreaves35_primary_article_recovery_shortlist_v01.json"
LABEL = "Meloetal_NewPhytol_2014"
DOI = "10.1111/nph.12779"
ORIGINAL_METHOD_COORDINATES = {
    "Dune": {"latitude": '28°19′54.66″S', "longitude": '153°34′17.04″E'},
    "Headland": {"latitude": '28°21′45.07″S', "longitude": '153°34′46.82″E'},
}
# Bernal-Franco 2014, University of Queensland PhD thesis, Chapter III Table 1
# (printed p. 122); rounded site locations, a different research sample.
# https://files01.core.ac.uk/download/pdf/43363386.pdf
THESIS_SITE_REFERENCE = {
    "Dune": {"population_code": "SD03", "place_name": "Cabarita Beach",
             "latitude": '28°19′54″S', "longitude": '153°34′17″E'},
    "Headland": {"population_code": "RH02", "place_name": "Hastings Point",
                 "latitude": '28°21′45″S', "longitude": '153°34′46″E'},
}
THESIS_SOURCE = "Bernal-Franco 2014 University of Queensland PhD, Chapter III Table 1"
PUBLISHED_APPROX_SEPARATION_KM = 0.8
RADIUS_KM = 6371.0088
DMS_RE = re.compile(r"^(\d+)°(\d+)′(\d+(?:\.\d+)?)″([NSEW])$")


def parse_dms(text, axis):
    m = DMS_RE.fullmatch(text)
    if m is None:
        raise ValueError("Published DMS syntax not parseable")
    deg, mins, secs, hem = m.groups()
    deg, mins, secs = int(deg), int(mins), float(secs)
    if not 0 <= mins < 60 or not 0 <= secs < 60:
        raise ValueError("Published DMS components out of range")
    allowed = "NS" if axis == "latitude" else "EW"
    limit = 90 if axis == "latitude" else 180
    if hem not in allowed or deg > limit or (deg == limit and (mins or secs)):
        raise ValueError("Published DMS hemisphere/degree invalid")
    return (deg + mins / 60 + secs / 3600) * (-1 if hem in ("S", "W") else 1)


def haversine(a, b):
    (lat1, lon1), (lat2, lon2) = a, b
    dlat, dlon = math.radians(lat2-lat1), math.radians(lon2-lon1)
    k = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return 2 * RADIUS_KM * math.asin(math.sqrt(min(1.0, max(0.0, k))))


def evaluate(shortlist):
    if shortlist.get("schema") != "ttf_Hargreaves35_original_primary_literature_recovery_shortlist_v01":
        raise ValueError("Frozen primary bibliography schema changed")
    matches = [r for r in shortlist["record_list"]
               if r["original_Hargreaves_study_label"] == LABEL]
    if len(matches) != 1:
        raise ValueError("Melo original label not unique")
    item = matches[0]
    if (item["primary_article_doi"] != DOI
        or item["distinct_original_source_ID_labels"] != 2
        or item["distinct_original_garden_ID_labels"] != 2
        or item["source_only_candidate_same_garden_cells"] != 10):
        raise ValueError("Melo original-vs-Hargreaves metadata drift")
    if any(item[k] for k in (
        "primary_experiment_identity_certified",
        "historical_climate_contrast_qualified",
        "genetic_adaptation_plasticity_speciation_effect_opened"
    )):
        raise ValueError("Original source-only admission gate changed")
    if set(ORIGINAL_METHOD_COORDINATES) != {"Dune", "Headland"}:
        raise ValueError("Original pair changed")
    pts = {label: [parse_dms(d["latitude"], "latitude"),
                   parse_dms(d["longitude"], "longitude")]
           for label, d in ORIGINAL_METHOD_COORDINATES.items()}
    thesis_site_concordance = {}
    if set(THESIS_SITE_REFERENCE) != {"Dune", "Headland"}:
        raise ValueError("Companion thesis pair changed")
    for label, point in pts.items():
        reference = THESIS_SITE_REFERENCE[label]
        secondary = [parse_dms(reference["latitude"], "latitude"),
                     parse_dms(reference["longitude"], "longitude")]
        delta = [round(abs(p-q)*3600.0, 6) for p,q in zip(point,secondary)]
        if any(x > 1.000001 for x in delta):
            raise ValueError("Cross-paper published coordinates no longer concordant")
        thesis_site_concordance[label] = {
            "thesis_population_code": reference["population_code"],
            "thesis_place_name": reference["place_name"],
            "thesis_latitude": reference["latitude"],
            "thesis_longitude": reference["longitude"],
            "latitude_difference_arcsec": delta[0],
            "longitude_difference_arcsec": delta[1],
            "matches_original_published_DMS_within_one_arcsecond": True,
            "same_2009_seed_harvest_or_2010_transplant_cohort_certified": False,
        }
    d_km = haversine(pts["Dune"], pts["Headland"])
    ratio = d_km / PUBLISHED_APPROX_SEPARATION_KM
    if not (3.4 < d_km < 3.6 and 4 < ratio < 4.6):
        raise AssertionError("Known published coordinate-vs-prose discrepancy drift")
    return {
        "schema": "ttf_Melo2014_original_Cabarita_coordinate_distance_conflict_v01",
        "article_doi": DOI,
        "source": "original Methods The system (published coordinates and c. 800 m separation)",
        "status": "PRIMARY_PUBLISHED_COORDINATES_CONFLICT_WITH_DESCRIBED_DISTANCE_HOLD",
        "original_hargreaves_source_ID_labels": 2,
        "original_hargreaves_garden_ID_labels": 2,
        "hargreaves_candidate_comparison_cells_not_independent_experiments": 10,
        "source_populations_for_focal_field_transplant": 2,
        "original_DMS_verbatim": ORIGINAL_METHOD_COORDINATES,
        "publisher_coordinates_DMS_syntax_valid": True,
        "companion_thesis_site_reference": THESIS_SOURCE,
        "companion_thesis_source_url": "https://files01.core.ac.uk/download/pdf/43363386.pdf",
        "published_site_name_and_seconds_concordance": thesis_site_concordance,
        "companion_thesis_labels_headland_as_Hastings_Point": True,
        "same_original_transplant_individuals_crosspaper_verified": False,
        "crosspaper_site_name_match_is_not_independent_original_GPS": True,
        "calculated_coordinate_points_wgs84_decimal_degrees": pts,
        "great_circle_distance_km_spherical": round(d_km, 6),
        "text_approx_separation_km": PUBLISHED_APPROX_SEPARATION_KM,
        "distance_difference_km": round(d_km-PUBLISHED_APPROX_SEPARATION_KM, 6),
        "coordinate_to_text_distance_ratio": round(ratio, 6),
        "published_points_might_be_representative_not_collection_GPS": True,
        "independently_verified_original_source_GPS": 0,
        "historical_climate_raster_extraction_qualified": False,
        "potential_coordinate_error_proven": False,
        "single_paper_standalone_modern_plus_LGM_effect_estimable": False,
        "fitness_genetics_adaptation_plasticity_speciation_LGM_response_opened": False,
        "old_Bontrager_44_of_50_stop_unchanged": True,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    r = evaluate(json.loads(SHORTLIST.read_text(encoding="utf-8")))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(r, sort_keys=True, indent=2, ensure_ascii=False)+"\n",
                        encoding="utf-8")
    print(json.dumps({k: r[k] for k in (
        "status", "great_circle_distance_km_spherical",
        "text_approx_separation_km", "distance_difference_km",
        "coordinate_to_text_distance_ratio", "independently_verified_original_source_GPS"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
