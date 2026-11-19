"""Kalske 2012 source-only DMS syntax AND geographic plausibility gate (v02).

Does not repair or substitute original coordinates. The conservative
screening box is *not* a Finland boundary polygon and is never proof of
a source GPS. No environmental rasters or biological outcomes are read.
"""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path
from kalske2012_publisher_dms_quarantine_v01 import PUBLISHED

# A deliberately generous bounding screen, NOT a polygon or verified site GPS.
FINLAND_SCREEN = {"lat_min": 59.0, "lat_max": 71.0,
                  "lon_min": 19.0, "lon_max": 32.0}
DMS = re.compile(r"^(\d+)°\s*(\d+)′\s*(\d+)″\s*([NSEW])$")


def parse_strict_dms(text: str, axis: str):
    match = DMS.fullmatch(text)
    if not match:
        return None, ["INVALID_FORMAT"]
    degrees, minutes, seconds, hemisphere = match.groups()
    degrees, minutes, seconds = int(degrees), int(minutes), int(seconds)
    problems = []
    if minutes >= 60:
        problems.append("MINUTES_OUT_OF_RANGE")
    if seconds >= 60:
        problems.append("SECONDS_OUT_OF_RANGE")
    if degrees > (90 if axis == "latitude" else 180):
        problems.append("DEGREES_OUT_OF_RANGE")
    allowed = "NS" if axis == "latitude" else "EW"
    if hemisphere not in allowed:
        problems.append("HEMISPHERE_MISMATCH")
    if problems:
        return None, problems
    value = degrees + minutes / 60.0 + seconds / 3600.0
    return (-value if hemisphere in "SW" else value), []


def assess_site(name: str, lat_text: str, lon_text: str):
    lat, lat_errors = parse_strict_dms(lat_text, "latitude")
    lon, lon_errors = parse_strict_dms(lon_text, "longitude")
    syntax_valid = not (lat_errors or lon_errors)
    inside_screen = (
        FINLAND_SCREEN["lat_min"] <= lat <= FINLAND_SCREEN["lat_max"]
        and FINLAND_SCREEN["lon_min"] <= lon <= FINLAND_SCREEN["lon_max"]
    ) if syntax_valid else False
    return {
        "site": name,
        "publisher_latitude_verbatim": lat_text,
        "publisher_longitude_verbatim": lon_text,
        "syntax_valid": syntax_valid,
        "syntax_issues": {"latitude": lat_errors, "longitude": lon_errors},
        "converted_latitude": round(lat, 6) if syntax_valid else None,
        "converted_longitude": round(lon, 6) if syntax_valid else None,
        "inside_conservative_finland_envelope": inside_screen,
        "original_source_GPS_certified": False,
        "admission": "QUARANTINE",
    }


def evaluate():
    if set(PUBLISHED) != {"Ruissalo", "Lammasluoto", "Ånskär", "Utö"}:
        raise ValueError("source study population set changed")
    entries = [assess_site(name, *coords) for name, coords in PUBLISHED.items()]
    if sum(x["syntax_valid"] for x in entries) != 1:
        raise AssertionError("publisher coordinate syntax status changed")
    if any(x["inside_conservative_finland_envelope"] for x in entries):
        raise AssertionError("geographic screen unexpected: re-audit source")
    if any(x["original_source_GPS_certified"] for x in entries):
        raise AssertionError("no source GPS qualification allowed")
    return {
        "schema": "ttf_Kalske2012_publisher_DMS_geographic_quarantine_v02",
        "source": "https://doi.org/10.1111/j.1365-2745.2012.02008.x",
        "source_scope": "publisher Study Populations coordinate strings, no outcome data",
        "screen": FINLAND_SCREEN,
        "screen_is_conservative_not_an_authoritative_boundary": True,
        "reported_study_region": "southwestern Finnish archipelago",
        "sites": entries,
        "site_count": 4,
        "syntax_invalid_count": 3,
        "syntax_valid_but_geographically_inconsistent_count": 1,
        "certified_original_source_GPS_count": 0,
        "climate_grid_extraction_allowed": False,
        "geographic_repair_or_other_study_substitution": False,
        "fitness_genetics_adaptation_plasticity_speciation_LGM_response_opened": False,
        "status": "ALL_FOUR_PUBLISHED_COORDINATES_QUARANTINED",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = evaluate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, sort_keys=True,
                                      indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in (
        "site_count", "syntax_invalid_count",
        "syntax_valid_but_geographically_inconsistent_count",
        "certified_original_source_GPS_count", "climate_grid_extraction_allowed",
        "status")}, sort_keys=True))


if __name__ == "__main__":
    main()
