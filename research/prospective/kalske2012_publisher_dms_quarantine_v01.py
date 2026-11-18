"""Fail-closed DMS validation of Kalske 2012 publisher-rendered site strings.

Original publisher HTML strings are preserved verbatim. Do not normalize
invalid minute/second fields, swap lat/lon, or borrow coordinates from
other studies of Vincetoxicum hirundinaria.
"""
import argparse,json,re
from pathlib import Path

PUBLISHED={
 "Ruissalo":("67° 02′ 57″ N","15° 63′ 73″ E"),
 "Lammasluoto":("66° 80′ 66″ N","15° 52′ 62″ E"),
 "Ånskär":("66° 76′ 06″ N","15° 38′ 90″ E"),
 "Utö":("66° 29′ 54″ N","15° 21′ 28″ E"),
}
PAT=re.compile(r"^(\d+)°\s*(\d+)′\s*(\d+)″\s*([NSEW])$")
def validate(s,axis):
    m=PAT.fullmatch(s)
    if not m:return ["FORMAT_INVALID"]
    deg,minute,second,hem=m.groups();deg=int(deg);minute=int(minute);second=int(second)
    errors=[]
    if minute>=60:errors.append("MINUTES_OUT_OF_RANGE")
    if second>=60:errors.append("SECONDS_OUT_OF_RANGE")
    if hem not in ("NS" if axis=="latitude" else "EW"):errors.append("HEMISPHERE_MISMATCH")
    if deg>(90 if axis=="latitude" else 180):errors.append("DEGREES_OUT_OF_RANGE")
    return errors

def evaluate():
    entries=[]
    for site,(lat,lon) in PUBLISHED.items():
        le=validate(lat,"latitude");oe=validate(lon,"longitude")
        entries.append({"site":site,"published_latitude":lat,"published_longitude":lon,
                        "latitude_issues":le,"longitude_issues":oe,
                        "valid_dms_pair":not(le or oe),"certified_source_GPS":False})
    if len(entries)!=4 or sum(x["valid_dms_pair"] for x in entries)!=0:
        raise AssertionError("publisher DMS qualification unexpectedly changed")
    return {"schema":"ttf_Kalske2012_publisher_DMS_quarantine_v01",
            "publisher_doi":"10.1111/j.1365-2745.2012.02008.x",
            "status":"ALL_FOUR_PUBLISHER_DMS_PAIRS_QUARANTINED",
            "published_site_count":4,
            "valid_publisher_DMS_pairs":0,
            "certified_original_source_GPS":0,
            "entries":entries,
            "no_coordinate_correction_or_cross_study_substitution":True,
            "no_fitness_genetics_LGM_or_adaptation_response_opened":True}
def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    a=p.parse_args();r=evaluate();a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":r["status"],"valid_publisher_DMS_pairs":0}))
if __name__=="__main__":main()
