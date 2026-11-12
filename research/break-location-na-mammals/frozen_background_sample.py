"""Stable capped within-survey background selection; NO response or climate used.

Input is the author's 20km grid after ≥50km uncertain break exclusion.
Raw sampled coordinates must remain ephemeral due IUCN-derived range licensing.
"""
from __future__ import annotations
import hashlib,math
from collections.abc import Iterable

PREFIX="jensen20km-background-v01|"
MAX_PER_SPECIES=1000
def hash_key(scientific_name, longitude,latitude):
    species=str(scientific_name).strip().lower()
    lon=float(longitude);lat=float(latitude)
    if len(species.split())<2 or not (-180<=lon<=180 and -90<=lat<=90):
        raise ValueError("invalid prequalified species or site")
    token=(PREFIX+species+"|"+format(lon,".6f")+"|"+format(lat,".6f"))
    return hashlib.sha256(token.encode("utf8")).digest()

def select_background(scientific_name,coordinates,max_sites=MAX_PER_SPECIES):
    """Treat exact 6-decimal coordinates as identity; freeze before 21ka climate."""
    if max_sites!=MAX_PER_SPECIES:raise ValueError("frozen 1000-site cap cannot be changed")
    candidates={}
    for longitude,latitude in coordinates:
        lon=float(longitude);lat=float(latitude)
        if not(math.isfinite(lon) and math.isfinite(lat)):
            raise ValueError("invalid coordinates")
        row=(round(lon,6),round(lat,6))
        if row in candidates:raise ValueError("duplicate frozen 20km grid site")
        candidates[row]=hash_key(scientific_name,lon,lat)
    return sorted(candidates, key=lambda p:(candidates[p],p))[:max_sites]
