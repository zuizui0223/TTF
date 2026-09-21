from ttf.relational_external import accepted_gbif_species_match, deterministic_greedy_thin, shard_items


def test_strict_gbif_species_match():
    ok, rule = accepted_gbif_species_match(
        "Bombus terrestris",
        {"usageKey": 1, "matchType": "EXACT", "rank": "SPECIES", "canonicalName": "Bombus terrestris"},
    )
    assert ok and rule == "exact_canonical_species"
    bad, _ = accepted_gbif_species_match(
        "Bombus terrestris",
        {"usageKey": 1, "matchType": "FUZZY", "rank": "SPECIES", "canonicalName": "Bombus terrestris"},
    )
    assert not bad


def test_deterministic_thinning_is_key_ordered_and_coordinate_deduplicated():
    rows = [
        {"source_key": 9, "latitude": 0.0, "longitude": 0.0},
        {"source_key": 3, "latitude": 0.0, "longitude": 0.0},
        {"source_key": 4, "latitude": 0.0, "longitude": 0.2},
    ]
    retained = deterministic_greedy_thin(rows, minimum_distance_km=10.0)
    assert [x["source_key"] for x in retained] == [3, 4]


def test_sharding_is_disjoint_and_complete():
    items = list(range(19))
    parts = [shard_items(items, shard_index=i, shards=4) for i in range(4)]
    assert sorted(x for part in parts for x in part) == items
    assert sum(len(p) for p in parts) == len(items)
