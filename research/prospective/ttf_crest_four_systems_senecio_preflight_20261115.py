"""Source-only Senecio preflight. Never decode phenotypes or genomic responses."""
from __future__ import annotations
import collections
import csv
import hashlib
import io
import json
import urllib.request

SOURCE_REPO = "OB-lab/James_et_al._2021-Evolution"
SOURCE_COMMIT = "33222710aca75f79153e4118b409e7d77787d045"
SOURCE_BLOBS = {
    "phenotype/data_files/allTraits.txt": "30f3edbdb486733306529e58135e6f1c205c0184",
    "phenotype/data_files/allTraitsPairs.txt": "ac8183378ea0778fa1a2fd42fd6a6d94d7008a85",
    "phenotype/data_files/traitMeansPairs.txt": "400f73bf21fefb6fcb775396f256c3447c755595",
}
GENOTYPE_OUTLIER_PAIRS = {
    "D00-H00", "D01-H01", "D02-H04", "D03-H02", "D04-H05", "D05-H06",
    "D12-H14", "D14-H15", "D32-H12",
}

def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()

def parse_source_keys(data: bytes, path: str) -> dict:
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")), delimiter="\t")
    fields = set(reader.fieldnames or ())
    individual = path.endswith(("allTraits.txt", "allTraitsPairs.txt"))
    required = {"Pair", "Ecotype", "Population.ID"} if individual else {"Pair", "Ecotype", "Group.1"}
    if not required <= fields:
        raise ValueError(f"missing metadata columns: {required - fields}")
    n = 0
    pops, pairs = collections.Counter(), collections.Counter()
    pair_ecotypes = collections.defaultdict(set)
    for row in reader:
        n += 1
        pair = (row.get("Pair") or "").strip()
        if pair:
            pairs[pair] += 1
            pair_ecotypes[pair].add((row.get("Ecotype") or "").strip())
        if individual:
            pop = (row.get("Population.ID") or "").strip()
            if not pop:
                raise ValueError("missing population ID")
            pops[pop] += 1
        # Deliberately do NOT inspect numeric outcome columns.
    return {
        "rows": n, "population_ids": sorted(pops),
        "pair_counts": dict(sorted(pairs.items())),
        "pair_ecotypes": {k: sorted(v) for k, v in sorted(pair_ecotypes.items())},
    }

def audit(fetch_fn=None) -> dict:
    if fetch_fn is None:
        def fetch_fn(path: str) -> bytes:
            u = f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{path}"
            with urllib.request.urlopen(u, timeout=35) as r:
                return r.read()
    parsed = {}
    for path, expected_sha in SOURCE_BLOBS.items():
        data = fetch_fn(path)
        if git_blob_sha(data) != expected_sha:
            raise ValueError(f"Git blob SHA mismatch: {path}")
        parsed[path] = parse_source_keys(data, path)
    full, paired, means = (parsed[k] for k in SOURCE_BLOBS)
    assert (full["rows"], len(full["population_ids"])) == (605, 20)
    assert (paired["rows"], len(paired["population_ids"])) == (481, 16)
    assert means["rows"] == 16
    if not set(paired["population_ids"]) <= set(full["population_ids"]):
        raise ValueError("paired population not present in source")
    if set(paired["pair_counts"]) != set(means["pair_counts"]):
        raise ValueError("pair means do not align with individual pair IDs")
    if len(paired["pair_counts"]) != 8:
        raise ValueError("expected 8 morphological habitat pairs")
    for pair, ecotypes in paired["pair_ecotypes"].items():
        if ecotypes != ["Dune", "Headland"]:
            raise ValueError(f"non-contrasting pair: {pair}")
    return {
        "status": "SOURCE_ONLY_NO_BIOLOGICAL_RESPONSE_OPENED",
        "source_commit": SOURCE_COMMIT,
        "source_blobs": SOURCE_BLOBS,
        "full_phenotype_individual_count": full["rows"],
        "full_phenotype_population_count": len(full["population_ids"]),
        "paired_phenotype_individual_count": paired["rows"],
        "paired_phenotype_population_count": len(paired["population_ids"]),
        "paired_locality_ids": sorted(paired["pair_counts"]),
        "paired_locality_counts": paired["pair_counts"],
        "excluded_allopatric_population_ids": sorted(set(full["population_ids"]) - set(paired["population_ids"])),
        "genotype_outlier_pair_file_ids": sorted(GENOTYPE_OUTLIER_PAIRS),
        "genotype_only_pair_ids": sorted(GENOTYPE_OUTLIER_PAIRS - set(paired["pair_counts"])),
        "no_individual_cross_modal_join_verified": True,
        "no_response_scores_or_genotypes_read": True,
    }

if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
