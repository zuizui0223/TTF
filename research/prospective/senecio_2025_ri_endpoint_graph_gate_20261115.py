"""Source-only TTF-Q graph-size bound for published Senecio 2025 crosses.

The source crossing outcomes are never decoded. A disjoint-endpoint matching
is an upper bound on pair comparisons without shared sampled populations,
NOT an estimate of independent evolutionary events or effective sample size.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import urllib.request
from functools import lru_cache

SOURCE="https://raw.githubusercontent.com/MaddieEJames/James_et_al._2025_MolecularEcology/fd8e8f075a72c1420267eecfd482dbc6a312caef/data/CrossesDH-ReciprocalsAveraged.txt"
SHA="ae6147df032baf492d692dea33aa03398b794d78"


def blob_sha(raw:bytes)->str:
    return hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()


def key_graph(raw:bytes)->dict:
    if blob_sha(raw)!=SHA: raise ValueError("original author Git blob mismatch")
    reader=csv.reader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t")
    header=next(reader)
    if len(header)!=9 or header[0]!="CrossType" or header[3]!="CrossType" or header[5]!="Clade":
        raise ValueError("source crossing data schema changed")
    edges=[]
    for cells in reader:
        if not cells:continue
        if len(cells)!=len(header):raise ValueError("unexpected crossing row")
        pair=cells[0]
        units=pair.split("-")
        if len(units)!=2 or units[0]==units[1]:
            raise ValueError("source crossing ID is not two distinct populations")
        edges.append({"pair":pair,"a":units[0],"b":units[1],"type":cells[3],"clade":cells[5]})
    if len(edges)!=18 or len({e["pair"] for e in edges})!=18:
        raise ValueError("source row count or duplicate pair key")
    nodes=sorted({v for e in edges for v in (e["a"],e["b"])})
    if len(nodes)!=8:raise ValueError("original source population count drift")
    # Dynamic programming over matched pairs, with no reuse of endpoint populations.
    nodeidx={v:i for i,v in enumerate(nodes)}
    @lru_cache(None)
    def solve(edge_start:int, mask:int)->tuple[str,...]:
        if edge_start==len(edges):return ()
        a=edges[edge_start]
        best=solve(edge_start+1,mask)
        bits=(1<<nodeidx[a["a"]])|(1<<nodeidx[a["b"]])
        if not mask&bits:
            yes=(a["pair"],)+solve(edge_start+1,mask|bits)
            if len(yes)>len(best):best=yes
        return best
    witness=solve(0,0)
    return {
       "status":"SOURCE_ONLY_TTF_RI_ENDPOINT_GEOMETRY",
       "source_git_blob":SHA,
       "crossing_comparison_rows":len(edges),
       "unique_population_nodes":len(nodes),
       "nonoverlapping_population_pair_maximum":len(witness),
       "example_maximum_disjoint_matching":list(witness),
       "species_evaluation_unit_count":1,
       "qualified_unseen_species_TTF":"NOT_EVALUABLE_SINGLE_SPECIES",
       "unseen_lineage_TTF":"NOT_EVALUABLE_ONLY_TWO_SOURCE_CLADES_AND_SHARED_ENDPOINTS",
       "graph_pair_independence_is_not_evolutionary_independence":True,
       "original_response_values_inspected":False,
    }


def main():
    data=urllib.request.urlopen(SOURCE,timeout=35).read(100_000)
    print(json.dumps(key_graph(data),sort_keys=True,indent=2))


if __name__=="__main__":main()
