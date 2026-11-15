"""Synthetic structural gate tests; no genetic or crossing outcome values."""
import importlib.util
from pathlib import Path
import pytest

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"senecio_2025_ri_endpoint_graph_gate_20261115.py"
spec=importlib.util.spec_from_file_location("ri_graph",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def generate():
    nodes=["D01","D04","D23","D32","H01","H05","H12","H21"]
    pairs=[(nodes[i],nodes[j]) for i in range(8) for j in range(i+1,8)][:18]
    header=["CrossType","RIseedSet-SobelChen","RIviability-SobelChen","CrossType",
            "Ecotype","Clade","Divergence time","EnvDist-MDS-5PCs","PhenoDist-MDS-7PCs"]
    rows=["\t".join(header)]
    for a,b in pairs: rows.append("\t".join([a+"-"+b,"0.5","0.3","DH","Different","Within","0","0","0"]))
    return ("\n".join(rows)+"\n").encode()

def audited(raw):
    old=m.SHA
    m.SHA=m.blob_sha(raw)
    try: return m.key_graph(raw)
    finally: m.SHA=old

def test_overlap_graph_is_never_18_independent_pairs():
    out=audited(generate())
    assert out["crossing_comparison_rows"]==18
    assert out["unique_population_nodes"]==8
    assert out["nonoverlapping_population_pair_maximum"]<=4
    assert out["qualified_unseen_species_TTF"]=="NOT_EVALUABLE_SINGLE_SPECIES"
    assert out["original_response_values_inspected"] is False

def test_modified_source_fails_closed():
    with pytest.raises(ValueError,match="Git blob"):
        m.key_graph(generate())

def test_duplicate_key_fails_closed():
    t=generate().decode("utf-8")
    lines=t.splitlines()
    lines[-1]=lines[1]
    with pytest.raises(ValueError,match="row count or duplicate"):
        audited(("\n".join(lines)+"\n").encode())
