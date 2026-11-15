"""Synthetic / frozen-aggregate constraints; never opens the 2021 VCF."""
from pathlib import Path
import importlib.util
import json
import pytest

PATH=Path(__file__).resolve().parents[1]/"research"/"prospective"/"senecio_p1_2025_ri_phylogeny_postoutcome_crosswalk_20261115.py"
modspec=importlib.util.spec_from_file_location("senecio_2025_pair",PATH)
m=importlib.util.module_from_spec(modspec)
modspec.loader.exec_module(m)


def test_phylogeny_header_parses_only_declared_taxa():
    n="#NEXUS\nbegin taxa;\n dimensions ntax=3;\n taxlabels D01 H01 D03;\nend;\nbegin trees;\n tree t1=(D01,H01,D03);\nend;"
    assert m.read_tree_tips(n)=={"D01","H01","D03"}
    assert "D0" not in m.read_tree_tips(n)
    with pytest.raises(ValueError):
        m.read_tree_tips(n.replace("ntax=3","ntax=4"))


def source_ri(e="DH"):
    header=[
        "CrossType","RIseedSet-SobelChen","RIviability-SobelChen",
        "CrossType","Ecotype","Clade","Divergence time",
        "EnvDist-MDS-5PCs","PhenoDist-MDS-7PCs"
    ]
    rows=[["D01-H01","-0.134","0.018",e,"Different","Within","100","12","2"]]
    for j in range(17):
        rows.append([f"X{j}-Y{j}","0.1","-0.1","DD","Same","Between","100","2","1"])
    return ("\t".join(header)+"\n"+"\n".join("\t".join(r) for r in rows)+"\n").encode()


def test_duplicate_cross_type_header_is_positionally_handled():
    data=m.read_ri_pair_table(source_ri())
    assert len(data)==18
    assert data["D01-H01"]["ecotype_cross"]=="DH"
    assert data["D01-H01"]["ri_seed_set"] == -0.134


def test_duplicate_or_schema_drift_is_refused():
    source=source_ri()
    with pytest.raises(ValueError):
        m.read_ri_pair_table(source.replace(b"X1-Y1",b"X0-Y0"))
    with pytest.raises(ValueError):
        m.read_ri_pair_table(source.replace(b"RIseedSet-SobelChen",b"RIseedSet"))


def test_unexpected_genotype_or_ri_join_is_refused():
    p1=json.loads(m.P1_FILE.read_text())
    crosslayer=json.loads(m.P0_P1_FILE.read_text())
    ri=m.read_ri_pair_table(source_ri())
    tips={"D01","H01","D03","H02","D04","H05","D32","H12"}
    out=m.join_sources(p1,ri,tips,crosslayer)
    assert out["p1_pair_count"]==8
    assert out["pairs_with_2025_tree_both_tips"]==4
    assert out["pairs_with_direct_2025_ri_cross"]==1
    assert out["no_new_original_vcf_or_genotype_rows_read"] is True
    alt=dict(ri)
    alt["D04-H05"]={"ecotype_cross":"DD","ecotype_contrast":"Same","ri_seed_set":0.3,"ri_viability":0.5}
    with pytest.raises(ValueError):
        m.join_sources(p1,alt,tips,crosslayer)
