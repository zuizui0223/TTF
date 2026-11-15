"""P1 offline known-truth guards before real genetic outcome opening."""
import importlib.util
from pathlib import Path

FILE = Path(__file__).resolve().parents[1] / "research" / "prospective" / "senecio_p1_genomic_holdout_20261115.py"
spec = importlib.util.spec_from_file_location("senecio_p1_genomic", FILE)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)


def example_loci(n=125):
    pop_names={p for pair in g.PAIRS for p in pair}
    loci=[]
    for j in range(n):
        freq={}
        calls={}
        for d,h in g.PAIRS:
            freq[d]=0.25+(j%4)*0.01
            freq[h]=0.75-(j%3)*0.01
            calls[d]=30
            calls[h]=30
        loci.append({"chrom":f"ctg_{j:04d}","pos":100+j,"freq":freq,"n":calls,"calls":bytearray()})
    return loci


def test_heldout_population_genotypes_never_select_markers():
    records=example_loci()
    original=g.train_markers(records,0)
    initial=[(r["chrom"],r["pos"]) for _,r,_,_ in original]
    assert len(initial)==100
    for r in records:
        r["freq"]["D00"]=1
        r["freq"]["H00"]=0
        r["n"]["D00"]=0
        r["n"]["H00"]=0
    revised=g.train_markers(records,0)
    assert [(r["chrom"],r["pos"]) for _,r,_,_ in revised] == initial
    assert all(direction==1 for _,_,direction,_ in revised)


def test_one_per_contig_and_fixed_minimum():
    records=example_loci(125)
    records.append(dict(records[0],pos=500000))
    marker=g.train_markers(records,2)
    assert len(marker)==100
    assert len({r["chrom"] for _,r,_,_ in marker})==100
    assert len(g.train_markers(example_loci(99),1))==99


def test_gt_parser_diploid_and_unknowns():
    assert g.parse_gt("1|0:35:60") == 1
    assert g.parse_gt("9:0/1:10",1) == 1
    assert g.parse_gt("0/0") == 0
    assert g.parse_gt("1/1") == 2
    assert g.parse_gt("./.")==255
    assert g.parse_gt("0/.")==255
    try:
        g.parse_gt("0/2")
    except ValueError:
        pass
    else:
        raise AssertionError("nonbinary allele should not silently contribute")


def test_auc_ties_and_pair_bootstrap_determinism():
    assert abs(g.auc_of([1,1],[0,1])-0.75)<1e-12
    assert abs(g.auc_of([1,2],[0,0])-1.0)<1e-12
    a=g.bootstrap_pair_auc([.1,.3,.5,.7],n=200,seed=33)
    assert a==g.bootstrap_pair_auc([.1,.3,.5,.7],n=200,seed=33)
    assert len(a)==2 and a[0]<=a[1]


def test_wrong_original_archive_is_not_genomic_outcome():
    try:
        g.read_original_vcf(b"bogus")
    except ValueError as e:
        assert "Git blob hash" in str(e)
    else:
        raise AssertionError("non-source archive accepted")
