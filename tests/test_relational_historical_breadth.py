import importlib.util
from pathlib import Path


RELATION = Path("scripts/build_relational_historical_relation.py")
OPPORTUNITY = Path("scripts/attach_relational_historical_opportunity.py")


def load(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gate():
    return {
        "largest_single_order_fraction_max":0.50,
        "order_fraction_threshold":0.05,
        "minimum_orders_at_or_above_fraction_threshold":4,
    }


def test_historical_relation_breadth_passes_balanced_four_orders():
    m=load(RELATION,"hist_relation")
    names=[f"sp{i}" for i in range(20)]
    orders=["A"]*8+["B"]*4+["C"]*4+["D"]*4
    metadata={name:{"order":order} for name,order in zip(names,orders)}
    out=m.taxonomic_breadth_summary(names,metadata,gate())
    assert out["pass"] is True
    assert out["largest_order_fraction"]==0.4


def test_historical_relation_breadth_fails_order_majority():
    m=load(RELATION,"hist_relation_fail")
    names=[f"sp{i}" for i in range(20)]
    orders=["A"]*11+["B"]*3+["C"]*3+["D"]*3
    metadata={name:{"order":order} for name,order in zip(names,orders)}
    out=m.taxonomic_breadth_summary(names,metadata,gate())
    assert out["pass"] is False


def test_historical_supported_breadth_fails_empty_group():
    m=load(OPPORTUNITY,"hist_opportunity")
    out=m.taxonomic_breadth_summary(set(),{},gate())
    assert out["pass"] is False
    assert out["species"]==0
