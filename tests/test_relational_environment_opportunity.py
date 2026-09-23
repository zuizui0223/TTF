import importlib.util
from pathlib import Path

import numpy as np
from ttf.relational_environment import nearest_coverage


def test_directed_coverage_is_target_fraction():
    source=np.asarray([[0.,0.,0.],[100.,0.,0.]])
    target=np.asarray([[1.,0.,0.],[99.,0.,0.],[1000.,0.,0.],[1200.,0.,0.]])
    assert np.isclose(nearest_coverage(target,source,radius=10),0.5)
    # Reversing direction changes the denominator and therefore the estimand.
    assert np.isclose(nearest_coverage(source,target,radius=10),1.0)


SCRIPT = Path("scripts/attach_relational_environment_opportunity.py")


def load_opportunity_script():
    spec=importlib.util.spec_from_file_location("envopp",SCRIPT)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_supported_taxonomic_breadth_uses_inherited_thresholds():
    m=load_opportunity_script()
    names={f"sp{i}" for i in range(20)}
    ordered=sorted(names)
    orders=["A"]*8+["B"]*4+["C"]*4+["D"]*4
    metadata={name:{"order":order} for name,order in zip(ordered,orders)}
    gate={
        "largest_single_order_fraction_max":0.50,
        "order_fraction_threshold":0.05,
        "minimum_orders_at_or_above_fraction_threshold":4,
    }
    out=m.taxonomic_breadth_summary(names,metadata,gate)
    assert out["pass"] is True
    assert out["largest_order_fraction"]==0.4


def test_supported_taxonomic_breadth_fails_empty_group():
    m=load_opportunity_script()
    gate={
        "largest_single_order_fraction_max":0.50,
        "order_fraction_threshold":0.05,
        "minimum_orders_at_or_above_fraction_threshold":4,
    }
    out=m.taxonomic_breadth_summary(set(),{},gate)
    assert out["pass"] is False
    assert out["species"]==0
