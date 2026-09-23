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


def test_compact_geometry_loader_reconstructs_exact_species_slices(tmp_path):
    m=load_opportunity_script()
    names=np.asarray([f"sp{i:04d}" for i in range(1000)],dtype="U160")
    offsets=np.arange(1001,dtype=np.int64)
    mids=np.column_stack((
        np.arange(1000,dtype=float),
        np.zeros(1000,dtype=float),
        np.ones(1000,dtype=float),
    ))
    path=tmp_path/"compact.npz"
    np.savez_compressed(
        path,
        species_order=names,
        edge_offsets=offsets,
        edge_midpoints_ecef_km=mids,
        n_localities=np.full(1000,12,dtype=np.int64),
        graph_k=np.full(1000,2,dtype=np.int64),
        n_edges=np.ones(1000,dtype=np.int64),
        min_endpoint_disjoint_training_edges=np.full(1000,5,dtype=np.int64),
    )
    edge_map,locality_n=m.load_compact_edges(path)
    assert set(edge_map)==set(map(str,names))
    assert edge_map["sp0007"].shape==(1,3)
    assert np.array_equal(edge_map["sp0007"][0],np.asarray([7.0,0.0,1.0]))
    assert locality_n["sp0007"]==12


def test_compact_geometry_loader_rejects_offset_edge_count_drift(tmp_path):
    m=load_opportunity_script()
    names=np.asarray([f"sp{i:04d}" for i in range(1000)],dtype="U160")
    offsets=np.arange(1001,dtype=np.int64)
    path=tmp_path/"bad.npz"
    np.savez_compressed(
        path,
        species_order=names,
        edge_offsets=offsets,
        edge_midpoints_ecef_km=np.zeros((1000,3),dtype=float),
        n_localities=np.full(1000,12,dtype=np.int64),
        graph_k=np.full(1000,2,dtype=np.int64),
        n_edges=np.r_[np.asarray([2],dtype=np.int64),np.ones(999,dtype=np.int64)],
        min_endpoint_disjoint_training_edges=np.full(1000,5,dtype=np.int64),
    )
    import pytest
    with pytest.raises(RuntimeError,match="edge-count/offset"):
        m.load_compact_edges(path)
