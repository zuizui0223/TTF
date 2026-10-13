"""External predictor/geometry only: no nucleotide sequences or genetic responses."""
import importlib.util
import numpy as np
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("host_sampling", ROOT / "scripts" / "audit_host_memory_sampling_geometry.py")


def lib():
    module = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(module)
    return module


def geom(seed=114):
    rng=np.random.default_rng(seed)
    xyz=rng.normal(size=(18,3))*130
    nodes=np.array([(a,b) for a in range(18) for b in range(a+1,18)
                    if (a+b) % 3 != 0],dtype=int)
    return xyz,nodes


def test_sampling_isolation_is_geographic_and_independent_of_edge_order():
    m=lib();xyz,nodes=geom()
    a=m.sampling_edge_features(xyz,nodes)
    b=m.sampling_edge_features(xyz,nodes[::-1])
    for k in ("chord_km","endpoint_nearest_neighbor_km","endpoint_graph_degree"):
        np.testing.assert_allclose(a[k],b[k][::-1],atol=1e-12)
    assert np.all(a["chord_km"]>0)
    assert np.all(a["endpoint_nearest_neighbor_km"]>0)


def test_sharp_conditional_overlap_agrees_with_explicit_OLS():
    m=lib();rng=np.random.default_rng(80);xyz,nodes=geom()
    g=m.sampling_edge_features(xyz,nodes);n=len(nodes)
    own=rng.normal(size=n);present=rng.normal(size=n)
    hist=0.40*own-0.20*present + 0.12*np.log1p(g["chord_km"]) + 2.4*np.log1p(g["endpoint_nearest_neighbor_km"]) +rng.normal(size=n)*0.17
    obs=m.partial_sampling_overlap(hist,own,present,g)
    standard=lambda a:(a-a.mean())/a.std(ddof=0)
    y=standard(hist);A=np.column_stack((np.ones(n),standard(own),standard(present),standard(np.log1p(g["chord_km"]))))
    iso=standard(np.log1p(g["endpoint_nearest_neighbor_km"]))
    ry=y-A@np.linalg.lstsq(A,y,rcond=None)[0]
    ri=iso-A@np.linalg.lstsq(A,iso,rcond=None)[0]
    r2=float((ry@ri)**2/((ry@ry)*(ri@ri)))
    assert obs["nearest_neighbor_partial_r2"]==pytest.approx(r2,abs=1e-11)
    assert obs["sampling_geometry_joint_partial_r2"]+1e-12>=r2
    assert 0<=obs["nearest_neighbor_partial_r2"]<=1


def test_zero_independent_geometry_signal_is_not_fake_zero():
    m=lib();rng=np.random.default_rng(9);xyz,nodes=geom();g=m.sampling_edge_features(xyz,nodes)
    n=len(nodes);own=rng.normal(size=n);present=rng.normal(size=n)
    iso=np.log1p(g["endpoint_nearest_neighbor_km"])
    # Make M_host equal a linear combination of controls: no uniquely identifiable predictor.
    hist=0.22*own+0.4*present+0.2*np.log1p(g["chord_km"])
    with pytest.raises(ValueError,match="history information"):
        m.partial_sampling_overlap(hist,own,present,g)


def test_duplicate_graph_edges_are_not_silently_repaired():
    m=lib();xyz,nodes=geom()
    bad=nodes.copy();bad[-1]=nodes[0][::-1]
    with pytest.raises(ValueError,match="duplicate"):
        m.sampling_edge_features(xyz,bad)


def test_zero_sampling_variation_is_explicitly_not_estimable():
    m=lib();rng=np.random.default_rng(10)
    n=60
    own=rng.normal(size=n);pres=rng.normal(size=n);host=rng.normal(size=n)
    g={"chord_km":np.linspace(10,100,n),"endpoint_nearest_neighbor_km":np.full(n,5.0),"endpoint_graph_degree":rng.integers(2,5,size=n)}
    with pytest.raises(ValueError,match="sampling isolation"):
        m.partial_sampling_overlap(host,own,pres,g)


def test_proxy_is_not_ancient_corridor_connectivity():
    m=lib()
    assert m.PROXY_SCOPE=="historical_host_climate_analogue_distance_not_ancient_host_occupancy_or_corridor"
    assert m.ANALYSIS_STATUS=="POST_OUTCOME_RESPONSE_BLIND_PREDICTOR_GEOMETRY_AUDIT"


def test_sampling_r2_degrees_of_freedom_guard_small_graph_bias():
    m=lib();rng=np.random.default_rng(234);xyz,nodes=geom()
    g=m.sampling_edge_features(xyz,nodes)
    n=len(nodes)
    y=rng.normal(size=n); own=rng.normal(size=n);current=rng.normal(size=n)
    output=m.partial_sampling_overlap(y,own,current,g)
    standard=lambda a:(a-a.mean())/a.std(ddof=0)
    C=np.column_stack([np.ones(n),standard(own),standard(current),standard(np.log1p(g['chord_km']))])
    iso=standard(np.log1p(g['endpoint_nearest_neighbor_km']))
    degree=standard(np.log1p(g['endpoint_graph_degree']))
    B=np.column_stack((C,iso,degree))
    yy=standard(y)
    rss_base=np.linalg.norm(yy-C@np.linalg.lstsq(C,yy,rcond=None)[0])**2
    rss_joint=np.linalg.norm(yy-B@np.linalg.lstsq(B,yy,rcond=None)[0])**2
    expected=1 - (rss_joint/(n-np.linalg.matrix_rank(B)))/(rss_base/(n-4))
    assert output['sampling_geometry_df_adjusted_partial_r2']==pytest.approx(expected,abs=1e-10)
    assert output['nearest_neighbor_df_adjusted_partial_r2']<=output['nearest_neighbor_partial_r2'] + 1e-10


def test_source_audit_reports_species_equal_weight_and_graph_size_strata():
    m=lib();rng=np.random.default_rng(43);xyz,all_nodes=geom()
    preds=[];edges=[];locations=[];roles=[]
    for s,panel,num in (("S1","development",25),("S2","confirmatory",75)):
        roles.append({"species":s,"panel":panel})
        for j,pt in enumerate(xyz):
            locations.append({"species":s,"locality_index":str(j),
                              "x_km":str(pt[0]),"y_km":str(pt[1]),"z_km":str(pt[2])})
        for j,(left,right) in enumerate(all_nodes[:num]):
            mh=float(rng.normal()); self_c=float(rng.normal());present=float(rng.normal())
            preds.append({"species":s,"edge_index":str(j),"M_host":str(mh),
                          "M_self":str(self_c),"D_host_0BP":str(present),"D_host_LGM":str(present+mh)})
            edges.append({"species":s,"edge_index":str(j),"node_left":str(left),"node_right":str(right)})
    report=m.audit(preds,edges,locations,roles)
    assert report['source_bound_species']==2
    assert report['response_blind_external_edge_rows']==100
    assert report['panels']['all']['nearest_neighbor_partial_r2']['species']==2
    assert len(report['original_graph_edge_count_strata'])==3
    assert sum(g['species'] for g in report['original_graph_edge_count_strata'])==2
    assert report['interpretation']['coefficient_or_p_value_tested'] is False
    assert 'S1' not in str(report) and 'S2' not in str(report)
