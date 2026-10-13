"""Data-free TDD fixtures for an exploratory CHELSA climate-only LCP pilot."""
import importlib.util
import heapq
from pathlib import Path

import numpy as np
import pytest
from affine import Affine

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'pilot_historical_host_climate_corridor.py'


def load():
    spec=importlib.util.spec_from_file_location('host_lcp_pilot', SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def brute_grid_cost(cost, length_east, length_north, a, b):
    """Independent explicit Dijkstra for small toy grids with fixed edge lengths."""
    rows,cols=cost.shape; q=[(0.0,*a)];best={a:0.0}
    while q:
        current,r,c=heapq.heappop(q)
        if (r,c)==b:return current
        if current>best[(r,c)]:continue
        for dr,dc,step in ((0,1,length_east),(0,-1,length_east),(1,0,length_north),(-1,0,length_north)):
            nr,nc=r+dr,c+dc
            if not (0<=nr<rows and 0<=nc<cols) or not np.isfinite(cost[nr,nc]):continue
            w=0.5*(cost[r,c]+cost[nr,nc])*step
            proposal=current+w
            if proposal<best.get((nr,nc),float('inf')):
                best[(nr,nc)]=proposal;heapq.heappush(q,(proposal,nr,nc))
    return float('inf')


def test_uniform_native_grid_lcp_is_symmetric_and_uses_real_kilometers():
    m=load()
    costs=np.ones((5,6),float)
    transform=Affine.translation(-85.0,11.0)*Affine.scale(.01,-.01)
    nodes=np.array([(0,0),(0,5),(4,5),(4,0)])
    result=m.cost_distance_between_pixels(costs,transform,nodes,np.array([(0,1),(1,2),(0,2),(0,0)]))
    assert result[0]>5 and result[0]<6
    assert result[1]>4 and result[1]<5
    assert result[2]==pytest.approx(result[0]+result[1],abs=.05)
    assert result[3]==0


def test_off_axis_corridor_changes_path_not_five_center_costs():
    m=load()
    bad=np.full((7,9),20.0)
    good=bad.copy(); good[2,:]=1.0
    for a in (bad,good):
        a[3,0]=1;a[3,8]=1
    transform=Affine.translation(0,0)*Affine.scale(.01,-.01)
    pts=np.array([(3,0),(3,8)])
    pairs=np.array([(0,1)])
    a=m.cost_distance_between_pixels(good,transform,pts,pairs)[0]
    b=m.cost_distance_between_pixels(bad,transform,pts,pairs)[0]
    assert 0<a<b
    assert np.allclose(good[3,1:8],bad[3,1:8])


def test_impossible_and_invalid_endpoints_fail_closed():
    m=load(); transform=Affine.translation(0,0)*Affine.scale(.01,-.01)
    cost=np.ones((5,5),float); cost[:,2]=np.inf
    with pytest.raises(m.ClimateCorridorNotEvaluable,match='unreachable'):
        m.cost_distance_between_pixels(cost,transform,np.array([(2,0),(2,4)]),np.array([[0,1]]))
    with pytest.raises(m.ClimateCorridorNotEvaluable,match='endpoint'):
        m.cost_distance_between_pixels(cost,transform,np.array([(2,0),(2,2)]),np.array([[0,1]]))
    with pytest.raises(ValueError,match='positive'):
        m.cost_distance_between_pixels(np.zeros((5,5)),transform,np.array([(1,0),(1,4)]),np.array([[0,1]]))


def test_bad_graph_geometry_and_row_length_cannot_be_repaired():
    m=load();c=np.ones((4,4));transform=Affine.translation(0,0)*Affine.scale(.01,-.01)
    with pytest.raises(ValueError,match='self|duplicate'):
        m.cost_distance_between_pixels(c,transform,np.array([(0,0),(1,1)]),np.array([[0,0]]),forbid_self_edges=True)
    with pytest.raises(ValueError,match='outside'):
        m.cost_distance_between_pixels(c,transform,np.array([(0,0),(5,1)]),np.array([[0,1]]))


def test_outward_native_window_and_bound_are_deterministic():
    m=load(); trans=Affine.translation(-90,20)*Affine.scale(.01,-.01)
    lon=np.array([-85.5100,-85.2804]);lat=np.array([10.455,11.05]);
    w=m.native_bounding_window(trans,lon,lat,0.12,shape=(1800,3600),maximum_cells=15000)
    assert int(w.width*w.height)<15000
    cropped=Affine.translation(0,0)*trans # identical full transform
    for x,y in zip(lon,lat):
        row,col=__import__('rasterio').transform.rowcol(trans,x,y)
        assert w.row_off<=row<w.row_off+w.height
        assert w.col_off<=col<w.col_off+w.width
    with pytest.raises(m.ClimateCorridorNotEvaluable,match='cell cap'):
        m.native_bounding_window(trans,lon,lat,0.12,shape=(1800,3600),maximum_cells=100)


def test_quantify_climate_cost_nearest_neighbor_finite_only():
    m=load()
    climate=np.array([[[0.,0.],[1.,1.]],[[2.,2.],[np.nan,0.]]])
    cloud=np.array([[0.,0.]])
    res=m.analogue_cell_cost(climate,cloud)
    assert res[0,0]==1
    assert res[0,1]==pytest.approx(1+2**.5)
    assert np.isinf(res[1,1]); assert np.all(np.isfinite(res[:1]))


def test_grid_cost_matches_independent_bruteforce_shortest_path():
    m=load()
    rng=np.random.default_rng(25)
    costs=1+rng.random((5,6))*3
    trans=Affine.translation(-85,11)*Affine.scale(.01,-.01)
    # Explicit local independent edge lengths at each latitude.
    lat0=11-.005
    east=np.array([m.haversine_km(-85+.005,lat0-r*.01,-85+.015,lat0-r*.01) for r in range(5)])
    north=m.haversine_km(-85+.005,lat0,-85+.005,lat0-.01)
    import heapq
    def explicit(a,b):
        q=[(0.,*a)];best={tuple(a):0.}
        while q:
            d,r,c=heapq.heappop(q)
            if (r,c)==tuple(b):return d
            if d>best[(r,c)]:continue
            for dr,dc in ((0,1),(0,-1),(1,0),(-1,0)):
                nr,nc=r+dr,c+dc
                if nr<0 or nr>=5 or nc<0 or nc>=6:continue
                step=east[r] if dr==0 else north
                nd=d+step*.5*(costs[r,c]+costs[nr,nc])
                if nd<best.get((nr,nc),float('inf')):
                    best[(nr,nc)]=nd;heapq.heappush(q,(nd,nr,nc))
        return float('inf')
    nodes=np.array([(0,0),(4,5),(2,1),(3,4)])
    edges=np.array([[0,1],[1,2],[3,0]])
    actual=m.cost_distance_between_pixels(costs,trans,nodes,edges)
    expected=[explicit(nodes[i],nodes[j]) for i,j in edges]
    np.testing.assert_allclose(actual,expected,atol=1e-10,rtol=1e-12)


def test_native_chelsa_patch_rejects_misaligned_variable_or_time(tmp_path):
    m=load()
    import rasterio
    from rasterio.transform import from_origin
    t=from_origin(-86.2,11.8,.01,.01)
    assets={}
    for ti in (-190,20):
        for v in m.VARS:
            p=tmp_path/f'CHELSA_TraCE21K_{v}_{ti}_V1.0.tif'
            with rasterio.open(p,'w',driver='GTiff',height=200,width=200,count=1,dtype='float32',crs='EPSG:4326',transform=t,nodata=-9999) as ds:
                ds.write(np.full((200,200),10,dtype='float32'),1)
            assets[(v,ti)]=p
    lo=np.array([-85.51,-85.28]);la=np.array([10.455,11.05])
    raw,transform,metadata=m.check_assets_and_read_native_window(assets,lo,la)
    assert set(raw)=={'0BP','LGM'}
    assert np.all(raw['LGM']==10)
    assert metadata['rows']*metadata['cols']<15000
    wrong=from_origin(-86.199,11.8,.01,.01)
    altered=assets[('bio15',-190)]
    with rasterio.open(altered,'w',driver='GTiff',height=200,width=200,count=1,dtype='float32',crs='EPSG:4326',transform=wrong,nodata=-9999) as ds:
        ds.write(np.full((200,200),10,dtype='float32'),1)
    with pytest.raises(m.ClimateCorridorNotEvaluable,match='different native grid'):
        m.check_assets_and_read_native_window(assets,lo,la)


def test_original_source_bindings_use_frozen_contract_field_names():
    m=load()
    import json
    rule=json.loads((SCRIPT.parents[1]/'docs/corridor_climate_pilot_contract_v0.1.json').read_text())
    assert m.verify_pilot_source_bindings(rule) is True
    damaged=json.loads(json.dumps(rule))
    damaged['source_sha256']['frozen_edge_predictors']='0'*64
    with pytest.raises(RuntimeError,match='source binding'):
        m.verify_pilot_source_bindings(damaged)
