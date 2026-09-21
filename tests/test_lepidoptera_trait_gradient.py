import numpy as np
from ttf.lepidoptera_trait_gradient import (
    composite_trait_similarity,
    equal_target_gradient,
    residualize_within_target,
)

def test_composite_is_symmetric_and_self_maximal():
    wing={"a":1.0,"b":3.0,"c":9.0}; vol={"a":"U","b":"B","c":"M"}; host={"a":1,"b":2,"c":5}
    hab={"a":("open","edge","xeric","dist"),"b":("open","edge","mesic","dist"),"c":("closed","avoid","mesic","avoid")}
    ab=composite_trait_similarity("a","b",wing_size=wing,voltinism=vol,host_breadth=host,habitat=hab)
    ba=composite_trait_similarity("b","a",wing_size=wing,voltinism=vol,host_breadth=host,habitat=hab)
    aa=composite_trait_similarity("a","a",wing_size=wing,voltinism=vol,host_breadth=host,habitat=hab)
    assert ab==ba and aa==1.0 and 0<=ab<=1

def test_geometry_residualization_removes_linear_geometry_signal_within_target():
    target=np.array(["a"]*8+["b"]*8); g=np.arange(16,dtype=float).reshape(-1,1)
    G=np.column_stack([g,g**2,np.sin(g),np.cos(g)])
    y=np.empty(16)
    for label in ("a","b"):
        i=np.flatnonzero(target==label); y[i]=1+G[i]@np.array([.3,-.01,.2,.1])
    r=residualize_within_target(y,G,target)
    assert np.max(np.abs(r))<1e-9

def test_equal_target_gradient_recovers_known_slope():
    target=np.array(["a"]*5+["b"]*5); x=np.tile(np.array([-2.,-1.,0.,1.,2.]),2); y=2.5*x
    stat,slopes=equal_target_gradient(y,x,target)
    assert np.isclose(stat,2.5)
    assert all(np.isclose(v,2.5) for v in slopes.values())
