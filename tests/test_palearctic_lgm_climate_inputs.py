import importlib.util
from pathlib import Path
import numpy as np
from shapely.geometry import box


SCRIPT=Path("scripts/prepare_palearctic_lgm_climate_inputs.py")


def load():
    spec=importlib.util.spec_from_file_location("climin",SCRIPT)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_frozen_global_centroid_origin_and_resolution():
    m=load()
    lon,lat=m.global_centres()
    assert len(lon)==3600 and len(lat)==1800
    assert np.isclose(lon[0],-179.95)
    assert np.isclose(lon[-1],179.95)
    assert np.isclose(lat[0],-89.95)
    assert np.isclose(lat[-1],89.95)


def test_trace_bio1_is_harmonized_from_kelvin_to_celsius_only():
    m=load()
    x=np.array([273.15,283.15])
    assert np.allclose(m.harmonize(x,"bio1","lgm"),[0,10])
    assert np.allclose(m.harmonize(x,"bio7","lgm"),x)
    assert np.allclose(m.harmonize(x,"bio1","current"),x)


def test_palearctic_grid_uses_fixed_global_centres():
    m=load()
    pts=m.palearctic_grid_points(box(0.0,0.0,0.2,0.2))
    got={tuple(np.round(x,2)) for x in pts}
    assert got=={(0.05,0.05),(0.05,0.15),(0.15,0.05),(0.15,0.15)}
