import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
try:
 import audit_lgm_external_snapshot as a
except ImportError:
 a=None

def test_module_exists(): assert a is not None, 'snapshot audit not implemented'
def test_independent_spherical_distance():
 assert a is not None
 assert a.nearest_haversine(np.array([[0,179.9],[0,0]]),np.array([[0,-179.9]]))[0]==pytest.approx(22.239016,abs=1e-5)
def point(lat,near=100,unc=100):
 return dict(species='A',gbif_id=lat,latitude=lat,longitude=0,nearest_evaluation_km=near,uncertainty_m=unc)
def test_near_reference_rejected():
 assert a is not None
 with pytest.raises(ValueError,match='separation'):
  a.audit_points([point(.1)],{'A':np.array([[0,0]])})
def test_bad_stored_distance_rejected():
 assert a is not None
 with pytest.raises(ValueError,match='stored'):
  a.audit_points([point(1,near=150)],{'A':np.array([[0,0]])})
def test_thinning_violation_rejected():
 assert a is not None
 rows=[point(1,111.1950802335),point(1.01,112.3070310359)]
 with pytest.raises(ValueError,match='thinning'):
  a.audit_points(rows,{'A':np.array([[0,0]])})
def test_valid_points_independent_audit():
 assert a is not None
 rows=[point(1,111.1950802335),point(2,222.3901604671)]
 q=a.audit_points(rows,{'A':np.array([[0,0]])})
 assert q['retained_points']==2 and q['species_with_retained_points']==1
 assert q['minimum_verified_clearance_km']>111
