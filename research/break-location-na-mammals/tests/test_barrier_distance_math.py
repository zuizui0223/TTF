import sys,unittest
from pathlib import Path
import numpy as np
from shapely.geometry import box,LineString
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from barrier_distance_math import distance_covariates_wgs84,boundary_network,original_distance_transform

class Tests(unittest.TestCase):
 def setUp(self):
  self.region1=[box(0,0,1,1),box(1,0,2,1)]
  self.region2=[box(0,0,2,0.5),box(0,0.5,2,1)]
 def test_internal_ecotone_not_dissolved(self):
  result=distance_covariates_wgs84(
    [(1,0.25),(.5,.25)], [box(.99,.20,1.01,.30)],
    [LineString([(0,0),(0,1)])], [box(1.8,.8,1.99,.99)],
    self.region1,self.region2)
  self.assertLess(result[0,2],.005)  # <5 m projection approximation at 0.5-degree edge
  self.assertGreater(result[1,2],40)
  self.assertEqual(result[0,0],0.)
  self.assertGreater(result[1,1],40)
 def test_distance_sqrt_nonnegative(self):
  np.testing.assert_allclose(original_distance_transform(np.array([0.,4.,100.])),[0.,2.,10.])
  with self.assertRaises(ValueError):original_distance_transform([-1])
 def test_bad_location_fails(self):
  with self.assertRaises(ValueError):
   distance_covariates_wgs84([(190,20)],[box(0,0,1,1)],
      [LineString([(0,0),(1,1)])],[box(0,0,1,1)],self.region1,self.region2)
if __name__=="__main__":unittest.main()
