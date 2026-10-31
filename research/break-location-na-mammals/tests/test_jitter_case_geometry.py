import unittest
import numpy as np
import shapely
from shapely.geometry import box
from jitter_case_geometry import jitter_cases

class Tests(unittest.TestCase):
 def test_fixed_circle_accepted_in_polygon(self):
  p=box(-200000,-200000,200000,200000)
  pts=np.array([[0,0],[2000,0],[0,1000]],float)
  moved,r=jitter_cases(pts,p,np.random.default_rng(17))
  d=np.linalg.norm(moved-pts,axis=1)
  self.assertTrue((d>0).all() and (d<=50000.000001).all())
  self.assertEqual(r["successful"],3)
  moved2,_=jitter_cases(pts,p,np.random.default_rng(17))
  np.testing.assert_array_equal(moved,moved2)
 def test_narrow_polygon_cannot_produce_invalid_move(self):
  p=box(0,0,1000,1000)
  moved,r=jitter_cases(np.array([[500.,500.]]),p,np.random.default_rng(19))
  self.assertTrue(bool(shapely.covers(p,shapely.points(moved)).all()))
  self.assertGreaterEqual(r["unmoved"],0)
 def test_frozen_radius_cannot_change(self):
  with self.assertRaises(ValueError):
   jitter_cases(np.array([[0.,0.]]),box(-1,-1,1,1),np.random.default_rng(1),radius_m=55000)
if __name__=="__main__":unittest.main()
