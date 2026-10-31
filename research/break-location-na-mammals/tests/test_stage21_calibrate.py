import unittest
from calibrate_spatial_null import wilson95
class Tests(unittest.TestCase):
 def test_wilson_calibration_bounds(self):
  self.assertLess(wilson95(2,200)[1],.1)
  self.assertGreater(wilson95(190,200)[0],.8)
 def test_wilson_zero(self):
  lo,hi=wilson95(0,200)
  self.assertLess(abs(lo),1e-12)
  self.assertGreater(hi,0.)
if __name__=="__main__":unittest.main()
