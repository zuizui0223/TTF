import unittest
from calibrate_spatial_null import wilson95,check_receipt_hashes
class Tests(unittest.TestCase):
 def test_wilson_calibration_bounds(self):
  self.assertLess(wilson95(2,200)[1],.1)
  self.assertGreater(wilson95(190,200)[0],.8)
 def test_wilson_zero(self):
  lo,hi=wilson95(0,200)
  self.assertLess(abs(lo),1e-12)
  self.assertGreater(hi,0.)
 def test_explicit_amended_feature_sha_is_applied(self):
  receipt={"feature_matrix_raw_sha256":"NEW","coordinates_raw_sha256":"UNCHANGED"}
  check_receipt_hashes(receipt,dict(receipt))
  with self.assertRaises(ValueError):
   check_receipt_hashes(receipt,{"feature_matrix_raw_sha256":"OLD","coordinates_raw_sha256":"UNCHANGED"})
if __name__=="__main__":unittest.main()
