import hashlib
import unittest
import numpy as np
import terminal_numeric_provenance as audit

class NumericProvenanceTests(unittest.TestCase):
    def test_hash_uses_c_order_little_endian_float64(self):
        x=np.arange(12,dtype=float).reshape(3,4)
        expected=hashlib.sha256(x.astype('<f8').tobytes()).hexdigest()
        self.assertEqual(audit.array_hash(x),expected)
        self.assertEqual(audit.array_hash(np.asfortranarray(x)),expected)

    def test_identical_arrays_report_no_difference(self):
        x=np.arange(12,dtype=float).reshape(3,4)
        r=audit.numeric_difference(x,x.copy())
        self.assertIsInstance(r,dict)
        if r is not None:
            self.assertEqual(r['bitwise_different_elements'],0)
            self.assertEqual(r['different_columns'],[])

    def test_one_ulp_is_not_accepted_as_exact(self):
        x=np.ones((3,4)); y=x.copy();y[1,2]=np.nextafter(1.,2.)
        r=audit.numeric_difference(x,y)
        self.assertIsInstance(r,dict)
        if r is not None:
            self.assertEqual(r['bitwise_different_elements'],1)
            self.assertEqual(r['different_columns'],[2])
            self.assertGreater(r['max_absolute_difference'],0)
            self.assertFalse(r['exact_byte_identity'])
            self.assertNotIn('values',r)

    def test_signed_zero_difference_is_preserved(self):
        x=np.array([0.,1.]);y=np.array([-0.,1.])
        r=audit.numeric_difference(x,y)
        self.assertIsInstance(r,dict)
        if r is not None:
            self.assertEqual(r['max_absolute_difference'],0)
            self.assertEqual(r['bitwise_different_elements'],1)
            self.assertFalse(r['exact_byte_identity'])

    def test_finite_and_shape_guards(self):
        with self.assertRaises(ValueError):audit.numeric_difference(np.zeros(3),np.zeros(4))
        with self.assertRaises(ValueError):audit.numeric_difference(np.array([np.nan]),np.array([np.nan]))

    def test_matching_both_matrices_required(self):
        self.assertFalse(audit.exact_pair_matches('a'*64,'b'*64,'a'*64,'c'*64))
        self.assertTrue(audit.exact_pair_matches('a'*64,'b'*64,'a'*64,'b'*64))

if __name__=='__main__':unittest.main()
