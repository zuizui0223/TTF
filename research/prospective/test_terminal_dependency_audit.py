import unittest
import terminal_dependency_audit as a
class Tests(unittest.TestCase):
    def test_cross_species_study_structure(self):
        v=a.summarize(['a','a','b','c','d'],['t','t','s','s','u'],['train','train','heldout','heldout','heldout'])
        self.assertEqual(v['heldout_studies_with_multiple_species'],1)
        self.assertEqual(v['heldout_species_pairs_sharing_a_study'],1)
        self.assertEqual(v['heldout_connected_component_sizes'],[2,1])
        self.assertAlmostEqual(v['synthetic_heldout_species_mean_covariance_offdiagonal_sum'],.5)
    def test_train_test_study_overlap_is_rejected(self):
        with self.assertRaises(ValueError):
            a.summarize(['a','b'],['s','s'],['train','heldout'])
    def test_no_cross_study_species_dependence_case(self):
        v=a.summarize(['a','b','c'],['a1','b1','c1'],['train','heldout','heldout'])
        self.assertEqual(v['heldout_species_pairs_sharing_a_study'],0)
        self.assertEqual(v['heldout_connected_component_sizes'],[1,1])
        self.assertEqual(v['synthetic_heldout_species_mean_covariance_offdiagonal_sum'],0)
