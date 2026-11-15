"""Synthetic predictor-geometry controls; absolutely no PROCEED rate outcomes."""
import sys
import unittest
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed70_plant_history_information import (
    deduplicate_study_species_site_pairs,haversine_pairs,
    equal_species_weights,weighted_center_and_std,features_for_original_pairs,
    source_only_information,
)

class PlantHistorySourceOnlyTests(unittest.TestCase):
    def test_duplicates_are_not_pseudoreplicates(self):
        a=("Plantus alpha","reference1",(-71.,43.),(-70.,43.))
        b=("Plantus alpha","reference1",(-70.,43.),(-71.,43.))
        c=("Plantus beta","reference2",(-71.,43.),(-70.,43.))
        self.assertEqual(deduplicate_study_species_site_pairs([a,a,a,b,c]),[a,b,c])

    def test_great_circle_zero_and_dateline_wrap(self):
        z=np.array([[179.5,0.],[-70.,10.]])
        y=np.array([[-179.5,0.],[-70.,10.]])
        distance=haversine_pairs(z,y)
        self.assertGreater(distance[0],100)
        self.assertLess(distance[0],112)
        self.assertEqual(distance[1],0)

    def test_equal_species_weights_not_rows(self):
        s=["A","A","A","B","C","C"]
        roles=np.array(["train","train","train","train","heldout","heldout"])
        w=equal_species_weights(s,roles)
        np.testing.assert_allclose(w[:4],[1/6,1/6,1/6,1/2])
        np.testing.assert_allclose(w[-2:],[0,0])
        np.testing.assert_allclose(w.sum(),1.)

    def test_historical_focal_formula(self):
        a=np.array([[2,4,6,8,0,0,0,0]],float)
        b=np.array([[4,8,12,16,0,0,0,0]],float)
        z1=np.array([[-100.,40.]])
        z2=np.array([[-100.,40.]])
        H,controls=features_for_original_pairs(a,b,z1,z2,np.array([1.,2.,3.,4.]))
        self.assertEqual(H.shape,(1,))
        self.assertAlmostEqual(H[0],6.)
        self.assertEqual(controls.shape,(1,11))
        self.assertEqual(controls[0,-1],0.)

    def test_predictor_fraction_is_response_blind(self):
        rng=np.random.default_rng(20261114)
        n=120
        current1=rng.normal(size=(n,4))*[4,3,500,20]+[10,15,1000,40]
        current2=current1+rng.normal(size=(n,4))*[2,2,200,8]
        history1=current1+rng.normal(size=(n,4))*[3,3,300,10]
        history2=current2+rng.normal(size=(n,4))*[3,3,300,10]
        a=np.column_stack([history1,current1])
        b=np.column_stack([history2,current2])
        z1=np.column_stack([rng.uniform(-130,-60,n),rng.uniform(5,60,n)])
        z2=z1+rng.normal(0,0.7,(n,2))
        species=[f"Testus species{i//5}" for i in range(n)]
        study=[f"study{i//5}" for i in range(n)]
        role=np.array(["train" if i<60 else "heldout" for i in range(n)])
        report=source_only_information(a,b,z1,z2,species,study,role)
        self.assertEqual(report["train_pair_count"],60)
        self.assertEqual(report["heldout_pair_count"],60)
        self.assertEqual(report["train_species_count"],12)
        self.assertEqual(report["heldout_species_count"],12)
        self.assertEqual(report["train_baseline_design_rank"],12)
        self.assertTrue(0<=report["train_unique_historical_predictor_variance_fraction"]<=1)
        self.assertLessEqual(report["train_species_effective_residual_signal"],12)
        self.assertTrue(report["no_biological_response_used"])

    def test_reject_underdetermined_training_source(self):
        x=np.ones((20,8))
        z=np.zeros((20,2))
        with self.assertRaisesRegex(ValueError,"36 unique"):
            source_only_information(x,x,z,z,["A"]*20,["s"]*20,np.array(["train"]*20))


if __name__=="__main__":
    unittest.main()
