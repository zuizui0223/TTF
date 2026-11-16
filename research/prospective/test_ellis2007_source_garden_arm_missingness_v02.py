"""Known-truth Ellis three-basin source-only missing-arm classification tests."""
import copy
import json
import sys
import unittest
from collections import defaultdict
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from ellis2007_source_garden_arm_missingness_v02 import (
    CONTRACT,arm_census,classify,CATEGORIES
)
class EllisMissingOriginalTransplantArmsTests(unittest.TestCase):
    def setUp(self):
        self.c=copy.deepcopy(json.loads(CONTRACT.read_text(encoding="utf-8")))
        self.study=self.c["previous_sha_verified_geometry"]["study_label"]
        self.taxon=self.c["previous_sha_verified_geometry"]["species"]
    def fixture(self,absent="same_basin_other_quartz"):
        positions={}
        for i,basin in enumerate(("Ariz","Sout","QK")):
            for h,j in (("dense",0),("sparse",1)):
                name=basin+"_"+h
                positions[name]=(18.41+i*0.05,-31.1+j*0.005)
        sourcekeys={(self.study,self.taxon,s) for s in positions}
        gardens=defaultdict(set)
        for k in sourcekeys:
            gardens[k].add(positions[k[2]])
        def category(origin,garden):
            a=origin.split("_");b=garden.split("_")
            return CATEGORIES[0] if a==b else CATEGORIES[1] if a[0]==b[0] else CATEGORIES[2] if a[1]==b[1] else CATEGORIES[3]
        groups={}
        for garden in positions:
            sources={s for s in positions if category(s,garden)!=absent and category(s,garden)!=CATEGORIES[3]}
            for outcome_type in ("survival","size","reproduction"):
                groups[self.study,self.taxon,"adult2000",garden,outcome_type]=sources
        return {"origin_keys":sourcekeys,"gardens":gardens,"source_groups":groups}
    def test_original_24_total_18_secondary_missing_6_same_basin_other_quartz(self):
        v=arm_census(self.fixture(),self.c)
        self.assertEqual(v["original_primary_expected_directed_source_garden_pairs"],24)
        self.assertEqual(v["observed_Hargreaves_distinct_directed_source_garden_pairs"],18)
        self.assertEqual(v["source_garden_pairs_absent_from_Hargreaves_extract"],6)
        counts={x["class"]:x for x in v["published_origin_site_pairs_by_class"]}
        self.assertEqual(counts[CATEGORIES[0]]["original_Hargreaves_extract_observed"],6)
        self.assertEqual(counts[CATEGORIES[1]]["not_represented_in_secondary_original"],6)
        self.assertEqual(counts[CATEGORIES[2]]["original_Hargreaves_extract_observed"],12)
        self.assertEqual(counts[CATEGORIES[3]]["original_Hargreaves_extract_observed"],0)
        self.assertTrue(v["source_pair_absence_is_not_original_2000_intervention_absence"])
        self.assertTrue(v["no_fitness_growth_survival_reproduction_or_AFLP_response_rows_read"])
    def test_alternative_curation_missing_six_local_pairs_not_force_known_cause(self):
        v=arm_census(self.fixture(absent="local_same_source_and_garden"),self.c)
        d={r["class"]:r for r in v["published_origin_site_pairs_by_class"]}
        self.assertEqual(d[CATEGORIES[0]]["not_represented_in_secondary_original"],6)
        self.assertEqual(d[CATEGORIES[1]]["original_Hargreaves_extract_observed"],6)
        self.assertFalse(v["coordinate_discrepancy_cause_certified"])
    def test_illegal_crossbasin_different_quartz_does_not_become_expected(self):
        m=self.fixture()
        k=(self.study,self.taxon,"adult2000","Ariz_dense","survival")
        m["source_groups"][k]=set(m["source_groups"][k])|{"Sout_sparse"}
        with self.assertRaisesRegex(ValueError,"not compatible"):
            arm_census(m,self.c)
    def test_old_experimental_raw_source_and_garden_labels_must_be_classifiable(self):
        m=self.fixture()
        target=(self.study,self.taxon,"adult2000","Ariz_dense","survival")
        m["source_groups"][target].add("unknown_population")
        with self.assertRaisesRegex(ValueError,"unrecognized original source"):
            arm_census(m,self.c)
    def test_original_meta_subset_underrepresented_is_not_full_2000_design(self):
        v=arm_census(self.fixture(),self.c)
        self.assertTrue(v["source_pair_absence_may_be_downstream_meta_analysis_extraction"])
        self.assertEqual(v["published_source_population_IDs_per_target_garden"],4)
        self.assertEqual(v["observed_secondary_source_IDs_per_target_garden"],3)
        self.assertEqual(v["source_garden_pairs_absent_from_Hargreaves_extract"],6)
        self.assertFalse(v["no_LGM_21ka_0BP_climate_effect_estimated"] is False)
        self.assertTrue(v["old_Bontrager_44of50_terminal_stop_unchanged"])
    def test_published_distance_is_descriptive_not_threshold_fitted(self):
        m=self.fixture()
        v=arm_census(m,self.c)
        self.assertEqual(set(v["within_basin_dense_sparse_pair_distance_km"]),{"Ariz","Sout","QK"})
        self.assertEqual(v["publisher_within_basin_approximate_distance_claim_km"],[0.4,1.2])
        self.assertFalse(v["coordinate_discrepancy_cause_certified"])
    def test_primary_design_finite_precommitted_guard(self):
        self.c["predeclared_primary_article_transplant_topology"]["expected_total_pair_edges"]=18
        with self.assertRaisesRegex(ValueError,"topology changed"):
            arm_census(self.fixture(),self.c)
    def test_no_fuzzy_label_assignment_or_binomial_changes(self):
        m=self.fixture()
        k=next(iter(m["origin_keys"]))
        m["origin_keys"].remove(k)
        m["origin_keys"].add((self.study,self.taxon,"unknown"))
        with self.assertRaisesRegex(ValueError,"classifiable"):
            arm_census(m,self.c)
if __name__=="__main__":
    unittest.main()
