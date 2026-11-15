"""Scientific erratum regression: external-site reference is Raabova 2008, NOT 2007.

No original fitness data, coordinates or biological outcomes are fetched.
Original numerical receipt remains intact and is invalidated scientifically.
"""
import json
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

class RaabovaCrosspaperErratumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.err=json.loads((HERE/"raabova_2007_2008_wrong_primary_reference_erratum_v01.json").read_text())
        cls.old=json.loads((HERE/"raabova_sudova_crosspaper_sites_contract_v01.json").read_text())
        cls.out=json.loads((HERE/"source_registry/raabova_sudova_primary_crosspaper_geography_v01.json").read_text())

    def test_original_immutable_2007_target_and_2010_source_preserved(self):
        self.assertEqual(self.old["original_transplant_study_doi"],"10.1016/j.biocon.2007.07.007")
        self.assertEqual(self.old["external_primary_source_doi"],"10.3732/ajb.1000114")
        self.assertEqual(self.err["previous"]["original_2007_target_DOI"],self.old["original_transplant_study_doi"])
        self.assertEqual(self.out["cross_paper_site_source_doi"],self.old["external_primary_source_doi"])

    def test_actual_publisher_bibliography_50_is_2008_not_2007(self):
        p=self.err["author_publication_referent"]
        self.assertEqual(p["referenced_primary_DOI"],"10.1007/s00442-008-1156-1")
        self.assertEqual(p["erroneously_assigned_2007_DOI"],self.old["original_transplant_study_doi"])
        self.assertTrue(p["wrong_referent_is_directly_verified_from_publisher_bibliography"])
        self.assertEqual(p["Sudova2010_Table1_ploidy_of_four_Aster_populations"],["2x","2x","6x","6x"])
        self.assertNotEqual(p["referenced_primary_DOI"],p["erroneously_assigned_2007_DOI"])

    def test_old_numeric_record_preserved_but_cannot_claim_accuracy_failure(self):
        p=self.err["previous"]
        self.assertEqual(p["original_numeric_sites"],4)
        self.assertEqual(p["original_numeric_1km_matches"],1)
        q=self.out["source_only_primary_published_site_correspondence"]
        self.assertEqual(q["published_sites_spatially_corroborated_uniquely_within_1km"],1)
        self.assertEqual(q["published_sites_ambiguous_or_unmatched"],3)
        self.assertTrue(p["numeric_result_preserved_but_scientific_test_invalid"])
        self.assertEqual(p["original_precommitted_distance_threshold_km_unchanged"],1.0)
        self.assertEqual(self.err["scientific_reassessment"]["previous_claim_external_validation_failed_3_of_4"],"RETRACT_AS_INVALID_COMPARISON")

    def test_wrong_paper_has_hard_stop_and_old_Bontrager_stop_remains(self):
        self.assertEqual(self.err["status"],"TERMINAL_INVALID_EXTERNAL_VALIDATION_WRONG_PRIMARY_PAPER_REFERENT")
        s=self.err["scientific_reassessment"]
        self.assertEqual(s["external_primary_georeference_validation_for_2007_study"],"NOT_EVALUATED_VALID_STUDY_MATCH_REQUIRED")
        self.assertTrue(self.err["downstream_policy"]["do_not_use_this_2008_paper_reference_as_2007_validation"])
        self.assertTrue(s["prior_Bontrager_44_of_50_gate_unchanged"])
        self.assertFalse(s["LGM_fitness_adaptation_plasticity_speciation_effect_fitted"])
        self.assertTrue(s["main_untouched"])

    def test_original_garden_proxy_counts_not_modified_by_erratum(self):
        receipt=json.loads((HERE/"source_registry/hargreaves_source_as_garden_georecovery_v01.json").read_text())
        s=receipt["summary"]
        self.assertEqual(s["calibration_true_source_comparisons"],212)
        self.assertEqual(s["unlinked_distinct_source_IDs_eligible_for_garden_derived_proxy"],76)
        self.assertEqual(s["unlinked_study_labels_with_two_or_more_pop_as_garden_recovered"],21)
        self.assertEqual(s["unlinked_design_cells_with_two_georecovered_origins"],142)

if __name__=="__main__":
    unittest.main()
