"""Guard author-confirmed same 2004 adult transplant paper lineage.

No original numerical fitness, genetic or individual-tag rows are touched.
This distinguishes experiment-level deduplication from individual-level matching.
"""
import json
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent

class AuthorConfirmed2004CohortTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current=json.loads((HERE/"raabova2007_pankova2014_author_confirmed_adult_cohort_v02.json").read_text(encoding="utf-8"))
        cls.previous=json.loads((HERE/"raabova2007_pankova2014_probable_shared_adult_cohort_v01.json").read_text(encoding="utf-8"))
        cls.geo=json.loads((HERE/"source_registry/raabova2007_pankova2014_six_site_exploratory_v01.json").read_text(encoding="utf-8"))

    def test_direct_author_confirmation_links_reference_21_to_original_2007_paper(self):
        x=self.current["evidence"]
        self.assertEqual(x["reference_number"],21)
        self.assertEqual(x["reference_doi"],self.previous["target_2007"]["doi"])
        self.assertEqual(x["followup_doi"],self.previous["followup_2014"]["doi"])
        self.assertEqual(x["short_source_quote"],"same reciprocal transplant experiments after two years")
        self.assertTrue(x["authors_explicitly_identify_same_experiment"])

    def test_two_papers_are_one_2004_adult_cohort_not_two_independent_trials(self):
        x=self.current["deduplication"]
        self.assertEqual(x["papers"],2)
        self.assertEqual(len(set(x["source_dois"])),2)
        self.assertEqual(x["shared_adult_cohort_count"],1)
        self.assertEqual(x["incremental_independent_adult_cohorts_from_followup"],0)
        self.assertTrue(x["not_an_additional_independent_experiment"])
        self.assertFalse(x["raw_individual_ID_crosspaper_join_verified"])
        self.assertIsNone(x["overlap_count_of_raw_measurement_rows"])

    def test_establishment_and_2009_garden_subset_do_not_change_2004_cohort_root(self):
        x=self.current["cohort"]
        self.assertEqual(x["field_transplant"],"October 2004")
        self.assertEqual(x["plants_per_origin_garden"],25)
        self.assertEqual(x["original_six_origin_and_garden_codes"],["K1","K2","K3","S1","S2","S3"])
        self.assertEqual(x["late_gardens_analysable_after_erosion"],["K1","K3","S1","S2"])
        self.assertTrue(x["late_individuals_marked_with_numbered_metal_tags"])
        self.assertFalse(x["individual_tags_from_early_paper_verified"])
        self.assertTrue(x["separate_year2003_year2004_seed_sowing_experiments_not_merged"])

    def test_geography_remains_exploratory_and_old_wrong_paper_invalid(self):
        x=self.current
        self.assertEqual(x["cohort"]["key"],"Raabova2007__Aster_amellus__adult_transplant__2004Oct__K1K2K3S1S2S3")
        self.assertEqual(self.geo["target_primary_paper_doi"],x["evidence"]["reference_doi"])
        self.assertEqual(self.geo["external_later_paper_doi"],x["evidence"]["followup_doi"])
        self.assertTrue(self.geo["old_invalid_Sudova_2010_vs_2007_receipt_remains_invalid"])
        self.assertEqual([z["distance_km"] for z in self.geo["site_details"]],
                         [0.799619,1.812728,1.928084,0.734489,1.311237,1.92598])
        self.assertTrue(x["deduplication"]["not_independent_GPS_verification"])

    def test_scientific_freeze_not_reopened(self):
        q=self.current["boundaries"]
        self.assertTrue(q["prior_geo_receipts_unchanged"])
        self.assertTrue(q["no_fitness_effect_or_LGM_adaptation_genetics_plasticity_speciation_fit"])
        self.assertTrue(q["PR154_44of50_stop_unchanged"])
        self.assertTrue(q["PROCEED_TTFQ_and_genetic_outcomes_unchanged"])
        self.assertTrue(q["PR155_draft"])
        self.assertTrue(q["main_untouched"])

if __name__=="__main__":
    unittest.main()
