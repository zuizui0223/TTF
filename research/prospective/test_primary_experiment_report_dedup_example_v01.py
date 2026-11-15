"""Known-truth physical-experiment vs publication-report counterexamples."""
import copy
import json
import sys
import unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from primary_experiment_report_dedup_example_v01 import deduplicate_cohorts, EVIDENCE

class ExperimentReportDedupTests(unittest.TestCase):
    def setUp(self):
        self.e=json.loads(EVIDENCE.read_text(encoding="utf-8"))
    def test_two_dois_one_same_adult_cohort_zero_new_independent_followup(self):
        r=deduplicate_cohorts(self.e)
        self.assertEqual(r["original_unique_article_DOIs"],2)
        self.assertEqual(r["adult_transplant_report_papers"],2)
        self.assertEqual(r["adult_transplant_underlying_cohorts"],1)
        self.assertEqual(r["independent_adult_experiments_added_by_later_2014_paper"],0)
        self.assertTrue(r["author_statement_direct_not_only_author_year_or_DOI_similarity"])
    def test_different_sowing_years_remain_distinct_from_adult(self):
        r=deduplicate_cohorts(self.e)
        self.assertEqual(r["article_by_experimental_component_records"],4)
        self.assertEqual(r["distinct_underlying_experiment_cohorts"],3)
        self.assertEqual(r["original_2007_seed_sowing_intervention_cohorts"],2)
        self.assertEqual(len(set(n["underlying_experiment_id"] for n in r["report_components"])),3)
        self.assertTrue(r["does_not_certify_three_statistically_independent_species_or_populations"])
    def test_missing_explicit_author_same_experiment_mention_abstains(self):
        self.e["evidence"]["authors_explicitly_identify_same_experiment"]=False
        with self.assertRaisesRegex(ValueError,"Author bibliography link"):
            deduplicate_cohorts(self.e)
    def test_citation_to_different_2008_paper_abstains(self):
        self.e["evidence"]["reference_doi"]="10.1007/s00442-008-1156-1"
        with self.assertRaisesRegex(ValueError,"Author bibliography link"):
            deduplicate_cohorts(self.e)
    def test_reworded_statement_not_assumed_same_underlying_experiment(self):
        self.e["evidence"]["short_source_quote"]="similar population gardens"
        with self.assertRaisesRegex(ValueError,"Author bibliography link"):
            deduplicate_cohorts(self.e)
    def test_seed_sowing_cohort_collapse_prohibited(self):
        self.e["cohort"]["separate_year2003_year2004_seed_sowing_experiments_not_merged"]=False
        with self.assertRaisesRegex(ValueError,"boundaries changed"):
            deduplicate_cohorts(self.e)
    def test_fake_plant_id_matching_claim_fails(self):
        self.e["deduplication"]["raw_individual_ID_crosspaper_join_verified"]=True
        with self.assertRaisesRegex(ValueError,"qualification drifted"):
            deduplicate_cohorts(self.e)
    def test_source_only_and_old_science_immutable(self):
        r=deduplicate_cohorts(self.e)
        self.assertFalse(r["original_individual_tag_IDs_directly_matched"])
        self.assertFalse(r["raw_fitness_or_genetic_observations_read"])
        self.assertTrue(r["previous_44of50_LGM_adaptation_source_gate_unchanged"])
        self.assertEqual(r["status"],"ONE_PHYSICAL_ADULT_EXPERIMENT_TWO_PAPERS_SEED_SOWING_COMPONENTS_DISTINCT")

if __name__=="__main__":
    unittest.main()
