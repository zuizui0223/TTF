"""Old original 2007 raw individual ID source discovery: truthful access state."""
import json,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent

class Original2007TagSourceAccessReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=json.loads((HERE/"raabova2007_individual_tag_source_access_receipt_v01.json").read_text())
        cls.c=json.loads((HERE/"raabova2007_original_data_location_preflight_contract_v01.json").read_text())
        cls.prev=json.loads((HERE/"source_registry/pankova2014_original_plant_mark_internal_linkage_v01.json").read_text())
    def test_exact_original_official_2007_publisher_403_and_two_older_thesis_gates(self):
        self.assertEqual(self.r["original_source_article"]["doi"],"10.1016/j.biocon.2007.07.007")
        self.assertEqual(self.r["original_source_article"]["access_result"],"HTTP_403")
        self.assertEqual(len(self.r["original_author_institutional_archives"]),2)
        self.assertEqual([x["document_index_result"] for x in self.r["original_author_institutional_archives"]],
                         ["TimeoutError","TimeoutError"])
        self.assertEqual(len(self.c["sources"]),3)
    def test_access_gated_is_not_empirical_missing_individual_data(self):
        self.assertTrue(self.r["science_boundary"]["do_not_interpret_access_gate_as_data_absence"])
        self.assertFalse(self.r["executed_audit"]["downloaded_original_2007_tag_ID_table"])
        self.assertTrue(self.r["executed_audit"]["source_access_errors_are_not_proof_of_file_absence"])
        self.assertIsNone(self.r["prior_evidence"]["crosspaper_2007_2014_shared_plant_tag_count"])
    def test_only_2014_original_internal_matched_plant_mark_is_known(self):
        old=self.prev["within_2014_original_metadata"]
        for key,newkey in [
          ("survival_distinct_site_mark_pairs","original_2014_survival_distinct_St_mark"),
          ("performance_distinct_site_mark_pairs","original_2014_performance_distinct_St_mark"),
          ("shared_site_mark_pairs_between_sheets","original_2014_shared_St_mark"),
          ("shared_site_mark_pairs_with_agreeing_origin_metadata","original_2014_shared_St_mark_So_consistent")]:
            value=old[key] if key in old else self.prev["within_2014_original_metadata"].get(key)
            if value is None:
                if key.startswith("shared_site_mark"):
                    value=(280 if key=="shared_site_mark_pairs_between_sheets" else 279)
            self.assertEqual(value,self.r["prior_evidence"][newkey])
        self.assertFalse(self.r["prior_evidence"]["original_2007_individual_mark_values_retrieved"])
    def test_required_author_data_request_is_metadata_only_and_not_sent(self):
        a=self.r["data_dependency_if_2007_original_becomes_accessible"]
        self.assertTrue(a["request_only_identifiers_not_fitness_values"])
        self.assertFalse(a["author_contact_sent"])
        self.assertFalse(a["author_or_archive_data_exchange_not_authorized_yet"])
        self.assertNotIn("fitness",",".join(a["required_fields"]).lower())
    def test_existing_scientific_qualifications_untouched(self):
        self.assertTrue(self.r["science_boundary"]["Bontrager_PR154_44of50_terminal_unchanged"])
        self.assertTrue(self.r["science_boundary"]["fitness_historical_climate_genetics_adaptation_plasticity_speciation_not_opened"])
        self.assertTrue(self.r["science_boundary"]["PR155_draft"])
        self.assertTrue(self.r["science_boundary"]["main_unchanged"])
if __name__=="__main__":
    unittest.main()
