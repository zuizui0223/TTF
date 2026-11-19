"""Source-only pinned XLS signature and individual metadata linkage receipt guards."""
import json
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
RECEIPT=HERE/"source_registry/pankova2014_original_plant_mark_internal_linkage_v01.json"

class OriginalPLOSPlantMarkReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=json.loads(RECEIPT.read_text(encoding="utf-8"))
        cls.c=json.loads((HERE/"pankova2014_s1_original_plant_mark_join_contract_v01.json").read_text(encoding="utf-8"))
    def test_source_original_and_artifact_hashes(self):
        a=self.r["source"]
        self.assertEqual(a["original_XLS_sha256"],self.c["original"]["sha256"])
        self.assertEqual(a["original_full_xls_bytes"],147456)
        self.assertEqual(a["original_codebook_artifact_sha256"],
                         "4a271585e9deb20f3d919e9f6a34a40a11ba1299b4bc90d1c79282a5678b8b2c")
        self.assertEqual(a["original_mark_St_So_link_artifact_sha256"],
                         "7f30d8c34fc317090782cb6f3696e6166c4ad99746d49c315ba7d4e8b0fd1d13")
    def test_publisher_mark_definition_overrides_inconclusive_header_heuristic(self):
        self.assertEqual(self.r["source_published_variable_definitions"]["mark"],
                         "The number of the plant in the field")
        self.assertEqual(self.r["source_published_variable_definitions"]["nohelp"],"Ordering Numbers")
        for sheet in ("survival","performance"):
            self.assertIn("mark",self.r["metadata_only_original_worksheet_columns"][sheet])
        self.assertEqual(self.r["precommitted_source_ID_allowlist"],["mark","St","So"])
    def test_original_composite_plant_keys_unique_in_both_sheets(self):
        v=self.r["within_2014_original_metadata"]
        self.assertEqual(v["survival_distinct_site_mark_pairs"],547)
        self.assertEqual(v["performance_distinct_site_mark_pairs"],281)
        self.assertEqual(v["duplicates_site_mark_within_survival"],0)
        self.assertEqual(v["duplicates_site_mark_within_performance"],0)
        self.assertEqual(v["target_gardens_survival"],4)
        self.assertEqual(v["target_gardens_performance"],4)
        self.assertLess(v["distinct_mark_numbers_unconditional_survival"],547)
        self.assertLess(v["distinct_mark_numbers_unconditional_performance"],281)
    def test_no_silent_metadata_conflict_or_2014_performance_only_key_repair(self):
        v=self.r["within_2014_original_metadata"]
        self.assertEqual(v["shared_site_mark_pairs_between_sheets"],280)
        self.assertEqual(v["shared_site_mark_pairs_with_origin_metadata_conflict"],1)
        self.assertEqual(v["shared_site_mark_pairs_with_agreeing_origin_metadata"],279)
        self.assertEqual(v["performance_site_mark_pairs_absent_in_survival"],1)
        self.assertEqual(v["survival_site_mark_pairs_absent_in_performance"],267)
        self.assertTrue(self.r["interpretation"]["do_not_silently_resolve_origin_conflict"])
        self.assertTrue(self.r["interpretation"]["do_not_silently_assign_unmatched_performance_key"])
    def test_2007_plant_tag_ID_join_and_fitness_still_closed(self):
        a=self.r["interpretation"]
        self.assertIsNone(a["2007_2014_shared_individual_ID_count"])
        self.assertFalse(a["2007_raw_individual_tag_file_located"])
        self.assertTrue(a["no_new_fitness_genotype_history_adaptation_plasticity_speciation_effect"])
        self.assertTrue(a["original_PR154_44of50_source_stop_unchanged"])
        self.assertTrue(self.r["within_2014_original_metadata"]["biological_survival_biomass_mycorrhizal_fitness_measurements_read"] is False)
if __name__=="__main__":
    unittest.main()
