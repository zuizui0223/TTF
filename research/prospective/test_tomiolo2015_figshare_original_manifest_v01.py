"""Synthetic exact DOI and no-original-data-opening Figshare manifest tests."""
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from tomiolo2015_figshare_original_manifest_v01 import (
    CONTRACT,doi_norm,files_from_article,valid_collection,run
)

class TomioloOfficialManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    def test_exact_stable_collection_DOI_and_id_qualified(self):
        original={"id":3307653,"doi":"https://doi.org/10.6084/m9.figshare.c.3307653.v1"}
        self.assertEqual(valid_collection(original,self.c),
                         "10.6084/m9.figshare.c.3307653.v1")
        self.assertEqual(doi_norm("doi:10.6084/M9.FIGSHARE.C.3307653.V1/"),
                         "10.6084/m9.figshare.c.3307653.v1")
    def test_wrong_id_and_title_only_never_qualified(self):
        with self.assertRaisesRegex(ValueError,"collection id"):
            valid_collection({"id":3307654,"doi":self.c["target"]["figshare_collection_doi"]},self.c)
        with self.assertRaisesRegex(ValueError,"DOI mismatch"):
            valid_collection({"id":3307653,"title":"same title"},self.c)
    def test_article_files_return_only_original_file_metadata(self):
        x={"files":[{"id":12,"name":"data.csv","size":440,
                     "computed_md5":"1234","download_url":"https://other.invalid/data"}]}
        result=files_from_article(x,self.c["figshare_source_manifest_policy"])
        self.assertEqual(result[0]["original_file_name"],"data.csv")
        self.assertEqual(result[0]["published_size_bytes"],440)
        self.assertNotIn("download_url",result[0])
        self.assertTrue(result[0]["no_source_file_body_opened"])
    def test_403_never_establishes_data_absence(self):
        def blocked(url,policy):
            raise HTTPError(url,403,"Forbidden",None,None)
        with patch("tomiolo2015_figshare_original_manifest_v01.get_official_metadata",
                   side_effect=blocked):
            r=run()
        self.assertEqual(r["status"],"ORIGINAL_FIGSHARE_COLLECTION_METADATA_NOT_EVALUABLE")
        self.assertFalse(r["all_data_file_bodies_or_biological_outcome_values_opened"])
        self.assertFalse(r["source_experiment_primary_independence_certified"])
    def test_fake_exact_record_but_missing_article_manifest_not_independence(self):
        def stub(url,policy):
            if url.endswith("/3307653"):
                return {"id":3307653,"doi":self.c["target"]["figshare_collection_doi"],"title":"Original"}
            if url.endswith("/3307653/articles"):
                return [{"id":7654321}]
            if url.endswith("/7654321"):
                return {"id":7654321,"doi":"10.6084/m9.figshare.7654321",
                        "title":"Original data article","files":[{"name":"table.csv","id":78,"size":100}]}
            raise AssertionError(url)
        with patch("tomiolo2015_figshare_original_manifest_v01.get_official_metadata",side_effect=stub):
            r=run()
        self.assertEqual(r["status"],"PUBLIC_FIGSHARE_SOURCE_MANIFEST_LISTED_NO_DATA_FILES")
        self.assertEqual(len(r["figshare_public_article_manifest"]),1)
        self.assertFalse(r["source_experiment_primary_independence_certified"])
        self.assertFalse(r["original_per_source_population_GPS_certified"])
        self.assertFalse(r["all_data_file_bodies_or_biological_outcome_values_opened"])
    def test_frozen_6_source_IDs_are_two_geographic_sites_not_six_climates(self):
        x=self.c["already_frozen_source_only_Hargreaves"]
        self.assertEqual(x["origin_population_ID_keys_across_all_species"],6)
        self.assertEqual(x["geographically_unique_garden_derived_coordinates"],2)
        self.assertEqual(x["comparable_same_garden_source_ID_cells"],12)
        self.assertFalse(self.c["qualification_gates"]["unique_original_source_population_GPS_verified"])
        self.assertTrue(self.c["qualification_gates"]["previous_Bontrager_44of50_gate_unchanged"])

if __name__=="__main__":
    unittest.main()
