"""Synthetic tests: never open transplant fitness outcomes or study coordinates."""
import unittest
from dryad_metadata_v01 import (parse_one_dataset,parse_files_page,
    safe_dryad_api_url,select_version_link,project_listing)

class TestDryadCatalog(unittest.TestCase):
    def test_official_metadata_api_host_only(self):
        correct="https://datadryad.org/api/v2/versions/44/files"
        self.assertEqual(safe_dryad_api_url(correct),correct)
        for bad in ("https://evil.example.com/api/v2/versions/44/files",
                    "http://datadryad.org/api/v2/versions/44/files",
                    "https://datadryad.org/dataset/doi:dryad.x"):
            with self.assertRaises(ValueError):safe_dryad_api_url(bad)

    def test_original_doi_and_version_identity(self):
        d={"identifier":"doi:10.5061/dryad.x","title":"transplant study",
           "_links":{"stash:version":{"href":"/api/v2/versions/32"}}}
        self.assertEqual(select_version_link(d),"https://datadryad.org/api/v2/versions/32")
        p=parse_one_dataset(d,"10.5061/dryad.x","primary_multispecies")
        self.assertEqual(p["current_version_id"],"32")
        self.assertFalse(p["phenotype_rows_inspected"])
        with self.assertRaises(ValueError):
            parse_one_dataset(d,"10.5061/dryad.other","test")

    def test_HAL_files_only_not_data_rows(self):
        d={"_embedded":{"stash:files":[{
            "path":"Data.csv","size":1024,"mimeType":"text/csv",
            "digest":"f"*64,"digestType":"sha-256",
            "_links":{"stash:download":{"href":"/api/v2/files/17/download"}}}]},
            "_links":{"next":{"href":"/api/v2/versions/32/files?page=2"}}}
        files,more=parse_files_page(d)
        self.assertEqual(files[0]["sha256_catalog"],"f"*64)
        self.assertEqual(files[0]["download_href"],"/api/v2/files/17/download")
        self.assertEqual(more,"https://datadryad.org/api/v2/versions/32/files?page=2")
        self.assertFalse(files[0]["data_records_read"])

    def test_missing_target_not_silently_replaced(self):
        c={"id":"first","doi":"10.5061/dryad.x","role":"primary",
           "expected_file_name":"transplant_data.csv",
           "why":"source-only","known_overlap":["other"]}
        d={"identifier":"doi:10.5061/dryad.x",
           "_links":{"stash:version":{"href":"/api/v2/versions/3"}}}
        p=project_listing(c,d,[{"path":"other.csv","size_bytes":200}])
        self.assertFalse(p["expected_file_present"])
        self.assertEqual(p["source_catalog_status"],"EXPECTED_SOURCE_FILE_NOT_FOUND")
        self.assertEqual(p["independent_data_status"],"UNKNOWN_POSSIBLE_OVERLAP")

if __name__=="__main__":unittest.main()
