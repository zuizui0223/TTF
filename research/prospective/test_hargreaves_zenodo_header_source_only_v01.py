"""Synthetic safety tests for exact-DOI official CSV header-only access."""
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hargreaves_zenodo_header_source_only_v01 import (
    RECORD_ID, TARGET_DOI, TARGET_FILENAME, original_file_info, parse_header, probe,
)

class HargreavesZenodoHeaderTests(unittest.TestCase):
    def record(self, url=None):
        if url is None:
            url = f"https://zenodo.org/api/records/{RECORD_ID}/files/{TARGET_FILENAME}/content"
        return {
            "id": RECORD_ID,
            "metadata": {"related_identifiers": [{"identifier": TARGET_DOI}]},
            "files": {"entries": {TARGET_FILENAME: {
                "key": TARGET_FILENAME, "size": 4_380_000,
                "links": {"content": url}, "checksum": "md5:placeholder"}}},
        }

    def test_accept_exact_metadata_and_original_file(self):
        d=original_file_info(self.record())
        self.assertIn("zenodo.org/api/records",d["content_url"])

    def test_reject_matching_title_without_correct_original_doi(self):
        r=self.record()
        r["metadata"]["related_identifiers"]=[]
        with self.assertRaisesRegex(ValueError, "both be exact"):
            original_file_info(r)

    def test_reject_nonofficial_external_redirect_source(self):
        with self.assertRaisesRegex(ValueError,"not verified official"):
            original_file_info(self.record(url="https://other.example.org/csv"))

    def test_reject_zenodo_other_record(self):
        with self.assertRaisesRegex(ValueError,"not verified official"):
            original_file_info(self.record(url=f"https://zenodo.org/api/records/{RECORD_ID+1}/files/{TARGET_FILENAME}/content"))

    def test_first_line_only_required(self):
        self.assertEqual(parse_header(b"study,taxon,source,site,fitness\n"),
                         ["study","taxon","source","site","fitness"])
        with self.assertRaises(ValueError):
            parse_header(b"study,taxon,source,site,fitness")

    def test_success_never_reads_data_rows(self):
        rec=self.record()
        class FakeRecordResponse:
            def __enter__(self):return self
            def __exit__(self,*a):return False
            def read(self,n):return json.dumps(rec).encode()
        class FakeCSVResponse:
            def __enter__(self):return self
            def __exit__(self,*a):return False
            def readline(self,n):return b"study,taxon,source,site,fitness\n"
            def read(self,*a):raise AssertionError("Original response rows must not be read")
        with patch("hargreaves_zenodo_header_source_only_v01.urlopen",side_effect=[FakeRecordResponse(),FakeCSVResponse()]):
            res=probe()
        self.assertEqual(res["status"],"OFFICIAL_PRESERVED_ORIGINAL_CSV_HEADER_OBTAINED_SOURCE_ONLY")
        self.assertEqual(res["original_csv_data_rows_read"],0)
        self.assertEqual(res["source_csv_column_count"],5)
        self.assertFalse(res["source_study_overlap_measured"])

    def test_403_does_not_infer_no_overlap(self):
        from urllib.error import HTTPError
        def reject(req,timeout=30):
            raise HTTPError(req.full_url,403,"Forbidden",None,None)
        with patch("hargreaves_zenodo_header_source_only_v01.urlopen",side_effect=reject):
            res=probe()
        self.assertEqual(res["status"],"OFFICIAL_ZENODO_PRESERVED_CSV_HEADER_ACCESS_NOT_EVALUABLE")
        self.assertEqual(res["independent_new_studies_certified"],0)
        self.assertEqual(res["fitness_rows_read"],0)

if __name__=="__main__":
    unittest.main()
