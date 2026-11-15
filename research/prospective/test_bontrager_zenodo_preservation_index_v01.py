"""Synthetic tests: never infer a Bontrager original archive from title similarity."""
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
from bontrager_zenodo_preservation_index_v01 import (
    TARGET_DOI, canonical_doi, exact_original_doi_in_record, original_public_files, probe
)

class BontragerPreservationSearchTests(unittest.TestCase):
    def rec(self, doi=TARGET_DOI, id=321):
        return {
            "id": id, "metadata": {"title":"Adaptation across geographic ranges",
                "related_identifiers":[{"identifier":"https://doi.org/"+doi}]},
            "files":{"entries":{"archive":{"key":"dryad.zip",
                                    "size":811204,
                                    "checksum":"md5:2de98c1bb2cd2ddffcb090eb65a74d9f"}}},
        }
    def test_exact_doi(self):
        self.assertEqual(canonical_doi("doi:10.5061/DRYAD.PNVX0K6MC/"),TARGET_DOI)
        self.assertTrue(exact_original_doi_in_record(self.rec()))
    def test_title_only_does_not_authorize(self):
        r=self.rec()
        r["metadata"]["related_identifiers"]=[]
        self.assertFalse(exact_original_doi_in_record(r))
    def test_similar_other_doi_does_not_authorize(self):
        self.assertFalse(exact_original_doi_in_record(self.rec("10.5061/dryad.pnvx0k6r0")))
    def test_read_file_list_only(self):
        files=original_public_files(self.rec())
        self.assertEqual(files[0]["filename"],"dryad.zip")
        self.assertEqual(files[0]["byte_length"],811204)
    def test_legacy_file_list(self):
        r=self.rec()
        r["files"]=[{"key":"dryad.zip", "size":811204}]
        self.assertEqual(len(original_public_files(r)),1)
    def test_exact_match_is_only_record_returned(self):
        class Fake:
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def read(self,n):
                return json.dumps({"hits":{"hits":[valid,invalid]}}).encode()
        valid=self.rec()
        invalid=self.rec("10.5061/dryad.other",322)
        with patch("bontrager_zenodo_preservation_index_v01.urlopen",return_value=Fake()):
            a=probe()
        self.assertEqual(len(a["matched_records"]),1)
        self.assertEqual(a["matched_records"][0]["zenodo_record_id"],"321")
        self.assertFalse(a["fitness_or_experiment_data_read"])
    def test_403_is_not_absence_proof(self):
        from urllib.error import HTTPError
        def block(req,timeout=35):
            raise HTTPError(req.full_url,403,"Forbidden",None,None)
        with patch("bontrager_zenodo_preservation_index_v01.urlopen",side_effect=block):
            a=probe()
        self.assertEqual(a["matched_records"],[])
        self.assertTrue(a["absence_in_index_NOT_proof_of_archive_absence"])
        self.assertFalse(a["source_study_experiment_overlap_measured"])

if __name__=="__main__":
    unittest.main()
