"""No-data tests for conservative exact-DOI official preservation index."""
import io
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parent))
from hargreaves_zenodo_preservation_index_v01 import (
  TARGET_DOI, TARGET_FILENAME, doi_norm, is_original_dryad_doimatch,
  metadata_inventory, probe,
)

class ZenodoDryadPreservationIndexTests(unittest.TestCase):
    def rec(self, doi=TARGET_DOI, file=TARGET_FILENAME):
        return {"id":1234,"metadata":{"related_identifiers":[
            {"identifier":f"https://doi.org/{doi}","relation":"isSupplementTo"}]},
            "title":"Local adaptation to biotic interactions"},
            "files":{"entries":{"f":{"key":file}}}}

    def test_doi_canonicalization(self):
        self.assertEqual(doi_norm("https://doi.org/10.5061/DRYAD.0VT4B8GV1/"), TARGET_DOI)

    def test_exact_original_dryad_match_from_related_identifier(self):
        self.assertTrue(is_original_dryad_doimatch(self.rec()))
        x=metadata_inventory(self.rec())
        self.assertTrue(x["target_original_csv_filename_in_public_metadata"])
        self.assertFalse(x["data_values_or_headers_read"])

    def test_title_only_never_qualifies(self):
        r=self.rec()
        r["metadata"]["related_identifiers"]=[]
        with self.assertRaisesRegex(ValueError,"No explicit"):
            metadata_inventory(r)
        self.assertFalse(is_original_dryad_doimatch(r))

    def test_other_dryad_DOI_is_not_the_target(self):
        self.assertFalse(is_original_dryad_doimatch(self.rec("10.5061/dryad.pnvx0k6mc")))

    def test_exact_doi_but_missing_file_is_source_only_candidate_not_data(self):
        r=self.rec(file="other.csv")
        x=metadata_inventory(r)
        self.assertTrue(x["source_original_doi_explicitly_matched"])
        self.assertFalse(x["target_original_csv_filename_in_public_metadata"])

    def test_legacy_list_file_shape(self):
        r=self.rec()
        r["files"]=[{"key":TARGET_FILENAME}]
        self.assertTrue(metadata_inventory(r)["target_original_csv_filename_in_public_metadata"])

    def test_404_or_403_is_not_archival_absence(self):
        from urllib.error import HTTPError
        def blocked(req, timeout=30):
            raise HTTPError(req.full_url,403,"Forbidden",None,None)
        with patch("hargreaves_zenodo_preservation_index_v01.urlopen",side_effect=blocked):
            x=probe()
        self.assertEqual(len(x["verified_exact_doi_candidates"]),0)
        self.assertTrue(x["no_candidate_is_NOT_proof_of_absence"])
        self.assertTrue(all(a["status"]=="ACCESS_ERROR" for a in x["attempts"]))
        self.assertFalse(x["Hargreaves_fitness_outcome_opened"])

    def test_search_payload_deduplicates_same_exact_doi_record(self):
        class Fake:
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def read(self, cap):return json.dumps({"hits":{"hits":[self_rec]}}).encode()
        self_rec=self.rec()
        with patch("hargreaves_zenodo_preservation_index_v01.urlopen",return_value=Fake()):
            x=probe()
        self.assertEqual(len(x["verified_exact_doi_candidates"]),1)
        self.assertEqual(x["source_index_status"],"EXACT_DOI_PUBLIC_PRESERVATION_RECORD_LOCATED")
        self.assertFalse(x["source_study_overlap_measured"])

if __name__=="__main__":unittest.main()
