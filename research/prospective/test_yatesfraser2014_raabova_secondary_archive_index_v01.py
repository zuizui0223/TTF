"""Synthetic DOI-exact secondary archive preflight tests, no raw fitness rows."""
import copy,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from yatesfraser2014_raabova_secondary_archive_index_v01 import (
    CONTRACT,doi,exact_original_doi,safe_files,zenodo_meta,dryad_meta,probe
)

class SecondaryRaabovaDryadMetadataTests(unittest.TestCase):
    def setUp(self):
        self.c=json.loads(CONTRACT.read_text())
    def rec(self,ident="10.5061/dryad.81k4v",rid=123):
        return {"id":rid,"metadata":{"title":"Does source population size affect performance",
                "related_identifiers":[{"identifier":"https://doi.org/"+ident}]},
                "files":{"entries":{"x":{"key":"survival_table.csv","size":555,
                    "checksum":"md5:abc"}}}}
    def test_structured_exact_doi_not_title_similarity(self):
        a=self.rec()
        self.assertTrue(exact_original_doi(a,"10.5061/dryad.81k4v"))
        a["metadata"]["related_identifiers"]=[]
        self.assertFalse(exact_original_doi(a,"10.5061/dryad.81k4v"))
        self.assertEqual(doi("doi:10.5061/DRYAD.81K4V/"),"10.5061/dryad.81k4v")
    def test_similar_but_different_DOI_must_not_qualify(self):
        self.assertFalse(exact_original_doi(self.rec("10.5061/dryad.81k4x"),
                                            "10.5061/dryad.81k4v"))
    def test_no_data_download_from_file_metadata(self):
        files=safe_files(self.rec())
        self.assertEqual(files[0]["filename"],"survival_table.csv")
        self.assertEqual(files[0]["size_bytes"],555)
        self.assertEqual(set(files[0]),{"filename","size_bytes","published_checksum"})
    def test_access_error_is_not_published_data_absence(self):
        def fail(req,timeout=30):raise HTTPError(req.full_url,403,"Forbidden",None,None)
        with patch("yatesfraser2014_raabova_secondary_archive_index_v01.urlopen",side_effect=fail):
            q=probe(self.c)
        self.assertEqual(q["status"],"SOURCE_METADATA_NO_CONFIRMED_INDEX_NOT_PROOF_OF_ABSENCE")
        self.assertEqual(q["dryad_original_metadata"]["http_status"],403)
        self.assertFalse(q["raw_data_files_or_spreadsheets_downloaded"])
        self.assertIsNone(q["2007_vs_2014_individual_overlap_count"])
    def test_zenodo_duplicate_matches_deduped_by_record_id(self):
        good=self.rec()
        wrong=self.rec("10.5061/dryad.invalid",124)
        body=json.dumps({"hits":{"hits":[good,wrong]}}).encode("utf-8")
        class Fake:
            def __enter__(self):return self
            def __exit__(self,*a):return False
            def read(self,n):return body
        with patch("yatesfraser2014_raabova_secondary_archive_index_v01.urlopen",return_value=Fake()):
            results=zenodo_meta(self.c)
        self.assertEqual(len(results["matched_exact_original_Dryad_doi_preservation_records"]),1)
        self.assertEqual(results["matched_exact_original_Dryad_doi_preservation_records"][0]["id"],123)
    def test_wrong_dryad_dataset_rejected(self):
        body=json.dumps({"identifier":"doi:10.5061/dryad.other","title":"same title"}).encode()
        class Fake:
            def __enter__(self):return self
            def __exit__(self,*a):return False
            def read(self,n):return body
        with patch("yatesfraser2014_raabova_secondary_archive_index_v01.urlopen",return_value=Fake()):
            result=dryad_meta(self.c)
        self.assertEqual(result["status"],"OFFICIAL_DRYAD_API_METADATA_ACCESS_NOT_EVALUABLE")
        self.assertEqual(result["error_type"],"ValueError")
    def test_original_meta_not_authorized_for_secondary_fitness(self):
        self.c["sources"]["data_file_download_authorized"]=True
        with self.assertRaisesRegex(ValueError,"not authorized"):
            probe(self.c)
if __name__=="__main__":
    unittest.main()
