"""Source listing only; fail-closed when original 2007 plant tag data not indexed."""
import copy,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from raabova2007_original_data_location_preflight_v01 import (
    CONTRACT,file_links,interpret_link_listing,preflight
)

class RaabovaOriginalSourceIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    def test_two_thesis_pdfs_only_no_structured_link_not_data_absence(self):
        html='<html><a href="/bitstream/handle/20.500.11956/94265/work.pdf?sequence=1">Text práce</a><a href="/bitstream/handle/20.500.11956/94265/abstract.pdf">Abstract</a></html>'
        found=file_links(html,self.c["sources"][0]["url"],set(self.c["file_extensions"]))
        a=interpret_link_listing(found,"official thesis")
        self.assertEqual(a["linked_report_or_text_files"],2)
        self.assertEqual(a["linked_structured_file_candidate_count"],0)
        self.assertFalse(a["candidate_2007_individual_tag_table_verified"])
        self.assertIn("NOT_DATASET_ABSENCE",a["status"])
    def test_xlsx_attachment_is_candidate_only_not_matched_tag_file(self):
        html='<html><a href="files/raw_2007_field_data.xlsx">Raw supplement</a></html>'
        f=file_links(html,self.c["sources"][0]["url"],set(self.c["file_extensions"]))
        x=interpret_link_listing(f,"official")
        self.assertEqual(x["linked_structured_file_candidate_count"],1)
        self.assertFalse(x["candidate_2007_individual_tag_table_verified"])
        self.assertEqual(x["files"][0]["extension"],"xlsx")
    def test_repeated_identical_hrefs_deduplicated(self):
        html='<html><a href="/bitstream/file.pdf">a</a><a href="/bitstream/file.pdf">b</a></html>'
        f=file_links(html,self.c["sources"][0]["url"],set(self.c["file_extensions"]))
        self.assertEqual(len(f),1)
    def test_2007_thesis_link_does_not_automatically_link_individual_IDs(self):
        class Fake:
            status=200
            headers={"Content-Type":"text/html"}
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def read(self,n):return b'<html><a href="thesis.pdf">Thesis</a></html>'
            def geturl(self):return "https://dspace.cuni.cz/handle/20.500.11956/94265"
        with patch("raabova2007_original_data_location_preflight_v01.urlopen",return_value=Fake()):
            x=preflight(self.c)
        self.assertEqual(len(x["sources"]),3)
        self.assertIsNone(x["2007_2014_same_individual_tag_join"])
        self.assertFalse(x["original_2007_individual_number_table_confirmed"])
        self.assertFalse(x["source_PDFs_or_data_files_downloaded"])
    def test_403_is_access_not_proof_dataset_absent(self):
        def blocked(req,timeout=28):
            raise HTTPError(req.full_url,403,"Forbidden",None,None)
        with patch("raabova2007_original_data_location_preflight_v01.urlopen",side_effect=blocked):
            x=preflight(self.c)
        self.assertTrue(all(v["status"]=="OFFICIAL_SOURCE_LIST_ACCESS_OR_SCHEMA_NOT_EVALUABLE" for v in x["sources"]))
        self.assertTrue(x["no_hit_not_proof_of_unavailability"])
        self.assertTrue(x["prior_Bontrager_44of50_stop_unchanged"])
    def test_unapproved_repository_or_unfrozen_policy_abstains(self):
        c=copy.deepcopy(self.c)
        c["sources"][0]["url"]="https://unknown-mirror.example.com/data.csv"
        with self.assertRaisesRegex(ValueError,"host"):
            preflight(c)
        c=copy.deepcopy(self.c)
        c["status"]="OUTCOME_OPENED"
        with self.assertRaisesRegex(ValueError,"not predeclared"):
            preflight(c)
if __name__=="__main__":
    unittest.main()
