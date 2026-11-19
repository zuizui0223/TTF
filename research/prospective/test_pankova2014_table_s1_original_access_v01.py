import sys
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError

sys.path.insert(0,str(Path(__file__).resolve().parent))
from pankova2014_table_s1_original_access_v01 import (
    OLE_SIGNATURE,ZIP_SIGNATURE,LINKS,recognized_storage_signature,probe
)

class OriginalPankovaTableS1AccessSafetyTests(unittest.TestCase):
    def test_ole_signature(self):
        self.assertEqual(recognized_storage_signature(OLE_SIGNATURE+b"\x00"*504),
                         "OLE2_BINARY_EXCEL_CANDIDATE")
    def test_zip_signature(self):
        self.assertEqual(recognized_storage_signature(ZIP_SIGNATURE+b"\x00"*508),
                         "ZIP_XLSX_STYLE_CANDIDATE")
    def test_html_gate(self):
        self.assertEqual(recognized_storage_signature(b"<!DOCTYPE html>\n"),
                         "HTML_INTERSTITIAL_NOT_WORKBOOK")
    def test_403_never_certifies_individual_ids(self):
        def blocked(req,timeout=30):
            raise HTTPError(req.full_url,403,"Denied",None,None)
        with patch("pankova2014_table_s1_original_access_v01.urlopen",side_effect=blocked):
            out=probe()
        self.assertEqual(out["status"],"ORIGINAL_TABLE_S1_ACCESS_NOT_EVALUABLE")
        self.assertTrue(all(r["status"]=="OFFICIAL_SUPPLEMENT_ACCESS_GATED" for r in out["attempts"]))
        self.assertFalse(out["worksheet_columns_read"])
        self.assertFalse(out["individual_plant_tag_numbers_read"])
        self.assertFalse(out["crosspaper_individual_tag_ID_join_completed"])
        self.assertEqual(out["new_independent_experiment_cohort_added"],0)
    def test_valid_ole_first_512_bytes_still_not_id_confirmation(self):
        class Fake:
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def read(self,n):return OLE_SIGNATURE+b"\0"*(n-8)
            def geturl(self):return LINKS[0]
            status=206
            @property
            def headers(self):return {"Content-Type":"application/vnd.ms-excel","Content-Length":"512"}
        with patch("pankova2014_table_s1_original_access_v01.urlopen",return_value=Fake()):
            out=probe()
        self.assertEqual(out["status"],"ORIGINAL_TABLE_S1_STORAGE_ACCESS_OBSERVED_ONLY")
        self.assertEqual(len(out["attempts"]),1)
        self.assertFalse(out["individual_plant_tag_numbers_read"])
        self.assertFalse(out["supplement_row_data_read"])
        self.assertTrue(out["actual_shared_2004_adult_experiment_confirmed_by_author"])
if __name__=="__main__":
    unittest.main()
