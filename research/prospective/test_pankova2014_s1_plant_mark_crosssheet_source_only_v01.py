"""Pure synthetic BIFF ID metadata examples, no published biological responses."""
import json
import math
import struct
import sys
import unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from pankova2014_s1_plant_mark_crosssheet_source_only_v01 import (
    CONTRACT,per_sheet_id_metadata,read_selected_id_cells,rk_number,scalar_id
)
def rec(kind,payload):
    return struct.pack("<HH",kind,len(payload))+payload
def number(row,col,value):
    return rec(0x0203,struct.pack("<HHHd",row,col,0,float(value)))
def as_record(rows):
    return rec(0x0809,b"\0"*8)+rows+rec(0x000A,b"")

class OriginalPlantMarkMetadataOnlyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    def test_numeric_plant_mark_and_target_site_only(self):
        cols={"mark":1,"St":5,"So":6}
        stream=as_record(
            number(1,1,11)+number(1,5,2)+number(1,6,3)+
            number(1,7,99999.)+number(2,1,12)+number(2,5,2)+number(2,6,3))
        rows,refs=read_selected_id_cells(stream,0,cols,2000)
        self.assertEqual(len(rows),2)
        self.assertEqual(rows[1],{"mark":11.,"St":2.,"So":3.})
        self.assertFalse(refs)
    def test_unique_site_plus_mark_join_not_mark_alone(self):
        r=per_sheet_id_metadata({
            "survival":{1:{"St":"K1","So":"S1","mark":1},2:{"St":"K3","So":"S1","mark":1}},
            "performance":{1:{"St":"K1","So":"S1","mark":1}}},self.c)
        self.assertEqual(r["per_sheet"]["survival"]["unique_mark_text_or_number_tokens"],1)
        self.assertEqual(r["per_sheet"]["survival"]["distinct_site_by_mark_keys"],2)
        self.assertEqual(r["within_2014_crosssheet"]["shared_site_by_mark_identifiers"],1)
        self.assertEqual(r["within_2014_crosssheet"]["survival_IDs_not_in_performance"],1)
    def test_duplicated_site_mark_fails_uniqueness(self):
        r=per_sheet_id_metadata({
            "survival":{1:{"St":"K1","mark":1,"So":"S1"},2:{"St":"K1","mark":1,"So":"S1"}},
            "performance":{1:{"St":"K1","mark":1,"So":"S1"}}},self.c)
        self.assertEqual(r["per_sheet"]["survival"]["duplicate_site_by_mark_row_count"],1)
        self.assertFalse(r["within_2014_crosssheet"]["same_site_mark_is_unique_in_both_sheets"])
    def test_conflicting_original_origin_is_not_silent_join(self):
        r=per_sheet_id_metadata({
            "survival":{1:{"St":"K1","mark":1,"So":"S1"}},
            "performance":{1:{"St":"K1","mark":1,"So":"K2"}}},self.c)
        self.assertEqual(r["within_2014_crosssheet"]["matched_site_by_mark_with_disagreeing_origin_ID_sets"],1)
        self.assertFalse(r["within_2014_crosssheet"]["crosssheet_individual_2014_metadata_link_present"])
    def test_missing_mark_or_target_site_not_imputed(self):
        r=per_sheet_id_metadata({
            "survival":{1:{"St":"K1","mark":None,"So":"S1"},2:{"St":None,"mark":10,"So":"S1"}},
            "performance":{}},self.c)
        self.assertEqual(r["per_sheet"]["survival"]["rows_with_target_site_and_mark"],0)
    def test_source_only_numbers_are_identifiers_not_fitness(self):
        self.assertEqual(scalar_id(42.0),"42")
        self.assertEqual(scalar_id("K1"),"K1")
        self.assertIsNone(scalar_id(float("nan")))
        self.assertIsNone(scalar_id(0.33))
        self.assertEqual(rk_number((10<<2)|2),10)
        self.assertEqual(rk_number((100<<2)|3),1)
    def test_unexpected_metadata_field_access_prohibited(self):
        c=json.loads(CONTRACT.read_text(encoding="utf-8"))
        c["allow_row_cell_values_only_for_columns"].append("surv2009")
        with self.assertRaisesRegex(ValueError,"allowlist"):
            per_sheet_id_metadata({"survival":{},"performance":{}},c)
    def test_frozen_2007_absence_and_44of50(self):
        self.assertFalse(self.c["biological_limits"]["original_2007_raw_individual_tag_data_available"])
        self.assertTrue(self.c["biological_limits"]["previous_Bontrager_44of50_source_qualification_unchanged"])
        self.assertTrue(self.c["biological_limits"]["same_2004_adult_cohort_author_confirmed"])
if __name__=="__main__":
    unittest.main()
