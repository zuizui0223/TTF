"""Synthetic BIFF8 only; no published outcome or individual data rows needed."""
import json
import struct
import sys
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from pankova2014_s1_xls_schema_only_v01 import (
    CONTRACT,extract_early_text_headers,records,sst_first_fragment_strings
)

def rec(sid,payload):
    return struct.pack("<HH",sid,len(payload))+payload

def fake_original_workbook(columns=("plant_id","site","survival"),include_numeric=False,late=False):
    sst=struct.pack("<II",len(columns),len(columns))
    for name in columns:
        s=name.encode("cp1252")
        sst+=struct.pack("<HB",len(s),0)+s
    name=b"Data"
    gs=rec(0x0809,b"\0"*8)
    book=rec(0x0085,struct.pack("<IBBB",0,0,0,len(name))+b"\0"+name)
    strrec=rec(0x00FC,sst)
    endglobals=rec(0x000A,b"")
    offset=len(gs)+len(book)+len(strrec)+len(endglobals)
    book=rec(0x0085,struct.pack("<IBBB",offset,0,0,len(name))+b"\0"+name)
    ws=rec(0x0809,b"\0"*8)
    for col in range(len(columns)):
        ws+=rec(0x00FD,struct.pack("<HHHI",0,col,0,col))
    if include_numeric:
        ws+=rec(0x0203,struct.pack("<HHHd",3,5,0,999999.))
    if late:
        ws+=rec(0x00FD,struct.pack("<HHHI",15,0,0,0))
    ws+=rec(0x000A,b"")
    return gs+book+strrec+endglobals+ws

class SourceOnlyLegacyXLSSchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    def test_original_synthetic_OLE_workbook_first_row_names(self):
        a=extract_early_text_headers(fake_original_workbook(),self.c)
        self.assertEqual(len(a["sheets"]),1)
        self.assertEqual(a["sheets"][0]["chosen_early_header_row_1_based"],1)
        self.assertEqual([v["column_name"] for v in a["sheets"][0]["column_names_from_early_header"]],
                         ["plant_id","site","survival"])
        self.assertTrue(a["any_explicit_tag_id_column_label_found"])
        self.assertFalse(a["sheets"][0]["any_numeric_or_fitness_values_selected"])
    def test_real_values_in_numeric_cells_never_interpreted(self):
        a=extract_early_text_headers(fake_original_workbook(include_numeric=True,late=True),self.c)
        self.assertEqual(a["sheets"][0]["chosen_early_header_row_1_based"],1)
        self.assertEqual(len(a["sheets"][0]["column_names_from_early_header"]),3)
    def test_no_id_header_does_not_mean_no_individual_ids_in_workbook(self):
        a=extract_early_text_headers(fake_original_workbook(("site","flower_count","biomass")),self.c)
        self.assertFalse(a["any_explicit_tag_id_column_label_found"])
    def test_cross_CONTINUE_sst_fails_instead_of_inventing_header(self):
        with self.assertRaisesRegex(ValueError,"CONTINUE"):
            sst_first_fragment_strings(struct.pack("<II",2,2)+struct.pack("<HB",20,0)+b"abc",1)
    def test_corrupted_sheet_offset_fails_closed(self):
        raw=bytearray(fake_original_workbook())
        # Offset begins at record payload for first BOUNDSHEET at 16.
        struct.pack_into("<I",raw,16,9999999)
        with self.assertRaisesRegex(ValueError,"offset"):
            extract_early_text_headers(bytes(raw),self.c)
    def test_frozen_rows_only_12_before_original_access(self):
        self.assertEqual(self.c["method"]["max_header_rows_per_sheet"],12)
        self.assertTrue(self.c["method"]["forbid_numeric_cell_value_field_read"])
        self.assertFalse(self.c["comparison_boundary"]["individual_ID_from_2014_alone_does_not_confirm_crosspaper_ID_match"] is False)
        self.assertTrue(self.c["comparison_boundary"]["existing_44of50_bontrager_gate_stays_terminal"])
if __name__=="__main__":unittest.main()
