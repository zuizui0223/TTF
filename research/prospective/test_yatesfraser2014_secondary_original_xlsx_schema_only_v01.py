"""Original XLSX ZIP-XML FIRST-TEXT-header metadata tests, no real outcomes."""
import io
import json
import sys
import unittest
import zipfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from yatesfraser2014_secondary_original_xlsx_schema_only_v01 import (
    CONTRACT,verified_record_download,schema_only
)

M="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
R="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
P="http://schemas.openxmlformats.org/package/2006/relationships"

def synthetic_xlsx(headers=("study","plant_id","survival"),extra_nonheader=True):
    wb=f'<workbook xmlns="{M}" xmlns:r="{R}"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
    rel=f'<Relationships xmlns="{P}"><Relationship Id="rId1" Type="{R}/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'
    strings='<sst xmlns="'+M+'" count="'+str(len(headers))+'" uniqueCount="'+str(len(headers))+'">'+''.join('<si><t>'+s+'</t></si>' for s in headers)+'</sst>'
    row='<row r="1">'+''.join(f'<c r="{chr(65+i)}1" t="s"><v>{i}</v></c>' for i in range(len(headers)))+'</row>'
    numeric='<row r="2"><c r="A2"><v>101</v></c><c r="C2"><v>0.98345</v></c></row>' if extra_nonheader else ""
    sheet=f'<worksheet xmlns="{M}"><sheetData>{row}{numeric}</sheetData></worksheet>'
    target=io.BytesIO()
    with zipfile.ZipFile(target,"w",compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/workbook.xml",wb)
        z.writestr("xl/_rels/workbook.xml.rels",rel)
        z.writestr("xl/sharedStrings.xml",strings)
        z.writestr("xl/worksheets/sheet1.xml",sheet)
    return target.getvalue()

class RaabovaSecondaryXLSXSchemaSafetyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    def verified_record(self,override_doi=None):
        original=self.c["original"]
        return {"id":original["official_zenodo_record_id"],
          "metadata":{"related_identifiers":[{"identifier":override_doi or original["original_data_doi"]}]},
          "files":[{"key":original["original_filename"],"size":original["original_bytes"],
                    "checksum":"md5:"+original["original_md5"],
                    "links":{"self":"https://zenodo.org/api/files/verifiedbucket/"+original["original_filename"]}}]}
    def test_exact_original_DOI_size_checksum_and_url_required(self):
        url=verified_record_download(self.verified_record(),self.c)
        self.assertTrue(url.startswith("https://zenodo.org/api/files/"))
        bad=self.verified_record("10.5061/dryad.other")
        with self.assertRaisesRegex(ValueError,"exact original"):
            verified_record_download(bad,self.c)
        bad=self.verified_record()
        bad["files"][0]["size"]+=1
        with self.assertRaisesRegex(ValueError,"size/MD5"):
            verified_record_download(bad,self.c)
    def test_numeric_response_rows_not_exported(self):
        s=schema_only(synthetic_xlsx(),self.c)
        self.assertTrue(s["any_explicit_plant_mark_candidate_header"])
        head=s["sheets"][0]
        self.assertEqual([x["early_text_label"] for x in head["potential_header_text_labels"]],
                         ["study","plant_id","survival"])
        self.assertEqual(head["first_potential_text_header_row"],1)
        self.assertEqual(head["explicit_field_plant_tag_heading_candidates"],["plant_id"])
        self.assertTrue(head["numeric_fitness_effect_values_opened"] is False)
        self.assertNotIn("0.98345",json.dumps(s))
    def test_missing_marker_in_first_text_header_not_dataset_absence(self):
        s=schema_only(synthetic_xlsx(("taxon","population","effect")),self.c)
        self.assertFalse(s["any_explicit_plant_mark_candidate_header"])
        self.assertEqual(s["sheets"][0]["explicit_field_plant_tag_heading_candidates"],[])
    def test_row_budget_and_no_outcome_authorization_fixed(self):
        lim=self.c["safe_XLSX_member_scope"]
        self.assertEqual(lim["first_text_rows_per_sheet"],3)
        self.assertTrue(lim["numeric_cells_or_formulas_NOT_SELECTED"])
        self.assertFalse(self.c["sources_allow_file_outcome_read"]) if "sources_allow_file_outcome_read" in self.c else None
        self.assertFalse(self.c["no_2007_to_2014_row_match_authorized"])
    def test_unexpected_workbook_relation_fails_closed(self):
        raw=synthetic_xlsx()
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            vals={n:z.read(n) for n in z.namelist()}
        vals["xl/_rels/workbook.xml.rels"]=vals["xl/_rels/workbook.xml.rels"].replace(
            b'worksheets/sheet1.xml',b'../sheet1.xml')
        target=io.BytesIO()
        with zipfile.ZipFile(target,"w") as z:
            for n,b in vals.items():z.writestr(n,b)
        with self.assertRaisesRegex(ValueError,"Unexpected XLSX workbook relation"):
            schema_only(target.getvalue(),self.c)
    def test_frozen_Dryad_meta_analysis_is_not_original_field_tag_evidence(self):
        self.assertIn("Meta-analysis",self.c["origin"]["published_data_scale"])
        self.assertFalse(self.c["infer_original_source_ID_only_if_explicit_tag_columns_and_primary_2007_study_provenance_cooccur"] is False)
        self.assertTrue(self.c["existing_44_of_50_stopped_Bontrager_gate_unchanged"])
if __name__=="__main__":
    unittest.main()
