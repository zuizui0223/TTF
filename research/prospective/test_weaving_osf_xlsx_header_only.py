"""Synthetic XLSX ZIP tests: no source biological response rows accessed."""
import io
import sys
import unittest
import zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from weaving_osf_xlsx_header_only import file_schema


class SourceHeaderOnlyTests(unittest.TestCase):
    def test_first_row_headers_only(self):
        xmlworkbook=b'''<?xml version="1.0" encoding="utf-8"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
<sheets><sheet name="Data" sheetId="1" r:id="rId1"/></sheets></workbook>'''
        xmlrels=b'''<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Target="worksheets/sheet1.xml"
Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"/>
</Relationships>'''
        xmlsheet=b'''<?xml version="1.0" encoding="utf-8"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<dimension ref="A1:D900"/>
<sheetData>
<row r="1"><c r="A1" t="inlineStr"><is><t>Species</t></is></c>
<c r="B1" t="inlineStr"><is><t>Latitude</t></is></c>
<c r="C1" t="inlineStr"><is><t>Longitude</t></is></c>
<c r="D1" t="inlineStr"><is><t>ARR</t></is></c></row>
<row r="2"><c r="A2" t="inlineStr"><is><t>Testus syntheticus</t></is></c>
<c r="B2"><v>11</v></c><c r="C2"><v>33</v></c>
<c r="D2"><v>987654321.1234</v></c></row></sheetData></worksheet>'''
        b=io.BytesIO()
        with zipfile.ZipFile(b,'w') as z:
            z.writestr("xl/workbook.xml",xmlworkbook)
            z.writestr("xl/_rels/workbook.xml.rels",xmlrels)
            z.writestr("xl/worksheets/sheet1.xml",xmlsheet)
        got=file_schema(b.getvalue())
        self.assertEqual(len(got),1)
        self.assertTrue(got[0]['contains_latitude_header'])
        self.assertTrue(got[0]['contains_longitude_header'])
        self.assertTrue(got[0]['contains_taxon_header'])
        self.assertEqual(got[0]['sheet_dimension_from_source'],'A1:D900')
        self.assertNotIn('987654321.1234',str(got))
        self.assertFalse(got[0]['row_level_biological_response_values_read'])

    def test_invalid_zip_rejected(self):
        with self.assertRaises(ValueError):
            file_schema(b'not an xlsx')


if __name__=='__main__':
    unittest.main()
