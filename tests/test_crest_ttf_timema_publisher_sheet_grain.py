"""Offline synthetic workbook schema checks, no real Timema outcomes."""
import importlib.util
from pathlib import Path
import hashlib
import io
import zipfile
import pytest

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_timema_publisher_sheet_grain_20261115.py"
s=importlib.util.spec_from_file_location("timema_sheet_grain", FILE)
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

def minimal():
    wb=b'''<?xml version="1.0"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Test" sheetId="1" r:id="rId1"/></sheets></workbook>'''
    rel=b'''<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'''
    shared=b'''<?xml version="1.0"?><sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><si><t>Pair</t></si><si><t>Species1</t></si></sst>'''
    sheet=b'''<?xml version="1.0"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><dimension ref="A1:C2"/><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row><row r="2"><c r="C2"><v>0.9876543</v></c></row></sheetData></worksheet>'''
    buf=io.BytesIO()
    with zipfile.ZipFile(buf,"w") as z:
        z.writestr("xl/workbook.xml",wb)
        z.writestr("xl/_rels/workbook.xml.rels",rel)
        z.writestr("xl/sharedStrings.xml",shared)
        z.writestr("xl/worksheets/sheet1.xml",sheet)
    return buf.getvalue()

def test_source_only_header_hints_and_shape():
    raw=minimal()
    result=m.parse_workbook(raw,hashlib.sha256(raw).hexdigest())
    assert result["worksheet_count"]==1
    t=result["sheets"][0]
    assert t["declared_range"]=="A1:C2"
    assert t["row_element_count"]==2
    assert {x["header_hint"] for x in t["first_four_rows_header_hints"]}=={"Pair","Species1"}
    assert "0.9876543" not in repr(result)
    assert t["original_numeric_cells_not_decoded"] is True

def test_hash_drift_stops():
    with pytest.raises(ValueError,match="SHA256"):
        m.parse_workbook(minimal(),"0"*64)
