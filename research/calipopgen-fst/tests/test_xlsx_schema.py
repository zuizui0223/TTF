import io,sys,unittest,zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from inspect_xlsx_schema import schema_only
N='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
class TestXlsxHeader(unittest.TestCase):
 def test_only_header_read(self):
  b=io.BytesIO()
  with zipfile.ZipFile(b,'w') as z:
   z.writestr('xl/sharedStrings.xml','<sst xmlns="'+N+'"><si><t>GenDist</t></si><si><t>HabitatType</t></si></sst>')
   z.writestr('xl/worksheets/sheet1.xml','<worksheet xmlns="'+N+'"><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row><row r="2"><c r="A2"><v>0.92829372819</v></c></row></sheetData></worksheet>')
  r=schema_only(b.getvalue())
  self.assertEqual(r["columns"],["GenDist","HabitatType"])
  self.assertEqual(r["rows_in_first_worksheet"],1)
  self.assertNotIn("0.92829372819",str(r))
if __name__=="__main__":unittest.main()
