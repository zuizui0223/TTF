"""Publisher-hosted Timema supplementary XLSX *file/sheet manifest only*.

Downloads authentic Nature Communications 2024 supplementary files, verifies ZIP
container structure, inspects workbook sheet names and raw file SHA256.
Does NOT read spreadsheet cell values, experiment outcomes, or original Dryad CSV.
"""
import hashlib
import io
import json
import urllib.error
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

BASE = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-46294-9/MediaObjects/"
FILES = {
    "source_data": "41467_2024_46294_MOESM6_ESM.xlsx",
    "supplementary_data_1_5": "41467_2024_46294_MOESM4_ESM.xlsx",
}
MAX_BYTES = 12_000_000


def file_manifest(name, filename):
    url = BASE + filename
    req = urllib.request.Request(url, headers={"User-Agent":"TTF-source-only-science/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=50) as response:
            raw = response.read(MAX_BYTES+1)
    except urllib.error.HTTPError as exc:
        return {"source":name,"filename":filename,
                "status":"PUBLISHER_ORIGINAL_NOT_ACCESSIBLE",
                "http_status":exc.code,"no_cells_read":True}
    if len(raw)>MAX_BYTES or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError(f"source {name}: invalid or overlarge original XLSX")
    with zipfile.ZipFile(io.BytesIO(raw)) as book:
        names=book.namelist()
        if "xl/workbook.xml" not in names:
            raise ValueError("workbook XML missing")
        root=ET.fromstring(book.read("xl/workbook.xml"))
        ns={"x":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        sheets=[{"name":a.get("name"),"sheetId":a.get("sheetId")} for a in root.findall(".//x:sheets/x:sheet",ns)]
        if not sheets: raise ValueError("no workbook sheets")
        worksheet_files=sorted(n for n in names if n.startswith("xl/worksheets/sheet") and n.endswith(".xml"))
    return {
        "status":"PUBLISHER_XLSX_FILE_AND_SHEET_METADATA_VERIFIED",
        "source":name,
        "url":url,"filename":filename,"byte_length":len(raw),
        "sha256":hashlib.sha256(raw).hexdigest(),
        "sheet_names":sheets,"worksheet_file_count":len(worksheet_files),
        "no_worksheet_cell_values_read":True,
        "row_count":"NOT_DETERMINED",
        "join_to_dryad":"NOT_VERIFIED",
    }


def main():
    rows=[file_manifest(key,value) for key,value in FILES.items()]
    print(json.dumps({"doi":"10.1038/s41467-024-46294-9",
                      "status":"PUBLISHER_SOURCE_MANIFEST_ONLY",
                      "results":rows,"no_experiment_outcome_opened":True},
                     indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
