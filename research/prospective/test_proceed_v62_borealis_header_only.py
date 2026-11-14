"""Synthetic-only source qualification tests; do not open PROCEED effect-size rows."""
from __future__ import annotations
import io
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))

from proceed_v62_borealis_header_only import stream_one_header, verify_frozen_catalog, CONTRACT


class FakeResponse:
    def __init__(self, header):
        self.data = io.BytesIO(header + b'\nS1\tS2\t99999999\n')
        self.headers = {"Content-Type": "text/tab-separated-values"}
        self.readline_calls = []
    def __enter__(self): return self
    def __exit__(self, *args): return None
    def readline(self, n):
        self.readline_calls.append(n)
        return self.data.readline(n)

class SchemaOnlyTests(unittest.TestCase):
    def test_only_first_tsv_header_consumed(self):
        response=FakeResponse(b'StudyID\tOriginLat\tOriginLon\tGeneticEffect')
        with patch("proceed_v62_borealis_header_only.urlopen",return_value=response) as opener:
            result=stream_one_header(1024206)
        self.assertEqual(len(response.readline_calls),1)
        self.assertEqual(result["header_columns"],["StudyID","OriginLat","OriginLon","GeneticEffect"])
        self.assertFalse(result["rows_two_plus_read"])
        self.assertNotIn("99999999",repr(result))
        self.assertIn("/api/access/datafile/1024206",opener.call_args[0][0].full_url)

    def test_source_file_manifest_is_immutable(self):
        frozen=json.loads(CONTRACT.read_text())["current_source"]
        files=[{"dataFile":{"id":f["id"],"filename":f["name"],
            "filesize":f["bytes"],"md5":f["md5"],"contentType":"text/tab-separated-values"},
            "restricted":False} for f in frozen["files"]]
        m={"status":"OK","data":{"id":frozen["borealis_dataset_id"],
            "latestVersion":{"versionNumber":4,"versionMinorNumber":0,
                             "files":files}}}
        self.assertEqual(len(verify_frozen_catalog(m,frozen)),5)
        m["data"]["latestVersion"]["files"][0]["dataFile"]["md5"]="WRONG"
        with self.assertRaisesRegex(ValueError,"Source file"):
            verify_frozen_catalog(m,frozen)

    def test_reject_wrong_dataset_and_changed_version(self):
        frozen=json.loads(CONTRACT.read_text())["current_source"]
        for wrong in [{"status":"ERROR"},{"status":"OK","data":{"id":5,"latestVersion":{}}},
                      {"status":"OK","data":{"id":frozen["borealis_dataset_id"],
                        "latestVersion":{"versionNumber":5,"versionMinorNumber":0}}}]:
            with self.assertRaises(ValueError):verify_frozen_catalog(wrong,frozen)


if __name__=="__main__": unittest.main()
