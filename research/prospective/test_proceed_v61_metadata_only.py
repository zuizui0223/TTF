"""Synthetic-only tests for PROCEED v6.1 Borealis catalog interpretation."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed_v61_borealis_source_metadata import source_catalog,api_metadata_get

class Tests(unittest.TestCase):
    def test_no_datatable_is_required_for_public_metadata(self):
        obj={"status":"OK","data":{"id":100,"persistentId":"doi:10.5683/SP3/NXSL3Q",
            "latestVersion":{"versionNumber":6,"versionMinorNumber":1,
                "files":[{"label":"proceed.csv","description":"Data",
                          "dataFile":{"id":50,"filename":"proceed.csv","contentType":"text/tab-separated-values",
                                      "filesize":400,"md5":"abcd"}}],
                "metadataBlocks":{"citation":{"fields":[{"typeName":"title","value":"PROCEED"}]}}}}}
        r=source_catalog(obj)
        self.assertEqual(r["file_count"],1)
        self.assertEqual(r["persistent_id"],"doi:10.5683/SP3/NXSL3Q")
        self.assertFalse(r["phenotypic_effect_size_rows_opened"])
        self.assertFalse(r["history_adaptation_statistical_effect_fitted"])
        self.assertEqual(r["files"][0]["file_id"],50)
        self.assertNotIn("rowvalues",r["files"][0])

    def test_malformed_fails(self):
        with self.assertRaises(ValueError):source_catalog({"status":"ERROR"})

    def test_not_a_dataset_file_download(self):
        with self.assertRaises(ValueError):
            api_metadata_get("https://example.org/api/datasets/1",10)
        with self.assertRaises(ValueError):
            api_metadata_get("https://borealisdata.ca/api/access/datafile/123",10)

if __name__=="__main__":unittest.main()
