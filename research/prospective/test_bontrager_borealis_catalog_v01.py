"""Test official alternate DOI identity and strict metadata-only output."""
import unittest
from bontrager_borealis_catalog_v01 import parse_catalog

class CatalogTests(unittest.TestCase):
    def test_catalog_only(self):
        src={"status":"OK","data":{"persistentId":"doi:10.5683/SP2/DU3ZTQ","id":12,
            "latestVersion":{"versionNumber":3,"versionMinorNumber":0,
                "files":[{"label":"data.zip","restricted":False,
                          "dataFile":{"id":4,"filename":"data.zip",
                                      "filesize":900,"md5":"a"*32,
                                      "contentType":"application/zip"}}]}}}
        x=parse_catalog(src)
        self.assertEqual(x["file_count"],1)
        self.assertEqual(x["files"][0]["md5"],"a"*32)
        self.assertFalse(x["fitness_records_read"])
        self.assertFalse(x["independent_from_Dryad_or_Hargreaves_confirmed"])
    def test_mismatched_doi_stops(self):
        with self.assertRaises(ValueError):
            parse_catalog({"status":"OK","data":{"persistentId":"doi:OTHER",
                "latestVersion":{"files":[]}}})
if __name__=="__main__":unittest.main()
