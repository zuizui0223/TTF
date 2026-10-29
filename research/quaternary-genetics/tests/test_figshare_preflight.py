import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from figshare_preflight import normalize_metadata

class TestResponseBlindPreflight(unittest.TestCase):
    def test_valid_metadata_without_he(self):
        m=normalize_metadata({"id":28407989,"doi":"10.6084/m9.figshare.28407989",
            "files":[{"id":2,"name":"Population.csv","size":12345,
             "download_url":"https://figshare.com/example","He":0.91}]})
        self.assertEqual(m["file_count"],1)
        self.assertFalse(m["biology_response_values_accessed"])
        self.assertNotIn("0.91",repr(m))
    def test_invalid_article_stops(self):
        with self.assertRaises(ValueError):
            normalize_metadata({"id":1,"files":[]})
if __name__=="__main__":unittest.main()
