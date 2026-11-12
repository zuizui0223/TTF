"""Synthetic-only OSF file listing tests; never download biological data."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from weaving_osf_metadata_only_probe import crawl,safe_osf_url,extract_entry

class OSFMetadataTests(unittest.TestCase):
    def test_only_metadata_api_hosts(self):
        safe_osf_url("https://api.osf.io/v2/nodes/cbhv4/files/osfstorage/")
        for url in ["https://osf.io/cbhv4/download","http://api.osf.io/v2/nodes/1/",
                    "https://attacker.example.com/v2/nodes/1/"]:
            with self.assertRaises(ValueError):
                safe_osf_url(url)

    def test_file_listing_does_not_follow_download_link(self):
        root="https://api.osf.io/v2/nodes/cbhv4/files/osfstorage/"
        data={
          "data":[{"id":"file123","attributes":{"kind":"file","name":"StudyData.csv",
                     "size":345,"extra":{"hashes":{"md5":"abc"}}},
                   "links":{"download":"https://osf.io/malicious-test/download"}},{"id":"folder123",
                    "attributes":{"kind":"folder","name":"Group"},
                    "relationships":{"files":{"links":{"related":{"href":"https://api.osf.io/v2/files/folder123/children/"}}}}}],
          "links":{"next":None}
        }
        nested={"data":[{"id":"file456","attributes":{"kind":"file","name":"README.txt",
                            "size":12},"links":{"download":"https://example.com/raw"}}],"links":{}}
        calls=[]
        def fake_fetch(url):
            calls.append(url)
            return data if url==root else nested
        result=crawl(root,fetch=fake_fetch)
        self.assertEqual(result["metadata_pages_retrieved"],2)
        self.assertEqual(result["files"][0]["id"],"file456")
        self.assertEqual(result["files"][1]["id"],"file123")
        self.assertEqual(len(calls),2)
        self.assertTrue(all("download" not in z for z in result["files"]))

if __name__=="__main__":unittest.main()
