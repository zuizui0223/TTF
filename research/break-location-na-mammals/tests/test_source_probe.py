import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from source_probe import Meta,probe_share
class Tests(unittest.TestCase):
 def test_metadata_parser(self):
  h=Meta();h.feed('<head><meta property="og:title" content="Genetic Breaks"/><link rel="canonical" href="https://figshare.com/articles/dataset/Test/123"></head>')
  self.assertEqual(h.title,["Genetic Breaks"]);self.assertEqual(len(h.canonical),1)
 def test_direct_hostname_reject(self):
  from source_probe import safe_request
  with self.assertRaises(ValueError):safe_request("https://example.com/test")
if __name__=="__main__":unittest.main()
