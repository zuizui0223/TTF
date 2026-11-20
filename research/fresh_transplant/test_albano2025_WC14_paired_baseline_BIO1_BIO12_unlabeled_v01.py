"""Synthetic CRC/identity tests for WorldClim 1.4 1960–1990 matched BIL source cells."""
import copy,io,json,sys,unittest,zlib
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import albano2025_WC14_paired_baseline_BIO1_BIO12_unlabeled_v01 as a

class SourceBaselineTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.g=json.loads(a.CELLS.read_text("utf-8"))
  cls.m=json.loads(a.MANIFEST.read_text("utf-8"))
 def test_official_four_members(self):
  self.assertEqual({v["filename"] for v in self.m["members"]},{"bio1.bil","bio1.hdr","bio12.bil","bio12.hdr"})
  self.assertEqual(self.m["original_archive_total_bytes"],129319755)
 def test_original_95_grid_cell_fingerprint(self):
  self.assertEqual(self.g["unique_grid_cell_count"],95)
  self.assertEqual(self.g["source_population_count"],96)
  self.assertNotIn("sites",self.g)
 def test_original_BIL_expanded_byte_count(self):
  b=[x for x in self.m["members"] if x["filename"].endswith(".bil")]
  self.assertEqual([x["uncompressed_bytes"] for x in b],[62208000]*2)
 def test_globe_cell_centers(self):
  self.assertAlmostEqual(a.cell_center(0,0)[0],-180+a.RES/2)
  self.assertAlmostEqual(a.cell_center(0,0)[1],90-a.RES/2)
 def test_unrecognized_source_zip_rejected(self):
  m=copy.deepcopy(self.m);m["original_archive_total_bytes"]=1
  with self.assertRaisesRegex(ValueError,"archive identity drift"):
   a.evaluate(self.g,m,fetcher=lambda *_:b"")
 def test_duplicate_grid_cell_rejected(self):
  g=copy.deepcopy(self.g);g["cells"][-1]=g["cells"][0]
  with self.assertRaisesRegex(ValueError,"95-cell"):
   a.evaluate(g,self.m,fetcher=lambda *_:b"")
 def test_four_synthetic_members_validate_CRC_and_no_effect(self):
  compressed={}
  sample_names=("bio1.bil","bio1.hdr","bio12.bil","bio12.hdr")
  entries=[]
  for offset,name in enumerate(sample_names):
   raw=(name+"-test").encode()*3
   z=zlib.compressobj(6,zlib.DEFLATED,-15);encoded=z.compress(raw)+z.flush()
   crc=zlib.crc32(raw)&0xffffffff
   entry={"filename":name,"local_header_offset":offset*100,
          "compressed_bytes":len(encoded),"uncompressed_bytes":len(raw),
          "CRC32_hex":f"{crc:08x}","method":8}
   sig=a.HEADER.pack(b"PK\x03\x04",20,0,8,0,0,crc,len(encoded),len(raw),len(name),0)
   compressed[offset*100]=sig+name.encode()+encoded;entries.append(entry)
  m=copy.deepcopy(self.m);m["members"]=entries
  def fetch(url,start,end):return compressed[start]
  def sample(bil,hdr,points):
   v=1.0 if bil.startswith(b"bio1") else 12.0
   return [v]*len(points),{"valid_values":len(points)}
  r=a.evaluate(self.g,m,fetcher=fetch,sampler=sample)
  self.assertEqual(len(r["cells"]),95)
  self.assertEqual(r["cells"][0]["BIO1_raw"],1)
  self.assertEqual(r["cells"][0]["BIO12_raw"],12)
  self.assertFalse(r["phenotype_fitness_genetic_responses_used"])
  self.assertFalse(r["LGM_22ka_effect_tested"])
 def test_corrupted_member_crc_rejected(self):
  raw=b"TEST";z=zlib.compressobj(6,zlib.DEFLATED,-15);enc=z.compress(raw)+z.flush()
  e={"filename":"bio1.bil","local_header_offset":0,
     "compressed_bytes":len(enc),"uncompressed_bytes":len(raw),
     "CRC32_hex":"00000000"}
  b=a.HEADER.pack(b"PK\x03\x04",20,0,8,0,0,zlib.crc32(raw),len(enc),len(raw),8,0)+b"bio1.bil"+enc
  with self.assertRaisesRegex(ValueError,"CRC32"):
   a.read_member(e,fetcher=lambda url,start,end:b)
if __name__=="__main__":unittest.main()
