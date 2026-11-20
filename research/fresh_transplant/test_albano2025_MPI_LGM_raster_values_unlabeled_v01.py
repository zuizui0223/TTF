"""Synthetic source-only regression tests for CRC-pinned MPI-ESM-P LGM samples."""
import copy,hashlib,json,struct,sys,unittest,zlib
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import albano2025_MPI_LGM_raster_values_unlabeled_v01 as m

def make_member(name,payload):
    comp=zlib.compressobj(6,zlib.DEFLATED,-15)
    encoded=comp.compress(payload)+comp.flush()
    crc=zlib.crc32(payload)&0xffffffff
    header=m.LOCAL.pack(b"PK\x03\x04",20,0,8,0,0,crc,len(encoded),len(payload),len(name),0)
    entry={"name":name,"local_header_offset":0,
           "compressed_bytes":len(encoded),"uncompressed_bytes":len(payload),
           "crc32":f"{crc:08x}"}
    return entry,header+name.encode()+encoded

class UnlabeledLGMTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.grid=json.loads(m.CELLS.read_text("utf-8"))
    def test_exact_95_cells_no_source_ids(self):
        g=self.grid
        self.assertEqual(len(g["cells"]),95)
        self.assertEqual(g["source_population_count"],96)
        self.assertEqual(set(g["cells"][0]),{"row","col"})
        self.assertNotIn("sites",g)
    def test_upper_left_grid_center(self):
        self.assertAlmostEqual(m.cell_center(0,0)[0],-180+m.RES/2)
        self.assertAlmostEqual(m.cell_center(0,0)[1],90-m.RES/2)
    def test_grid_out_of_bounds_rejected(self):
        for row,col in [(-1,10),(4320,1),(12,8640)]:
            with self.assertRaisesRegex(ValueError,"outside WorldClim"):m.cell_center(row,col)
    def test_pinned_MPIO_members_have_expected_zipped_offsets(self):
        self.assertEqual(m.ARCHIVE_BYTES,223064411)
        self.assertEqual(m.ENTRIES["BIO1"]["local_header_offset"],0)
        self.assertEqual(m.ENTRIES["BIO12"]["local_header_offset"],35777679)
        self.assertEqual(m.ENTRIES["BIO12"]["crc32"],"1f15f5db")
    def test_crc32_original_deflate_roundtrip(self):
        e,b=make_member("synthetic.tif",b"A Test GeoTIFF sample")
        self.assertEqual(m.remote_member_bytes(e,fetcher=lambda u,a,c:b),b"A Test GeoTIFF sample")
    def test_bad_compressed_member_crc_rejected(self):
        e,b=make_member("synthetic.tif",b"A Test GeoTIFF sample")
        e["crc32"]="00000000"
        with self.assertRaisesRegex(ValueError,"CRC32 drift"):
            m.remote_member_bytes(e,fetcher=lambda u,a,c:b)
    def test_local_header_identity_rejected(self):
        e,b=make_member("synthetic.tif",b"A Test GeoTIFF sample")
        data=b"BADD"+b[4:]
        with self.assertRaisesRegex(ValueError,"header/member geometry"):
            m.remote_member_bytes(e,fetcher=lambda u,a,c:data)
    def test_partial_http_payload_rejected(self):
        e,b=make_member("synthetic.tif",b"A Test GeoTIFF sample")
        with self.assertRaisesRegex(ValueError,"Partial/oversized"):
            m.remote_member_bytes(e,fetcher=lambda u,a,c:b[:-1])
    def test_different_grid_source_sha_rejected(self):
        d=copy.deepcopy(self.grid)
        d["source_csv_sha256"]="invalid"
        with self.assertRaisesRegex(ValueError,"provenance hash"):
            m.evaluate(d,fetcher=lambda u,a,c:b"")
    def test_two_synthetic_raster_members_never_make_effect(self):
        e1,b1=make_member("x.tif",b"A"*10)
        e2,b2=make_member("y.tif",b"B"*10)
        e2["local_header_offset"]=100
        with patch.dict(m.ENTRIES,{"BIO1":e1,"BIO12":e2},clear=True):
            def fake(url,start,end):
                return b1 if start==0 else b2
            def sample(payload,cells):
                value=10 if payload.startswith(b"A") else 20
                return [value]*len(cells),{"sampled_nominal_cells":len(cells),"valid_cells":len(cells)}
            x=m.evaluate(self.grid,fetcher=fake,sampler=sample)
        self.assertEqual(len(x["cell_profiles"]),95)
        self.assertTrue(all(k["BIO1_raw"]==10 and k["BIO12_raw"]==20 for k in x["cell_profiles"]))
        self.assertFalse(x["fitness_genetics_adaptation_responses_opened"])
        self.assertFalse(x["LGM_independent_predictor_gain_tested"])
if __name__=="__main__":unittest.main()
