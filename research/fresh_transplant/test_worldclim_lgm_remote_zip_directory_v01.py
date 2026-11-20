"""Ten offline tests for bounded WorldClim ZIP central-directory inspection."""
import copy,io,json,struct,sys,unittest,zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import worldclim_lgm_remote_zip_directory_v01 as a

def fixture():
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w",zipfile.ZIP_DEFLATED) as z:
        z.writestr("bio1.tif",b"0123"*20)
        z.writestr("bio12.tif",b"4567"*30)
        z.writestr("bio12.tfw",b"0.04166667")
    return b.getvalue()

class ZipRangeTests(unittest.TestCase):
    def setUp(self):
        self.data=fixture()
        self.contract=json.loads(a.CONTRACT.read_text("utf-8"))
    def parse(self,raw=None):
        data=raw if raw is not None else self.data
        return a.parse_zip_central(data,len(data),lambda start,end:data[start:end+1])
    def test_only_official_frozen_urls(self):
        self.assertEqual(len(a.source_urls(self.contract)),4)
    def test_three_members_parsed(self):
        t=self.parse()
        self.assertEqual([x["filename"] for x in t],["bio1.tif","bio12.tif","bio12.tfw"])
    def test_exact_size_and_deflated_method(self):
        t=self.parse()
        self.assertEqual([x["uncompressed_bytes"] for x in t],[80,120,10])
        self.assertEqual([x["method"] for x in t],[8,8,8])
    def test_archive_local_header_offsets_have_sources(self):
        t=self.parse()
        self.assertTrue(all(self.data[x["local_header_offset"]:].startswith(b"PK\x03\x04") for x in t))
    def test_corrupt_eocd_rejected(self):
        data=self.data[:-22]
        with self.assertRaisesRegex(ValueError,"EOCD"):self.parse(data)
    def test_excess_comment_rejected(self):
        data=self.data+b"extra"
        with self.assertRaisesRegex(ValueError,"central directory"):self.parse(data)
    def test_corrupt_central_signature_rejected(self):
        data=bytearray(self.data)
        ix=data.find(b"PK\x01\x02");data[ix:ix+4]=b"BAD!"
        with self.assertRaisesRegex(ValueError,"signature"):self.parse(bytes(data))
    def test_changed_public_GCM_url_rejected(self):
        changed=copy.deepcopy(self.contract)
        changed["models"][0]["official_url"]="http://example.com/file.zip"
        with self.assertRaisesRegex(ValueError,"GCM/source list"):a.source_urls(changed)
    def test_file_metadata_not_climate_values(self):
        url=next(iter(a.source_urls(self.contract).values()))
        def fake_ranger(url,start,end):
            if start is None:return self.data,len(self.data)
            return self.data[start:end+1],len(self.data)
        receipt=a.inspect_one(url,ranger=fake_ranger)
        self.assertFalse(receipt["raster_data_values_read"])
        self.assertEqual(receipt["archive_file_count"],3)
    def test_access_failure_gated(self):
        x=a.evaluate(self.contract,inspector=lambda url: (_ for _ in ()).throw(ValueError("No HTTP range")))
        self.assertEqual(len(x["archives"]),4)
        self.assertTrue(all(s["status"]=="OFFICIAL_ZIP_INDEX_ACCESS_NOT_CONFIRMED" for s in x["archives"].values()))
        self.assertFalse(x["all_LGM_climate_values_extracted"])
if __name__=="__main__":unittest.main()
