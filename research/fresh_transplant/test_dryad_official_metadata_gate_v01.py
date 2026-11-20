"""12 offline original Dryad API metadata-only contract tests."""
import copy,json,sys,unittest
from pathlib import Path
from urllib.error import HTTPError
sys.path.insert(0,str(Path(__file__).resolve().parent))
import dryad_official_metadata_gate_v01 as a

class DryadManifestTests(unittest.TestCase):
    def setUp(self):self.d=json.loads(a.CONTRACT.read_text("utf-8"))
    @staticmethod
    def fake_get(url):
        if "/versions?" in url:
            return {"_embedded":{"stash:versions":[
                {"id":10,"versionNumber":1},{"id":20,"versionNumber":2}]}}
        if "/versions/20/files?" in url:
            return {"_embedded":{"stash:files":[
                {"id":4337649,"path":"Albano_et_al_2025_EM_Experiment_Data.csv","size":421800,"digest":"abc","digestType":"md5"},
                {"id":4692223,"path":"latlongs.csv","size":375},
                {"id":4692224,"path":"tplantlatlongs.csv","size":157},
                {"id":911,"path":"Clarkia-LTREB-transplant-archive-R2.zip","size":7500000}]}}
        raise ValueError("unrecognized synthetic API URL")
    def test_three_original_dois_and_version_selector(self):
        r=a.evaluate(self.d,fetcher=self.fake_get)
        self.assertEqual(len(r["studies"]),3)
        self.assertTrue(all(s["latest_public_version_id"]==20 for s in r["studies"]))
    def test_all_original_files_matched_synthetic(self):
        r=a.evaluate(self.d,fetcher=self.fake_get)
        self.assertTrue(all(s["status"]=="ALL_PINNED_FILES_REGISTERED_IN_OFFICIAL_METADATA" for s in r["studies"]))
        self.assertEqual([s["original_pinned_files"] for s in r["studies"]],[1,2,1])
    def test_version_drift_detected_not_autopromoted(self):
        def swapped(url):
            r=self.fake_get(url)
            if "/versions/20/files?" in url:
                r["_embedded"]["stash:files"][0]["id"]=444
            return r
        s=a.check_study(self.d["candidates"][0],fetcher=swapped)
        self.assertEqual(s["status"],"ORIGINAL_FILE_METADATA_PARTIALLY_UNQUALIFIED")
        self.assertFalse(s["files_metadata"][0]["file_id_matches_frozen_page"])
    def test_401_transport_fails_closed(self):
        def blocked(url):
            raise HTTPError(url,401,"unauthorized",{},None)
        s=a.check_study(self.d["candidates"][0],fetcher=blocked)
        self.assertEqual(s["status"],"PUBLIC_VERSION_MANIFEST_NOT_CONFIRMED")
        self.assertEqual(s["version_request_http_error_code"],401)
        self.assertFalse(s["original_source_GPS_verified"])
    def test_missing_file_no_false_match(self):
        def absent(url):
            r=self.fake_get(url)
            if "/versions/20/files?" in url:
                r["_embedded"]["stash:files"]=[]
            return r
        s=a.check_study(self.d["candidates"][1],fetcher=absent)
        self.assertEqual(s["status"],"ORIGINAL_FILE_METADATA_PARTIALLY_UNQUALIFIED")
    def test_original_api_is_only_https_dryad(self):
        self.assertTrue(a.API.startswith("https://datadryad.org/api/v2/"))
    def test_flatten_Hal_versions(self):
        got=a.flatten_api_records({"_embedded":{"stash:versions":[{"id":7}]}},"version")
        self.assertEqual(got,[{"id":7}])
    def test_flatten_unknown_schema_rejected(self):
        with self.assertRaisesRegex(ValueError,"Unknown"):
            a.flatten_api_records({"_embedded":{}},"version")
    def test_frozen_selection_fails_closed(self):
        d=copy.deepcopy(self.d)
        d["candidates"].reverse()
        with self.assertRaisesRegex(ValueError,"Frozen candidate"):
            a.evaluate(d,fetcher=self.fake_get)
    def test_no_raw_original_data_or_LGM_proven(self):
        r=a.evaluate(self.d,fetcher=self.fake_get)
        self.assertFalse(r["fitness_or_genetics_response_opened"])
        self.assertFalse(r["LGM21ka_contrast_qualified"])
        self.assertEqual(r["source_GPS_independently_certified"],0)
        self.assertTrue(all(not s["file_or_response_bytes_downloaded"] for s in r["studies"]))
    def test_dryad_HAL_version_record_without_inline_integer_id(self):
        r={"_links":{"self":{"href":"/api/v2/versions/261740"},
                     "stash:files":{"href":"/api/v2/versions/261740/files"}},
           "versionNumber":2}
        self.assertEqual(a.linked_record_id(r,"versions"),261740)
    def test_dryad_HAL_file_record_without_inline_integer_id(self):
        r={"path":"latlongs.csv","_links":{"stash:download":{
             "href":"/api/v2/files/4692223/download"}}}
        self.assertEqual(a.linked_record_id(r,"files"),4692223)
    def test_complete_version_and_file_Hal_only(self):
        def got(url):
            if "/versions?" in url:
                return {"_embedded":{"stash:versions":[{"versionNumber":2,
                   "_links":{"self":{"href":"/api/v2/versions/361800"},
                             "stash:files":{"href":"/api/v2/versions/361800/files"}}}]}}
            if "/versions/361800/files?" in url:
                return {"_embedded":{"stash:files":[
                   {"path":"latlongs.csv","size":375,
                    "_links":{"stash:download":{"href":"/api/v2/files/4692223/download"}}},
                   {"path":"tplantlatlongs.csv","size":157,
                    "_links":{"stash:download":{"href":"/api/v2/files/4692224/download"}}}
                ]}}
            raise ValueError("Wrong fake URL")
        s=a.check_study(self.d["candidates"][1],fetcher=got)
        self.assertEqual(s["latest_public_version_id"],361800)
        self.assertEqual(s["status"],"ALL_PINNED_FILES_REGISTERED_IN_OFFICIAL_METADATA")
        self.assertTrue(all(x["file_id_matches_frozen_page"] for x in s["files_metadata"]))
    def test_lexical_DOI_held_to_Dryad(self):
        d=copy.deepcopy(self.d["candidates"][0])
        d["doi"]="https://evil.example/data"
        with self.assertRaisesRegex(ValueError,"Dryad scope"):
            a.check_study(d,fetcher=self.fake_get)
    def test_original_dataset_manifest_unchanged(self):
        r=a.evaluate(self.d,fetcher=self.fake_get)
        self.assertTrue(r["old_Bontrager_44of50_stop_unchanged"])
        self.assertTrue(r["frozen_selection_before_biological_response"])

if __name__=="__main__":unittest.main()
