"""Offline known-truth regression tests for frozen original Dryad file receipt."""
import copy,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_frozen_official_dryad_identity_v01 import RECEIPT,verify

def synthetic_observed(r):
    return {
      "schema":"ttf_fresh_official_Dryad_DOI_file_manifest_v01",
      "studies":[{
        "study_id":d["study_id"],"dataset_doi":d["doi"],
        "latest_public_version_id":d["version_id"],
        "status":"ALL_PINNED_FILES_REGISTERED_IN_OFFICIAL_METADATA",
        "files_metadata":[{
          "expected_original_path":f["path"],
          "official_file_record_id":f["id"],
          "official_file_size_bytes":f["size_bytes"],
          "official_file_digest":f["sha256_from_repository_metadata"],
          "official_file_digest_type":"sha-256",
          "status":"ORIGINAL_FILE_METADATA_MATCHES_FROZEN_SOURCE",
          "file_id_matches_frozen_page":True,
        } for f in d["files"]]
      } for d in r["datasets"]],
      "fitness_or_genetics_response_opened":False,
      "source_GPS_independently_certified":0,
      "LGM21ka_contrast_qualified":False,
    }
class FrozenDryadIdentityTests(unittest.TestCase):
    def setUp(self):
        self.receipt=json.loads(RECEIPT.read_text("utf-8"))
        self.live=synthetic_observed(self.receipt)
    def test_exact_three_datasets_four_original_files(self):
        r=verify(self.receipt,self.live)
        self.assertEqual((r["original_dataset_count"],r["original_file_count"]),(3,4))
    def test_SHA256_version_and_size_are_official_metadata_only(self):
        r=verify(self.receipt,self.live)
        self.assertIn("NOT_FILE_BYTE_HASH",r["frozen_receipt_provenance"])
    def test_no_new_GPS_or_fitness_claim(self):
        r=verify(self.receipt,self.live)
        self.assertEqual(r["new_independent_GPS_populations_certified"],0)
        self.assertFalse(r["new_21ka_LGM_fitness_effect_tested"])
    def test_official_version_drift_fails(self):
        d=copy.deepcopy(self.live);d["studies"][1]["latest_public_version_id"]+=1
        with self.assertRaisesRegex(ValueError,"version identity drift"):verify(self.receipt,d)
    def test_filename_drift_fails(self):
        d=copy.deepcopy(self.live);d["studies"][0]["files_metadata"][0]["expected_original_path"]="nonoriginal.csv"
        with self.assertRaisesRegex(ValueError,"file identity/digest drift"):verify(self.receipt,d)
    def test_hash_drift_fails(self):
        d=copy.deepcopy(self.live);d["studies"][2]["files_metadata"][0]["official_file_digest"]="0"*64
        with self.assertRaisesRegex(ValueError,"file identity/digest drift"):verify(self.receipt,d)
    def test_download_metadata_not_accessible(self):
        d=copy.deepcopy(self.live);d["studies"][0]["status"]="PUBLIC_VERSION_MANIFEST_NOT_CONFIRMED"
        with self.assertRaisesRegex(ValueError,"not fully accessible"):verify(self.receipt,d)
    def test_record_addition_fails(self):
        d=copy.deepcopy(self.live);d["studies"].append(copy.deepcopy(d["studies"][0]))
        with self.assertRaisesRegex(ValueError,"Study family changed"):verify(self.receipt,d)
    def test_file_id_drift_fails(self):
        d=copy.deepcopy(self.live);d["studies"][1]["files_metadata"][0]["official_file_record_id"]=2
        with self.assertRaisesRegex(ValueError,"file identity/digest drift"):verify(self.receipt,d)
    def test_biological_response_promotion_fails(self):
        d=copy.deepcopy(self.live);d["fitness_or_genetics_response_opened"]=True
        with self.assertRaisesRegex(ValueError,"unexpectedly promoted"):verify(self.receipt,d)
    def test_published_sha256_is_64_hex(self):
        for d in self.receipt["datasets"]:
            for f in d["files"]:
                self.assertEqual(len(bytes.fromhex(f["sha256_from_repository_metadata"])),32)
    def test_wrong_digest_type_rejected(self):
        d=copy.deepcopy(self.live);d["studies"][2]["files_metadata"][0]["official_file_digest_type"]="md5"
        with self.assertRaisesRegex(ValueError,"digest format"):verify(self.receipt,d)
if __name__=="__main__":unittest.main()
