"""Synthetic tests for SHA-pinned mirror ZIP central directory only."""
import hashlib,io,unittest,zipfile
from bontrager_mirror_zip_v01 import verify_original_zip

def fake_zip():
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:
        z.writestr("original/fitness.csv","fitness\n9999999\n")
    return b.getvalue()

class ArchiveIdentity(unittest.TestCase):
    def test_sha_and_md5_match_only_file_names_visible(self):
        raw=fake_zip()
        lock={"expected_file_bytes":len(raw),"dryad_sha256":hashlib.sha256(raw).hexdigest(),
              "mirror_md5":hashlib.md5(raw).hexdigest()}
        result=verify_original_zip(raw,lock)
        self.assertEqual(result["member_names"],["original/fitness.csv"])
        self.assertFalse(result["zip_member_contents_read"])
        self.assertNotIn("9999999",str(result))
    def test_mismatch_does_not_override_publisher_sha(self):
        raw=fake_zip()
        lock={"expected_file_bytes":len(raw),"dryad_sha256":"0"*64,
              "mirror_md5":hashlib.md5(raw).hexdigest()}
        with self.assertRaises(ValueError):verify_original_zip(raw,lock)
if __name__=="__main__":unittest.main()
