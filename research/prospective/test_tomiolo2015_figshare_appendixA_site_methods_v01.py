"""Known-truth source-only Figshare Appendix A original identity and location snippets."""
import copy
import hashlib
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from tomiolo2015_figshare_appendixA_site_methods_v01 import (
    CONTRACT,verified_file_url,source_site_only_excerpts,probe
)
class TomioloOriginalAppendixASafetyTests(unittest.TestCase):
    def setUp(self):
        self.c=copy.deepcopy(json.loads(CONTRACT.read_text(encoding="utf-8")))
        self.src=self.c["source"]
    def article(self,which=None):
        f={"id":self.src["original_file_id"],"name":which or self.src["original_filename"],
           "size":self.src["published_bytes"],
           "computed_md5":self.src["published_md5"],
           "download_url":"https://ndownloader.figshare.com/files/"+str(self.src["original_file_id"])}
        return {"id":self.src["article_id"],"doi":self.src["article_doi"],"files":[f]}
    def test_exact_original_A_metadata_approved(self):
        url=verified_file_url(self.article(),self.c)
        self.assertTrue(url.startswith("https://ndownloader.figshare.com/files/"))
    def test_wrong_original_Appendix_B_or_DOI_never_approved(self):
        with self.assertRaisesRegex(ValueError,"Appendix A file absent"):
            verified_file_url(self.article("appendix-B.htm"),self.c)
        bad=self.article();bad["doi"]="10.6084/m9.figshare.3562071.v1"
        with self.assertRaisesRegex(ValueError,"same frozen"):
            verified_file_url(bad,self.c)
    def test_published_MD5_and_length_gate(self):
        bad=self.article()
        bad["files"][0]["computed_md5"]="00"*16
        with self.assertRaisesRegex(ValueError,"MD5 metadata"):
            verified_file_url(bad,self.c)
    def test_site_method_only_html_excerpts_not_survival_or_biomass(self):
        html=(b"<html><style>ignore this text</style><p>Study sites in Israel: Mediterranean M Matta"
              b" and semi-arid SA Lahav, both visited during the field trial. "+b" calcareous landscape"*25+b"</p>"
              b"<p>Survival data were 18 and performance was 5 in a site.</p></html>")
        self.src["published_bytes"]=len(html)
        self.src["published_md5"]=hashlib.md5(html).hexdigest()
        result=source_site_only_excerpts(html,self.c)
        self.assertGreater(result["text_context_count"],0)
        joined=json.dumps(result).lower()
        self.assertNotIn("survival data",joined)
        self.assertTrue(result["numeric_fitness_or_experimental_response_not_extracted"])
    def test_unverified_bytes_rejected_before_any_text_extraction(self):
        html=b"<p>Matta M site</p>"
        with self.assertRaisesRegex(ValueError,"Unverified original"):
            source_site_only_excerpts(html,self.c)
    def test_403_is_unavailable_not_population_gps_negative(self):
        def denied(req,timeout=30):
            raise HTTPError(req.full_url,403,"Forbidden",None,None)
        with patch("tomiolo2015_figshare_appendixA_site_methods_v01.urlopen",side_effect=denied):
            out=probe(self.c)
        self.assertEqual(out["status"],"OFFICIAL_APPENDIX_A_ACCESS_OR_IDENTITY_NOT_EVALUABLE")
        self.assertFalse(out["original_file_downloaded"])
        self.assertFalse(out["source_site_population_ID_crosswalk_certified"])
    def test_unfrozen_appendixB_access_not_allowed(self):
        self.c["authorization"]["no_appendixB_or_appendixC_fetch"]=False
        with self.assertRaisesRegex(ValueError,"contract drift"):
            probe(self.c)

if __name__=="__main__":
    unittest.main()
