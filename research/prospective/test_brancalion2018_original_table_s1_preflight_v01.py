"""Offline tests for official Brancalion 2018 Table S1 source-only recovery."""
import io,sys,tarfile,unittest,zipfile
from pathlib import Path
from urllib.error import HTTPError
sys.path.insert(0,str(Path(__file__).resolve().parent))
import brancalion2018_original_table_s1_preflight_v01 as audit

def synthetic_docx():
    body=('<?xml version="1.0"?><w:document '
          'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
          '<w:body><w:p><w:r><w:t>Caetetus Carlos Botelho Ilha do Cardoso'
          ' latitude longitude mother -22.3835933 -49.6793777</w:t></w:r></w:p>'
          '</w:body></w:document>').encode()
    obj=io.BytesIO()
    with zipfile.ZipFile(obj,"w") as z:
        z.writestr("[Content_Types].xml","<Types/>")
        z.writestr("word/document.xml",body)
    return obj.getvalue()

class BrancalionTableS1Tests(unittest.TestCase):
    def test_docx_metadata_source_only(self):
        x=audit.docx_metadata(synthetic_docx())
        self.assertEqual(len(x["original_site_names_mentioned"]),3)
        self.assertGreaterEqual(x["geographic_or_parentage_field_keyword_count"],3)
        self.assertEqual(x["coordinate_shaped_text_count"],1)
        self.assertFalse(x["document_text_extracted_to_receipt"])
        self.assertNotIn("latitude longitude",str(x))

    def test_fake_publisher_success_never_promotes_gps(self):
        x=audit.probe(fetcher=lambda url,n: synthetic_docx())
        self.assertIn("STORAGE_CONFIRMED",x["status"])
        self.assertEqual(x["independent_original_maternal_tree_GPS_certified"],0)
        self.assertEqual(x["independent_source_population_GPS_certified"],0)
        self.assertFalse(x["historical_climate_contrast_qualified"])

    def test_no_public_access_yields_honest_hold(self):
        def fail(url,n):
            raise HTTPError(url,403,"Forbidden",{},None)
        x=audit.probe(fetcher=fail)
        self.assertEqual(x["status"],"ORIGINAL_TABLE_S1_BYTES_NOT_RECOVERED")
        self.assertEqual(len(x["attempts"]),4)
        self.assertTrue(all(a["http_code"]==403 for a in x["attempts"]))
        self.assertFalse(x["fitness_SNP_LGM_adaptation_plasticity_speciation_response_opened"])

    def test_official_oa_api_url_validation(self):
        xml=b'<records><record id="PMC6106193"><link format="tgz" href="ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/test.tgz"/></record></records>'
        self.assertEqual(audit.oai_archive_url(xml),
                         "https://ftp.ncbi.nlm.nih.gov/pub/pmc/test.tgz")

    def test_unrelated_pmc_record_rejected(self):
        xml=b'<records><record id="PMC000000"><link format="tgz" href="ftp://ftp.ncbi.nlm.nih.gov/test.tgz"/></record></records>'
        self.assertIsNone(audit.oai_archive_url(xml))

    def test_official_oa_archive_success(self):
        raw=synthetic_docx();tarbuf=io.BytesIO()
        with tarfile.open(fileobj=tarbuf,mode="w:gz") as z:
            header=tarfile.TarInfo("PMC6106193/"+audit.SUPPLEMENT)
            header.size=len(raw);z.addfile(header,io.BytesIO(raw))
        self.assertEqual(audit.select_from_oa_archive(tarbuf.getvalue()),raw)

    def test_unapproved_host_rejected_before_network(self):
        with self.assertRaisesRegex(ValueError,"approved official HTTPS"):
            audit.fetch("https://example.invalid/evil.docx",1024)

    def test_not_a_docx_rejected(self):
        with self.assertRaisesRegex(ValueError,"DOCX ZIP"):
            audit.docx_metadata(b"<html>login</html>")

    def test_no_false_independent_studies_from_new_bytes(self):
        x=audit.probe(fetcher=lambda url,n: synthetic_docx())
        self.assertFalse(x["cross_assay_maternal_ID_join_completed"])
        self.assertTrue(x["prior_Bontrager_44_of_50_stop_unchanged"])

if __name__=="__main__":
    unittest.main()
