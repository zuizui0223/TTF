"""Offline official WorldClim paired LGM source URL and source-only tests."""
import copy,json,sys,unittest
from pathlib import Path
from urllib.error import HTTPError
sys.path.insert(0,str(Path(__file__).resolve().parent))
import worldclim_lgm_source_transport_preflight_v01 as a

class WorldClimPaleoSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract=json.loads(a.CONTRACT.read_text("utf-8"))
    def test_four_fixed_official_archives(self):
        x=a.probe(self.contract,runner=lambda url:{
            "source_head_classification":"ZIP_LOCAL_HEADER_OBSERVED_NOT_FULL_FILE_VERIFIED"})
        self.assertEqual(set(x["responses"]),{"CCSM4","MIROC-ESM","MPI-ESM-P","WC1.4_MODERN"})
        self.assertTrue(all(t["status"]=="OFFICIAL_ZIP_HEADER_PRESENT" for t in x["responses"].values()))
    def test_zip_head_not_full_integrity(self):
        self.assertEqual(a.classify_head(b"PK\x03\x04hello"),
                         "ZIP_LOCAL_HEADER_OBSERVED_NOT_FULL_FILE_VERIFIED")
    def test_html_not_valid_archive(self):
        self.assertEqual(a.classify_head(b"<html>Forbidden</html>"),"HTML_NOT_ZIP")
        self.assertEqual(a.classify_head(b""),"EMPTY_RESPONSE")
    def test_http_access_failure_gated(self):
        def block(url):
            raise HTTPError(url,403,"Forbidden",{},None)
        result=a.probe(self.contract,runner=block)
        self.assertEqual(len(result["responses"]),4)
        self.assertEqual({v["http_error"] for v in result["responses"].values()},{403})
    def test_gcm_url_mutation_rejected(self):
        d=copy.deepcopy(self.contract)
        d["models"][0]["official_url"]="https://not-worldclim.example.org/test.zip"
        with self.assertRaisesRegex(ValueError,"source URL"):a.probe(d,runner=lambda u:{})
    def test_climate_era_not_worldclim_v2_1970_2000(self):
        self.assertEqual(self.contract["historical_product"].split()[0],"WorldClim")
        self.assertEqual(self.contract["matched_baseline_product"].split()[0],"WorldClim")
        self.assertIn("1960-1990",self.contract["matched_baseline_product"])
    def test_21ka_not_mislabeled_as_precise_date(self):
        d=a.probe(self.contract,runner=lambda url:{"source_head_classification":"UNRECOGNIZED_CONTENT"})
        self.assertIn("not exactly 21000",d["LGM_snapshot_period"])
    def test_no_LGM_values_or_fitness_in_access_preflight(self):
        d=a.probe(self.contract,runner=lambda url:{"source_head_classification":"UNRECOGNIZED_CONTENT"})
        self.assertFalse(d["any_original_climate_zip_fully_downloaded"])
        self.assertFalse(d["any_LGM_or_current_bioclim_raster_values_extracted"])
        self.assertFalse(d["biological_response_opened_or_tested"])
        self.assertFalse(d["actual_source_GPS_records_submitted_to_CI"])
    def test_frozen_continent_and_location_grain(self):
        self.assertEqual(self.contract["source_population_count"],96)
        self.assertEqual(self.contract["nominal_unique_land_grid_cells_at_2_5_arcminutes"],95)
    def test_fixed_primary_and_sensitivity_GCM_order(self):
        self.assertEqual([(x["role"],x["code"]) for x in self.contract["models"]],
            [("primary","CCSM4"),("sensitivity","MIROC-ESM"),("sensitivity","MPI-ESM-P")])
if __name__=="__main__":unittest.main()
