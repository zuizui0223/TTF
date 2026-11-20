"""Offline known-truth Dryad original-source access / hold tests (no internet)."""
import copy, json, sys, unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import dryad_original_source_access_preflight_v01 as a

class FreshSourceTests(unittest.TestCase):
    def setUp(self):self.c=json.loads(a.CONTRACT.read_text(encoding="utf-8"))
    def test_three_candidates_and_three_source_files(self):
        r=a.evaluate(self.c,runner=lambda file_id,name,geo: {
            "file":name,"geography_only":geo,"id":file_id,
            "status":"SYNTHETIC_SOURCE_ONLY","response_values_parsed_or_exported":False
        })
        self.assertEqual(r["study_candidates_before_response"],3)
        self.assertEqual([z["result"]["id"] for z in r["files"]],[4337649,4692223,4692224])
        self.assertFalse(r["biological_response_rows_parsed_or_modeled"])
    def test_original_experiment_header_schema_present(self):
        fields=["Garden","Position","Population","Rep","ContOrigin","ContGarden","Latitude","Longitude",
                "Cyanogenesis","Survived1st","SeededFull"]
        self.assertTrue(a.SOURCE_COLUMNS.issubset(a.inspect_csv_header((",".join(fields)+"\n").encode())))
    def test_missing_original_geo_header_does_not_qualify(self):
        with patch.object(a,"attempt",return_value=(
            {"status":"BYTES_RECEIVED"},b"Garden,Position,Population,Rep,ContOrigin\n"
        )):
            r=a.preflight_file(4337649,"Albano_et_al_2025_EM_Experiment_Data.csv")
        self.assertEqual(r["status"],"ORIGINAL_BYTES_OR_SOURCE_SCHEMA_NOT_CONFIRMED")
        self.assertFalse(r["original_source_GPS_independently_certified"])
    def test_source_only_small_location_coordinates(self):
        sample=b"site,latitude,longitude\nE01,34.3,-119.4\nE02,33.1,-117.5\n"
        x=a.source_geography_only_rows(sample,["site","latitude","longitude"])
        self.assertEqual(x["unique_location_pairs"],2)
        self.assertFalse(x["independent_original_source_GPS_certified"])
    def test_duplicate_coordinates_not_counted_as_independent_locations(self):
        sample=b"site,lat,lon\nA,33,-117\nB,33,-117\n"
        x=a.source_geography_only_rows(sample,["site","lat","lon"])
        self.assertEqual((x["location_row_count"],x["unique_location_pairs"]),(2,1))
    def test_garden_geography_file_not_promoted_to_source_GPS(self):
        with patch.object(a,"attempt",return_value=(
            {"status":"BYTES_RECEIVED"},b"site,latitude,longitude\nG,34,-119\n"
        )):
            r=a.preflight_file(4692224,"tplantlatlongs.csv",geography_only=True)
        self.assertTrue(r["location_file_summary"]["coordinates_parsed"])
        self.assertFalse(r["original_source_GPS_independently_certified"])
    def test_no_HTML_acceptance(self):
        with self.assertRaisesRegex(ValueError,"not CSV"):
            a.inspect_csv_header(b"<html>Login required</html>")
    def test_official_https_only(self):
        for url in ("https://attacker.example.com/a.csv","http://datadryad.org/a.csv"):
            with self.subTest(url=url):
                with self.assertRaisesRegex(ValueError,"official Dryad"):
                    a.trusted(url)
    def test_identity_drift_fails_closed(self):
        c=copy.deepcopy(self.c)
        c["candidates"][0]["study_id"]="ReplacementAfterFitness"
        with self.assertRaisesRegex(ValueError,"Candidate order"):
            a.evaluate(c,runner=lambda *args: {})
    def test_prior_frozen_gate_fails_closed(self):
        c=copy.deepcopy(self.c)
        c["old_bontrager_44_of_50_terminal_stop_unchanged"]=False
        with self.assertRaisesRegex(ValueError,"frozen science"):
            a.evaluate(c,runner=lambda *args: {})
    def test_official_file_ids_exact(self):
        self.assertEqual(a.dryad_file_urls(4692223),(
            "https://datadryad.org/api/v2/files/4692223/download",
            "https://datadryad.org/downloads/file_stream/4692223"))
    def test_location_coordinate_range_rejected(self):
        with self.assertRaisesRegex(ValueError,"out of range"):
            a.source_geography_only_rows(b"pop,lat,lon\nX,120,30\n",["pop","lat","lon"])

if __name__=="__main__":
    unittest.main()
