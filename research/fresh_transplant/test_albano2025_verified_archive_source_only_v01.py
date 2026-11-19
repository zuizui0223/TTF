"""Synthetic fixture tests for genuine-file Albano source-only ZIP audit.

No original biological responses are present. Real uploaded ZIP is NOT
committed to GitHub; it was hashed and audited in the user container.
"""
import csv,io,sys,unittest,zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import albano2025_verified_archive_source_only_v01 as audit

HEADERS=list(audit.SOURCE_FIELDS)+["Survived1st","SeededFull","Ac_ac"]
def record(garden,position,pop,rep,origin,continent,lat,lon):
    return [garden,str(position),str(pop),str(rep),origin,continent,
            str(lat),str(lon),"1","0","Ac"]

def zip_fixture(rows,headers=None,name=None):
    headers=headers or HEADERS
    stream=io.StringIO(newline="")
    w=csv.writer(stream)
    w.writerow(headers)
    w.writerows(rows)
    buff=io.BytesIO()
    with zipfile.ZipFile(buff,"w",compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr(name or audit.ORIGINAL_CSV,stream.getvalue().encode("utf-8"))
    buff.seek(0)
    return buff

BASE=[
    record("Lafayette",1,1,1,"NAM","America",30.0,-90.0),
    record("Lafayette",2,1,1,"NAM","America",30.0,-90.0),
    record("Mississauga",1,1,1,"NAM","America",30.0,-90.0),
    record("Lafayette",3,2,1,"NAM","America",30.001,-90.0),
    record("Montpellier",1,56,1,"EUR","Europe",48.0,2.0),
]

class AlbanoSourceOnlyTests(unittest.TestCase):
    def test_whitelisted_source_schema_and_unique_individual(self):
        r,s=audit.audit_zip(zip_fixture(BASE),verify_pinned=False)
        self.assertEqual(r["rows_original_experiment"],5)
        self.assertEqual(r["original_source_population_IDs"],3)
        self.assertEqual(r["distinct_exact_coordinate_pairs"],3)
        self.assertTrue(r["individual_Garden_Position_unique"])
        self.assertFalse(r["phenotype_fitness_genotype_or_LGM_values_modeled"])
        self.assertFalse(r["source_id_to_reported_coordinate_mapping_verified_in_genuine_file"])
        self.assertEqual(s["site_count"],3)

    def test_rep_not_unique_is_detected_but_no_invented_individuals(self):
        r,_=audit.audit_zip(zip_fixture(BASE),verify_pinned=False)
        self.assertEqual(r["Garden_Population_Rep_colliding_key_groups"],1)
        self.assertEqual(r["Garden_Population_Rep_excess_rows_on_colliding_keys"],1)

    def test_incomplete_population_garden_design_not_filled(self):
        r,_=audit.audit_zip(zip_fixture(BASE),verify_pinned=False)
        self.assertEqual(r["population_garden_pairs_observed"],4)
        self.assertEqual(r["population_garden_pairs_full_cross_possible"],9)
        self.assertEqual(r["population_garden_pairs_absent_from_original_rows"],5)
        self.assertFalse(r["four_gardens_full_population_factorial"])

    def test_near_geographic_locations_detected(self):
        r,_=audit.audit_zip(zip_fixture(BASE),verify_pinned=False)
        self.assertEqual(r["same_continent_site_pairs_under_1km"],1)
        self.assertAlmostEqual(r["same_continent_nearest_source_pair_km"],0.111195,places=5)

    def test_bad_csv_sha_cannot_authenticate(self):
        with self.assertRaisesRegex(ValueError,"SHA256"):
            audit.audit_zip(zip_fixture(BASE),verify_pinned=True)

    def test_missing_required_source_column_fails(self):
        headers=[h for h in HEADERS if h!="Longitude"]
        rows=[[value for i,value in enumerate(row) if HEADERS[i]!="Longitude"] for row in BASE]
        with self.assertRaisesRegex(ValueError,"field names"):
            audit.audit_zip(zip_fixture(rows,headers),verify_pinned=False)

    def test_population_coordinate_drift_fails(self):
        rows=[row[:] for row in BASE]
        rows[1][HEADERS.index("Latitude")]="31.0"
        with self.assertRaisesRegex(ValueError,"location/origin drift"):
            audit.audit_zip(zip_fixture(rows),verify_pinned=False)

    def test_population_origin_drift_fails(self):
        rows=[row[:] for row in BASE]
        rows[2][HEADERS.index("ContOrigin")]="EUR"
        with self.assertRaisesRegex(ValueError,"location/origin drift"):
            audit.audit_zip(zip_fixture(rows),verify_pinned=False)

    def test_garden_continent_drift_fails(self):
        rows=[row[:] for row in BASE]
        rows[1][HEADERS.index("ContGarden")]="Europe"
        with self.assertRaisesRegex(ValueError,"Garden-to-continent provenance"):
            audit.audit_zip(zip_fixture(rows),verify_pinned=False)

    def test_garden_position_duplicate_fails(self):
        rows=[row[:] for row in BASE]
        rows[1][HEADERS.index("Position")]=rows[0][HEADERS.index("Position")]
        with self.assertRaisesRegex(ValueError,"not a unique"):
            audit.audit_zip(zip_fixture(rows),verify_pinned=False)

    def test_invalid_or_nan_geography_fails(self):
        for invalid in ["120.0","nan",""]:
            rows=[row[:] for row in BASE]
            rows[0][HEADERS.index("Latitude")]=invalid
            with self.subTest(invalid=invalid):
                with self.assertRaisesRegex(ValueError,"coordinate|Missing source"):
                    audit.audit_zip(zip_fixture(rows),verify_pinned=False)

    def test_unauthenticated_original_filename_fails(self):
        with self.assertRaisesRegex(ValueError,"original CSV missing"):
            audit.audit_zip(zip_fixture(BASE,name="fake_outcomes.csv"),verify_pinned=False)

    def test_site_registry_never_exports_response_values(self):
        r,sites=audit.audit_zip(zip_fixture(BASE),verify_pinned=False)
        x=str((r,sites))
        self.assertNotIn("Survived1st",x)
        self.assertNotIn("SeededFull",x)
        self.assertNotIn("Ac_ac",x)
        self.assertFalse(r["historical_21ka_exposure_obtained"])

if __name__=="__main__":unittest.main()
