"""Frozen genuine-file Albano source-only receipt contract tests.

The 421,805-byte CSV was SHA-matched against the official Dryad digest
in the user container. GitHub Actions intentionally does not carry raw
plant data; it regression-checks this aggregate receipt and synthetic
source-only parser behavior. No biological response inference.
"""
import json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import albano2025_verified_archive_source_only_v01 as audit
HERE=Path(__file__).resolve().parent
RECEIPT=HERE/"source_registry/albano2025_verified_source_only_receipt_v01.json"
OFFICIAL=HERE/"source_registry/fresh_dryad_original_file_identity_v01.json"

class OriginalAlbanoReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=json.loads(RECEIPT.read_text("utf-8"))
        cls.source=json.loads(OFFICIAL.read_text("utf-8"))
    def test_exact_file_id_version_size_and_sha(self):
        albano=next(d for d in self.source["datasets"] if d["study_id"]=="Albano2025_Trifolium_repens")
        file=albano["files"][0]
        self.assertEqual(self.r["official_Dryad_version_id"],albano["version_id"])
        self.assertEqual(self.r["official_original_file_id"],file["id"])
        self.assertEqual(self.r["original_csv_path"],file["path"])
        self.assertEqual(self.r["original_csv_bytes"],file["size_bytes"])
        self.assertEqual(self.r["sha256_computed_from_uploaded_original_CSV"],file["sha256_from_repository_metadata"])
        self.assertEqual(self.r["sha256_computed_from_uploaded_original_CSV"],audit.ORIGINAL_SHA256)
    def test_96_source_population_47NAM_49EUR(self):
        self.assertEqual(self.r["original_source_population_IDs"],96)
        self.assertEqual(self.r["source_populations_by_continent"],{"EUR":49,"NAM":47})
        self.assertEqual(self.r["distinct_exact_coordinate_pairs"],96)
    def test_original_row_and_garden_counts(self):
        self.assertEqual(self.r["rows_original_experiment"],1858)
        self.assertEqual(sum(self.r["individual_rows_by_garden"].values()),1858)
        self.assertEqual(set(self.r["gardens"]),audit.EXPECTED_GARDENS)
    def test_exact_372_of_384_links_not_invented(self):
        self.assertEqual(self.r["population_garden_pairs_observed"],372)
        self.assertEqual(self.r["population_garden_pairs_full_cross_possible"],384)
        self.assertEqual(self.r["population_garden_pairs_absent_from_original_rows"],12)
        self.assertFalse(self.r["four_gardens_full_population_factorial"])
    def test_91_pop_all_gardens(self):
        cov=self.r["population_garden_coverage_frequency"]
        self.assertEqual(cov,{"1":2,"2":3,"4":91})
        self.assertEqual(sum(int(k)*v for k,v in cov.items()),372)
        self.assertEqual(sum(cov.values()),96)
    def test_composite_individual_key_only(self):
        self.assertTrue(self.r["individual_Garden_Position_unique"])
        self.assertTrue(self.r["Garden_Population_Rep_not_unique"])
        self.assertEqual(self.r["Garden_Population_Rep_colliding_key_groups"],107)
        self.assertEqual(self.r["Garden_Population_Rep_excess_rows_on_colliding_keys"],131)
    def test_source_coordinates_provenance_not_invented(self):
        self.assertTrue(self.r["source_id_to_reported_coordinate_mapping_verified_in_genuine_file"])
        self.assertTrue(self.r["official_repository_sha256_MATCHED"])
        self.assertFalse(self.r["actual_independently_surveyed_source_GPS_provenance_verified"])
        self.assertEqual(self.r["source_population_coordinate_conflicts"],0)
        self.assertEqual(self.r["source_population_continent_conflicts"],0)
    def test_nearby_population_pair_not_independent_lgm_cells(self):
        self.assertAlmostEqual(self.r["same_continent_nearest_source_pair_km"],0.388079,places=6)
        self.assertEqual(self.r["same_continent_site_pairs_under_1km"],1)
        self.assertFalse(self.r["all_96_sources_occupy_independent_21ka_climate_cells_proven"])
    def test_no_response_or_historical_effect(self):
        self.assertFalse(self.r["phenotype_fitness_genotype_or_LGM_values_modeled"])
        self.assertFalse(self.r["historical_21ka_exposure_obtained"])
        self.assertFalse(self.r["independent_LGM_fitness_effect_tested"])
        self.assertTrue(self.r["old_bontrager_44of50_stop_unchanged"])
    def test_only_source_columns_selected(self):
        self.assertEqual(tuple(self.r["source_geographic_provenance_columns_verified"]),
                         audit.SOURCE_FIELDS)
        self.assertNotIn("Survived1st",str(self.r))
    def test_no_original_coordinate_rows_committed(self):
        self.assertNotIn("sites",self.r)
        self.assertNotIn("latitude",self.r)
        self.assertNotIn("longitude",self.r)
        self.assertEqual(self.r["original_csv_columns_total"],35)
if __name__=="__main__":unittest.main()
