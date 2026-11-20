"""Offline known-truth / frozen original-receipt tests for Albano modern covariates.

These tests intentionally DO NOT claim that GitHub Actions reruns the
user-supplied ZIP. Reproduce with --zip locally when the original is present.
"""
import inspect,json,math,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import albano2025_modern_climate_difference_integrity_v01 as a
FROZEN=ROOT/"source_registry/albano2025_modern_climate_integrity_combined_v01.json"
OFFICIAL=ROOT/"source_registry/fresh_dryad_original_file_identity_v01.json"

class ModernDifferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.f=json.loads(FROZEN.read_text("utf-8"))
        cls.d=json.loads(OFFICIAL.read_text("utf-8"))
    def test_original_sha_from_official_manifest(self):
        origin=next(v for v in self.d["datasets"] if v["study_id"].startswith("Albano"))
        self.assertEqual(a.HASHES[a.RAW],origin["files"][0]["sha256_from_repository_metadata"])
        self.assertEqual(self.f["original_csv_sha256"],origin["files"][0]["sha256_from_repository_metadata"])
    def test_six_bioclim_variables(self):
        self.assertEqual(a.BIO,("BIO1","BIO5","BIO6","BIO12","BIO13","BIO14"))
        self.assertEqual(len(self.f["Montpellier_omitted_garden_z_vector"]),6)
    def test_correct_source_minus_garden(self):
        self.assertTrue(all(abs(x)<1e-14 for x in a.modern_error(
            [10,8,6,4,2,0],[1,2,3,4,5,6],[9,6,3,0,-3,-6])))
    def test_source_only_input_missing_garden_subtraction(self):
        src=[2,3,4,5,6,7]
        garden=[0.81,0.5,0.92,-0.95,-0.31,-1.44]
        self.assertTrue(all(abs(x-y)<1e-12 for x,y in zip(
            a.modern_error(src,garden,src),garden)))
    def test_frozen_470_rows_and_2820_values(self):
        f=self.f
        self.assertEqual(f["Montpellier_rows_with_missing_garden_subtraction"],470)
        self.assertEqual(f["Montpellier_bioclim_values_affected"],470*6)
        self.assertEqual(f["non_Montpellier_rows_with_correct_6_signed_diff_fields"],1388)
        self.assertEqual(f["total_individuals"],1858)
    def test_three_gardens_correct_one_offset(self):
        d=self.f["garden_rows_and_diff_integrity"]
        self.assertEqual(set(d),set(a.GARDENS))
        self.assertTrue(all(d[g]["correct"] for g in a.GARDENS if g!="Montpellier"))
        self.assertFalse(d["Montpellier"]["correct"])
        self.assertEqual(d["Montpellier"]["max_abs_error"],1.447067214)
    def test_94_cross_garden_source_checks_two_readme_only(self):
        self.assertEqual(self.f["source_rows_independently_matched_via_other_gardens"],94)
        self.assertEqual(self.f["Montpellier_only_source_rows_supported_by_readme_order"],2)
        self.assertEqual(self.f["source_populations"],96)
    def test_pca_weights_known_truth(self):
        pc=[[1,2,3,4,5,6],[2,4,6,8,10,12]]
        expected=[1/91,4/91,9/91]
        self.assertTrue(all(abs(x-y)<1e-12 for x,y in zip(a.pca_weights(pc),expected)))
    def test_weighted_euclidean_distance_known_truth(self):
        self.assertAlmostEqual(a.weighted_distance([2,3,4],[1,1,1],[0.5,0.3,0.2]),
                               math.sqrt(0.5*1+0.3*4+0.2*9),12)
    def test_weightedCD_all_original_rows_reproduced(self):
        self.assertEqual(self.f["WeightedCD_reproduced_rows"],1858)
        self.assertLess(self.f["WeightedCD_max_abs_error"],2e-9)
        self.assertAlmostEqual(self.f["WeightedCD_pc_explained_weight_sum"],0.970009806,8)
        self.assertFalse(self.f["fitness_genetics_or_biological_response_modeled"])
    def test_different_bio_shifts_are_absorbed_by_linear_garden_intercepts(self):
        # The same prediction has different intercept+linear coefficients.
        # For quadratic x: alpha+beta*x+gamma*x^2, x=y+c.
        for x,c in [(-2.5,0.811598694),(0.0,-1.447067214),(1.4,0.5)]:
            alpha,beta,gamma=1.0,0.7,-0.4
            y=x-c
            lhs=alpha+beta*x+gamma*x*x
            rhs=(alpha+beta*c+gamma*c*c)+(beta+2*gamma*c)*y+gamma*y*y
            self.assertAlmostEqual(lhs,rhs,12)
    def test_nominal_grid_95_cells_not_96_independent_lgm_exposures(self):
        self.assertEqual(self.f["nominal_worldclim_2_5m_cells_from_geometry"],95)
        self.assertEqual(self.f["nominal_worldclim_2_5m_duplicate_cell_groups"],1)
        self.assertFalse(self.f["actual_LGM_climate_values_extracted"])
    def test_current_signal_not_historical_or_fitness(self):
        for flag in ("official_LGM_rasters_downloaded","history_incremental_predictor_tested",
                     "fitness_genetics_or_biological_response_modeled",
                     "new_confirmatory_LGM_effect_asserted"):
            self.assertFalse(self.f[flag])
        self.assertTrue(self.f["stopped_44of50_Bontrager_unchanged"])
    def test_original_data_zip_not_part_of_public_repository(self):
        self.assertFalse((ROOT/"doi_10_5061_dryad_bk3j9kdqk__v20250905.zip").exists())
        self.assertFalse((ROOT/"Albano_et_al_2025_EM_Experiment_Data.csv").exists())
if __name__=="__main__":
    unittest.main()
