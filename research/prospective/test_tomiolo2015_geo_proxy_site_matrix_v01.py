"""Synthetic frozen source-only two Hargreaves proxy sites vs three pub sites."""
import copy,json,sys,unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from tomiolo2015_geo_proxy_site_matrix_v01 import CONTRACT,proxy_geography

class TomioloTwoSiteGeographyTests(unittest.TestCase):
    def setUp(self):
        self.c=copy.deepcopy(json.loads(CONTRACT.read_text(encoding="utf-8")))
        self.study=self.c["source"]["Hargreaves_original_study_label"]
    def original_metadata(self):
        keys=set()
        sites={}
        for taxon in ("Biscutella didyma","Brachypodium distachyon","Stipa capensis"):
            for source,xy in (("M",(35.05,31.7)),("SA",(34.9,31.3833333333))):
                key=(self.study,taxon.lower(),f"{taxon}-{source}")
                keys.add(key);sites[key]={xy}
        return {"origin_keys":keys,"gardens":sites}
    def test_six_population_keys_are_two_not_six_site_locations(self):
        x=proxy_geography(self.original_metadata(),self.c)
        self.assertEqual(x["original_study_source_population_ID_keys"],6)
        self.assertEqual(x["original_study_unique_garden_derived_site_points"],2)
        self.assertEqual(x["original_study_binomial_species_labels"],3)
        self.assertEqual(x["strict_binomial_species_with_two_distinct_garden_proxy_locations"],3)
        self.assertFalse(x["source_population_GPS_directly_verified"])
        self.assertFalse(x["study_site_code_to_source_population_ID_directly_certified"])
        self.assertTrue(x["no_fitness_phenotype_genotype_LGM_model_opened"])
    def test_full_three_published_site_distances_are_descriptive_only(self):
        x=proxy_geography(self.original_metadata(),self.c)
        matrix=x["anonymous_proxy_point_to_published_site_distance_km"]
        self.assertEqual(len(matrix),2)
        self.assertEqual(set(matrix["anonymous_proxy_point_1"]),{"M","SA","A"})
        self.assertEqual(matrix["anonymous_proxy_point_1"]["SA"],0)
        self.assertEqual(matrix["anonymous_proxy_point_2"]["M"],0)
        self.assertTrue(x["A_is_reference_aridity_gradient_site_not_automatically_an_experiment_source"])
    def test_missing_garden_proxy_never_guess_site_by_paper_coordinates(self):
        m=self.original_metadata()
        k=next(iter(m["origin_keys"]))
        del m["gardens"][k]
        with self.assertRaisesRegex(ValueError,"not also"):
            proxy_geography(m,self.c)
    def test_two_different_sites_with_one_malformed_point_rejected(self):
        m=self.original_metadata()
        k=next(iter(m["origin_keys"]))
        m["gardens"][k].add(None)
        with self.assertRaisesRegex(ValueError,"missing or ambiguous"):
            proxy_geography(m,self.c)
    def test_wrong_site_cardinality_does_not_admit_third_dry_gradient(self):
        m=self.original_metadata()
        k=next(iter(m["origin_keys"]))
        m["gardens"][k]={(34.7667,30.8667)}
        with self.assertRaisesRegex(ValueError,"geometry drifted"):
            proxy_geography(m,self.c)
    def test_no_adjustment_of_distance_band_after_source_results(self):
        self.c["matching_rules"]["descriptive_distance_bands_km"]=[1,2,5,20]
        with self.assertRaisesRegex(ValueError,"no-retune"):
            proxy_geography(self.original_metadata(),self.c)
    def test_prior_source_and_genetic_stops_remain_unchanged(self):
        x=proxy_geography(self.original_metadata(),self.c)
        self.assertTrue(x["old_Bontrager_44of50_source_gate_unchanged"])
        self.assertFalse(x["independent_primary_experiment_cohort_verified"])
        self.assertFalse(x["Hargreaves_garden_proxy_coords_exported"])
if __name__=="__main__":unittest.main()
