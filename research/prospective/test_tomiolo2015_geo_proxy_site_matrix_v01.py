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
    def test_frozen_original_action_artifact_persisted_byte_identical(self):
        import hashlib
        original=HERE/"source_registry/tomiolo2015_geo_proxy_matrix_v01.json"
        self.assertEqual(hashlib.sha256(original.read_bytes()).hexdigest(),
            "78f1096d0e3ba985a6bfc55f2d7f2911f8b8187322b57ef0243fa2c35bf8e98c")
        data=json.loads(original.read_text(encoding="utf-8"))
        self.assertEqual(data["original_study_source_population_ID_keys"],6)
        self.assertEqual(data["original_study_unique_garden_derived_site_points"],2)
        self.assertEqual(data["strict_binomial_species_with_two_distinct_garden_proxy_locations"],3)
        self.assertFalse(data["source_population_GPS_directly_verified"])
    def test_source_manifest_three_appendices_are_not_three_origins(self):
        r=json.loads((HERE/"source_registry/tomiolo2015_figshare_appendixA_source_provenance_v01.json").read_text())
        self.assertEqual(r["original_paper"]["official_collection_id"],3307653)
        self.assertEqual(len(r["original_paper"]["original_Figshare_articles"]),3)
        self.assertEqual(r["original_paper"]["field_origin_site_names"],["M","SA"])
        self.assertEqual(r["original_paper"]["appendix_gradient_context_only_site"],"A")
        self.assertEqual(r["reference_full_original"]["Hargreaves_SHA256"],
                         self.c["source"]["Hargreaves_exact_sha256"])
        self.assertEqual(r["reference_full_original"]["original_geo_artifact_SHA256"],
                         "78f1096d0e3ba985a6bfc55f2d7f2911f8b8187322b57ef0243fa2c35bf8e98c")
        self.assertEqual(r["reference_full_original"]["population_ID_keys"],6)
        self.assertEqual(r["reference_full_original"]["distinct_source_derived_garden_proxy_locations"],2)
        self.assertEqual(r["reference_full_original"]["same_original_study_origin_GPS_independently_confirmed"],0)
    def test_not_independent_GPS_or_new_fitness_adaptation_support(self):
        r=json.loads((HERE/"source_registry/tomiolo2015_figshare_appendixA_source_provenance_v01.json").read_text())
        self.assertTrue(r["descriptive_geo_distance"]["this_is_NOT_independent_original_field_GPS_measurement"])
        self.assertTrue(r["qualifications"]["original_Appendix_A_source_context_contained_qualitative_mortality_language_in_initial_excerpt"])
        self.assertTrue(r["qualifications"]["no_numerical_fitness_response_data_modeled"])
        self.assertTrue(r["qualifications"]["no_LGM_21ka_0BP_raster_sampled"])
        self.assertTrue(r["qualifications"]["Bontrager_previous_44of50_gate_unchanged"])
    def test_hargreaves35_source_only_firstwave_not_19_certified_experiments(self):
        import hashlib
        root=HERE.parents[1]
        rec=json.loads((HERE/"source_registry/hargreaves35_independent_source_firstwave_eligibility_v01.json").read_text())
        inp=rec["input_original_source_evidence"]
        orig1=(root/inp["hargreaves_source_only_study_geometry_json"]).read_bytes()
        orig2=(root/inp["hargreaves_source_as_garden_proxy_json"]).read_bytes()
        self.assertEqual(hashlib.sha256(orig1).hexdigest(),inp["source_only_geometry_original_SHA256"])
        self.assertEqual(hashlib.sha256(orig2).hexdigest(),inp["source_garden_geo_original_SHA256"])
        old1=json.loads(orig1);old2=json.loads(orig2)
        ids=old2["per_unlinked_study"]
        qualified=[v for v in ids if v["multi_source_cells_with_two_distinct_georecovered_coordinates"]>0]
        duplicate={"Raabovaetal_BiolConservation_2007","PankovaandRaabova_PLoSOne_2014"}
        pending=[v for v in qualified if v["original_study_label"] not in duplicate]
        a=rec["source_only_counts"]
        self.assertEqual(a["total_original_unlinked_study_labels"],35)
        self.assertEqual(a["with_at_least_two_distinct_garden_derived_origin_proxy_coordinates"],len(qualified))
        self.assertEqual(a["original_geo_proxy_candidate_cells"],sum(v["multi_source_cells_with_two_distinct_georecovered_coordinates"] for v in qualified))
        self.assertEqual(a["remaining_non_Raabova_study_labels_with_two_proxy_locations"],len(pending))
        self.assertEqual(a["remaining_non_Raabova_comparison_cells"],sum(v["multi_source_cells_with_two_distinct_georecovered_coordinates"] for v in pending))
        self.assertEqual((len(qualified),len(pending)),(21,19))
        self.assertEqual(a["remaining_non_Raabova_comparison_cells"],114)
        self.assertEqual(a["independently_certified_original_experiment_cohorts"],0)
        self.assertEqual(a["independently_certified_source_population_GPS_studies"],0)
        self.assertTrue(rec["science_boundary"]["not_19_independently_certified_trials"])
        self.assertTrue(rec["science_boundary"]["Bontrager_44_of_50_source_gate_frozen"])
    def test_new_orchid_study_is_biotic_endpoint_not_inferred_climate_adaptation(self):
        rec=json.loads((HERE/"source_registry/hargreaves35_independent_source_firstwave_eligibility_v01.json").read_text())
        sun=rec["known_study_specific_source_geography"]["Sunetal_AnnBot_2014"]
        self.assertEqual(sun["original_DOI"],"10.1093/aob/mct219")
        self.assertIn("pollinator-guild",sun["original_research_scope"])
        self.assertEqual(sun["unique_proxy_points"],9)
        self.assertFalse(sun["source_population_original_GPS_qualified"])
    def test_prior_source_and_genetic_stops_remain_unchanged(self):
        x=proxy_geography(self.original_metadata(),self.c)
        self.assertTrue(x["old_Bontrager_44of50_source_gate_unchanged"])
        self.assertFalse(x["independent_primary_experiment_cohort_verified"])
        self.assertFalse(x["Hargreaves_garden_proxy_coords_exported"])
if __name__=="__main__":unittest.main()
