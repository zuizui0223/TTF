"""Known-truth six-source, six-garden three-basin habitat design source-only tests."""
import copy
import json
import sys
import unittest
from collections import defaultdict
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from ellis2007_six_sites_basin_habitat_geometry_v01 import (
    CONTRACT, geometry, source_id_category
)

class EllisOriginalBasinQuartzDesignTests(unittest.TestCase):
    def setUp(self):
        self.c=copy.deepcopy(json.loads(CONTRACT.read_text(encoding="utf-8")))
        self.study=self.c["study"]["Hargreaves_original_study_label"]
        self.taxon=self.c["study"]["original_organism"]
    def metadata(self):
        names=("Ariz","Sout","QK")
        positions={}
        for i,b in enumerate(names):
            for hab,lat_add in (("dense",0),("sparse",0.006)):
                suffix=b+"_"+hab
                positions[suffix]=(18.41+0.047*i,-31.1+lat_add)
        origin={(self.study,self.taxon,s) for s in positions}
        gardens=defaultdict(set)
        for suffix,coord in positions.items():
            gardens[(self.study,self.taxon,suffix)].add(coord)
        cells={}
        for target in positions:
            basin,hab=target.split("_")
            for metric in ("size","survival","reproduction"):
                # one local origin + one within basin other habitat +
                # two same habitat other basins, 4 sources per garden
                sources={s for s in positions if
                         s.split("_")[0]==basin or s.split("_")[1]==hab}
                cells[self.study,self.taxon,"adult2000",target,metric]=sources
        return {"origin_keys":origin,"gardens":gardens,"source_groups":cells}
    def test_original_3_basins_2_quartz_habitats_six_locations(self):
        r=geometry(self.metadata(),self.c)
        self.assertEqual(r["source_population_ID_keys"],6)
        self.assertEqual(r["unique_garden_derived_proxy_location_points"],6)
        self.assertTrue(r["original_study_basin_habitat_labels_directly_classifiable_from_ID_text"])
        self.assertEqual(r["all_15_unique_source_point_pair_distances_km"]["count"],15)
        self.assertEqual(r["within_basin_3_pairs_if_literal_IDs_qualified"]["count"],3)
        self.assertEqual(r["between_basin_12_pairs_if_literal_IDs_qualified"]["count"],12)
        self.assertEqual(r["original_distinct_source_plus_garden_directed_ID_pairs"],24)
        self.assertEqual(r["per_garden_distinct_source_ID_count_distribution"],{"4":6})
        self.assertEqual(r["original_source_metadata_comparison_cells_not_independent"],18)
        self.assertEqual(r["between_basin_different_quartz_habitat_edges_observed_if_codes_qualified"],0)
        self.assertFalse(r["no_21ka_0BP_climate_predictor_or_fitness_effect_fitted"] is False)
    def test_original_basin_and_habitat_names_cannot_be_inferred_from_positions(self):
        m=self.metadata()
        transformed={}
        coords={}
        mapping={x[2]:"unknownID"+str(i) for i,x in enumerate(sorted(m["origin_keys"]))}
        for key in m["origin_keys"]:
            transformed[key[:2]+(mapping[key[2]],)]=m["gardens"][key]
        newkeys=set(transformed)
        newcells={}
        for (study,taxon,exp,garden,metric),sources in m["source_groups"].items():
            newcells[(study,taxon,exp,mapping[garden],metric)]={mapping[s] for s in sources}
        revised={"origin_keys":newkeys,"gardens":transformed,"source_groups":newcells}
        r=geometry(revised,self.c)
        self.assertEqual(r["all_15_unique_source_point_pair_distances_km"]["count"],15)
        self.assertFalse(r["original_study_basin_habitat_labels_directly_classifiable_from_ID_text"])
        self.assertIsNone(r["within_basin_3_pairs_if_literal_IDs_qualified"])
        self.assertIsNone(r["between_basin_12_pairs_if_literal_IDs_qualified"])
    def test_source_has_two_different_origins_with_one_same_coordinate_rejected(self):
        m=self.metadata()
        a,b=sorted(m["origin_keys"])[:2]
        m["gardens"][b]=set(m["gardens"][a])
        with self.assertRaisesRegex(ValueError,"six distinct"):
            geometry(m,self.c)
    def test_original_source_ID_not_exact_garden_no_coordinate_imputation(self):
        m=self.metadata()
        k=sorted(m["origin_keys"])[0]
        m["gardens"].pop(k)
        with self.assertRaisesRegex(ValueError,"not an exact"):
            geometry(m,self.c)
    def test_two_garden_coordinates_for_same_source_abstain(self):
        m=self.metadata()
        k=sorted(m["origin_keys"])[0]
        m["gardens"][k].add((18.59,-31.1))
        with self.assertRaisesRegex(ValueError,"missing or ambiguous"):
            geometry(m,self.c)
    def test_new_unfrozen_source_qualifier_or_pair_distance_tolerance_fails(self):
        self.c["method"]["distance_bins_km"]=[1,3,5,20]
        with self.assertRaisesRegex(ValueError,"field or distance"):
            geometry(self.metadata(),self.c)
        self.c=copy.deepcopy(json.loads(CONTRACT.read_text(encoding="utf-8")))
        self.c["method"]["allowed_original_csv_fields"].append("fitness")
        with self.assertRaisesRegex(ValueError,"field or distance"):
            geometry(self.metadata(),self.c)
    def test_published_original_drainage_basin_not_synonym_for_LGM_effect(self):
        self.assertEqual(self.c["study"]["published_original_methods"]["physical_experiment"],
            "72 adult plant individuals collected from each population; each destination garden 18 from each of four origin populations; cross-basin same quartz habitat and within-basin other quartz habitat comparisons.")
        r=geometry(self.metadata(),self.c)
        self.assertTrue(r["species_contrasts_are_edaphic_vs_basin_not_direct_LGM_selection_evidence"])
        self.assertTrue(r["source_ID_pair_geometry_is_not_15_independent_historical_exposures"])
        self.assertTrue(r["original_fitness_or_genetic_outcomes_not_read_from_Hargreaves"])
        self.assertTrue(r["Bontrager_terminal_44_of_50_gate_unchanged"])
    def test_basin_label_can_be_explicit_but_habitat_missing_then_abstain(self):
        x=source_id_category("Sout-sample1",self.c)
        self.assertEqual(x["basin"],"Sout")
        self.assertIsNone(x["habitat"])
        self.assertIsNone(source_id_category("unknown000",self.c)["basin"])
if __name__=="__main__":
    unittest.main()
