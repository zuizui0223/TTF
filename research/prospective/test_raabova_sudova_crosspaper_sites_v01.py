"""Original-paper published DMS cross-paper QA, synthetic geography only."""
import copy,json,sys,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from raabova_sudova_crosspaper_sites_v01 import CONTRACT, dms, primary_sites, colocate

class RaabovaSudovaCrosspaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=json.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_published_DMS_original_site_geometry(self):
        sites=primary_sites(self.c)
        self.assertEqual([x[0] for x in sites],[3,1,5,6])
        self.assertAlmostEqual(sites[0][1][0],14+15/60+32/3600,places=7)
        self.assertAlmostEqual(sites[0][1][1],50+31/60+45/3600,places=7)

    def test_invalid_DMS_and_lat_lon_swap_fail(self):
        for d,axis in [("50:61:20N","lat"),("50:30:61N","lat"),
                       ("50:31:45E","lat"),("14:15:32N","lon"),
                       ("181:00:00E","lon"),("0,0","lon")]:
            with self.assertRaises(ValueError):
                dms(d,axis)

    def test_frozen_sites_unique_real_same_coordinates_match(self):
        primary=primary_sites(self.c)
        candidates=sorted([xy for _,xy in primary]+[(15.,49.),(16.,49.)])
        res=colocate(primary,candidates,self.c)
        self.assertEqual(res["published_sites_spatially_corroborated_uniquely_within_1km"],4)
        self.assertEqual(res["published_sites_ambiguous_or_unmatched"],0)
        self.assertFalse(res["primary_paper_exact_population_ID_to_Hargreaves_source_ID_crosswalk_certified"])

    def test_large_location_error_not_mapped(self):
        original=primary_sites(self.c)
        candidates=[xy for _,xy in original[1:]]+[(15.,49.),(16.,49.)]
        res=colocate(original,candidates,self.c)
        self.assertEqual(res["published_sites_spatially_corroborated_uniquely_within_1km"],3)
        self.assertEqual(res["published_sites_ambiguous_or_unmatched"],1)

    def test_close_second_choice_is_ambiguous_even_if_best_matches(self):
        primary=primary_sites(self.c)
        x,y=primary[0][1]
        nearby=(x+0.0005,y)
        others=[xy for _,xy in primary[1:]]
        c=colocate(primary,[primary[0][1],nearby]+others,self.c)
        self.assertFalse(c["external_original_site_match_details"][0]["unique_nearest_population_1km_corroboration"])
        self.assertTrue(c["external_original_site_match_details"][0]["nearest_second_choice_ambiguous"])

    def test_duplicate_candidates_rejected(self):
        primary=primary_sites(self.c)
        same=[primary[0][1],primary[0][1]]
        with self.assertRaisesRegex(ValueError,"Nonunique"):
            colocate(primary,same,self.c)

    def test_tuning_distance_radius_requires_new_contract(self):
        a=copy.deepcopy(self.c)
        a["prespecified_distance_rule"]["maximum_distance_for_cross_paper_corroboration_km"]=5.
        with self.assertRaisesRegex(ValueError,"tuning"):
            colocate(primary_sites(a),[xy for _,xy in primary_sites(a)],a)

if __name__=="__main__":
    unittest.main()
