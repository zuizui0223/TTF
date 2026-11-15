import copy,json,sys,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from raabova2007_castro2011_named_sites_exploratory_v01 import (
    CONTRACT,match_by_verified_primary_population_codes,dms
)
ERR=HERE/"raabova_2007_2008_wrong_primary_reference_erratum_v01.json"

class CorrectStudyPublishedSiteTests(unittest.TestCase):
    def setUp(self):
        self.c=copy.deepcopy(json.loads(CONTRACT.read_text(encoding="utf-8")))
        self.err=copy.deepcopy(json.loads(ERR.read_text(encoding="utf-8")))
        self.study=self.c["primary_original"]["Hargreaves_study_label"]
        self.taxon=self.c["primary_original"]["taxon"]
    def source(self):
        codes=self.c["primary_original"]["published_2007_population_codes"]
        origin={(self.study,self.taxon,self.study+"_"+code) for code in codes}
        gardens={}
        for code in codes:
            gardens[(self.study,self.taxon,self.study+"_"+code)]={(13.,50.)}
        for s in self.c["corroborating_external_publication"]["Table1_diploid_populations"]:
            code=s["primary_2007_population_code"]
            gardens[(self.study,self.taxon,self.study+"_"+code)]={
                (dms(s["longitude_DMS"],"lon"),dms(s["latitude_DMS"],"lat"))}
        return {"origin_keys":origin,"gardens":gardens}
    def test_known_truth_unique_same_named_sites_match(self):
        a=match_by_verified_primary_population_codes(self.source(),self.c,self.err)
        self.assertEqual(a["published_named_diploid_localities_compared"],3)
        self.assertEqual(a["within_1km_descriptive_site_count"],3)
        self.assertTrue(a["old_invalid_1_of_4_external_test_remains_invalid"])
        self.assertEqual(a["independently_verified_raw_field_GPS_coordinates"],0)
        self.assertTrue(a["no_fitness_genetic_LGM_adaptation_plasticity_speciation_outcomes_fitted"])
    def test_wrong_paper_referent_never_accepted(self):
        self.err["author_publication_referent"]["referenced_primary_DOI"]=self.c["primary_original"]["paper_doi"]
        with self.assertRaisesRegex(ValueError,"literature identity"):
            match_by_verified_primary_population_codes(self.source(),self.c,self.err)
    def test_incomplete_study_population_IDs_abstain(self):
        a=self.source()
        a["origin_keys"]={key for key in a["origin_keys"] if not key[2].endswith("_S1")}
        with self.assertRaisesRegex(ValueError,"six source"):
            match_by_verified_primary_population_codes(a,self.c,self.err)
    def test_wrong_site_name_is_blocked_before_geometry(self):
        a=self.c["corroborating_external_publication"]["Table1_diploid_populations"]
        a[0]["site_name"]="Another population"
        with self.assertRaisesRegex(ValueError,"site-name evidence"):
            match_by_verified_primary_population_codes(self.source(),self.c,self.err)
    def test_false_nearest_location_cannot_rescue_labelled_codes(self):
        a=self.source()
        s1=(self.study,self.taxon,self.study+"_S1")
        # Intentionally make the S1 code far away; a different code may
        # happen to be geographically close but literal mapping is fixed.
        a["gardens"][s1]={(15.,49.)}
        r=match_by_verified_primary_population_codes(a,self.c,self.err)
        self.assertEqual(r["within_1km_descriptive_site_count"],2)
        self.assertFalse(r["exploratory_published_name_matched_sites"][0]["descriptive_distance_le_original_1km_cutoff"])
    def test_no_automatic_threshold_change(self):
        self.c["matching_rule"]["descriptive_cutoff_km"]=5
        with self.assertRaisesRegex(ValueError,"distance/matching"):
            match_by_verified_primary_population_codes(self.source(),self.c,self.err)
    def test_site_coordinate_missing_or_ambiguous_abstains(self):
        a=self.source()
        k=(self.study,self.taxon,self.study+"_S2")
        a["gardens"][k]={(14.23,50.53),None}
        with self.assertRaisesRegex(ValueError,"Ambiguous/missing"):
            match_by_verified_primary_population_codes(a,self.c,self.err)

if __name__=="__main__":unittest.main()
