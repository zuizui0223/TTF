"""Synthetic known-truth tests of fresh same-code, six-population source-site QA."""
import copy
import json
import sys
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from raabova2007_pankova2014_six_site_exploratory_v01 import (
    CONTRACT,strict_dms,fixed_published_sites,
    coded_garden_proxy_coordinates,exploratory_comparison
)

class SixPublishedCodeSiteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base=json.loads(CONTRACT.read_text(encoding="utf-8"))
        cls.err=json.loads((HERE/"raabova_2007_2008_wrong_primary_reference_erratum_v01.json").read_text())
    def setUp(self):
        self.c=copy.deepcopy(self.base)
    def metadata(self):
        pubs=fixed_published_sites(self.c)
        original_study=self.c["target"]["original_study"]
        taxa=self.c["target"]["strict_species"]
        keys={(original_study,taxa,"Raabova2007_"+k) for k in pubs}
        gardens={key:{pubs[key[2].rsplit("_",1)[1]]} for key in keys}
        return {"origin_keys":keys,"gardens":gardens}

    def test_2014_exact_six_site_order_and_DMS(self):
        sites=fixed_published_sites(self.c)
        self.assertEqual(list(sites),["K1","K2","K3","S1","S2","S3"])
        self.assertAlmostEqual(sites["K1"][0],14+7/60+29/3600,places=9)
        self.assertAlmostEqual(sites["K1"][1],49+56/60+1/3600,places=9)
        self.assertAlmostEqual(sites["S3"][1],50+31/60+46/3600,places=9)
        self.assertEqual(strict_dms("50:32:24N","lat"),50+32/60+24/3600)

    def test_six_exact_sites_are_descriptively_near_in_known_truth(self):
        r=exploratory_comparison(self.metadata(),self.c,self.err)
        self.assertEqual(r["exact_code_site_count"],6)
        self.assertEqual(r["descriptive_sites_within_1km"],6)
        self.assertTrue(all(x["distance_bucket_predeclared"]=="0-0.1km" for x in r["site_details"]))
        self.assertFalse(r["independent_population_GPS_validation_certified"])
        self.assertTrue(r["no_adaptation_plasticity_speciation_genetic_history_climate_response_opened"])

    def test_cross_study_population_codes_are_not_joined(self):
        m=self.metadata()
        key=next(x for x in m["origin_keys"] if x[2].endswith("_S2"))
        coords=m["gardens"].pop(key)
        m["gardens"][("Another_study",key[1],key[2])]=coords
        with self.assertRaisesRegex(ValueError,"No unique original"):
            coded_garden_proxy_coordinates(m,self.c)

    def test_two_original_code_tokens_in_source_id_rejected(self):
        m=self.metadata()
        origin=next(x for x in m["origin_keys"] if x[2].endswith("_K1"))
        m["origin_keys"].remove(origin)
        other=(origin[0],origin[1],"Raabova2007_K1_S2")
        m["origin_keys"].add(other)
        m["gardens"][other]=m["gardens"][origin]
        with self.assertRaisesRegex(ValueError,"Ambiguous or nonunique"):
            coded_garden_proxy_coordinates(m,self.c)

    def test_duplicate_locations_for_same_source_id_rejected(self):
        m=self.metadata()
        key=next(x for x in m["origin_keys"] if x[2].endswith("_S3"))
        m["gardens"][key].add((13.0,50.0))
        with self.assertRaisesRegex(ValueError,"No unique"):
            coded_garden_proxy_coordinates(m,self.c)

    def test_2007_2008_bad_reference_cannot_be_resurrected(self):
        err=copy.deepcopy(self.err)
        err["author_publication_referent"]["referenced_primary_DOI"]=self.c["target"]["original_paper_doi"]
        with self.assertRaisesRegex(ValueError,"invalid referent"):
            exploratory_comparison(self.metadata(),self.c,err)

    def test_no_postresult_tolerance_expansion(self):
        self.c["rules"]["precommitted_descriptive_distance_bins_km"]=[0.1,5,10]
        with self.assertRaisesRegex(ValueError,"rules drifted"):
            exploratory_comparison(self.metadata(),self.c,self.err)

    def test_wrong_external_DOI_or_population_table_fails(self):
        self.c["external_primary_author_source"]["doi"]="10.1371/journal.pone.invalid"
        with self.assertRaisesRegex(ValueError,"PLOS"):
            fixed_published_sites(self.c)
        self.c=copy.deepcopy(self.base)
        self.c["external_primary_author_source"]["published_site_rows"][0]["code"]="K2"
        with self.assertRaisesRegex(ValueError,"site rows"):
            fixed_published_sites(self.c)

    def test_invalid_DMS_and_region_letter_fails(self):
        with self.assertRaises(ValueError):
            strict_dms("50:61:24N","lat")
        with self.assertRaises(ValueError):
            strict_dms("14:05:16N","lon")
        self.c["external_primary_author_source"]["published_site_rows"][0]["region"]="S"
        with self.assertRaisesRegex(ValueError,"region"):
            fixed_published_sites(self.c)

    def test_misplaced_same_code_site_yields_descriptive_distance_not_automated_rescue(self):
        m=self.metadata()
        key=next(x for x in m["origin_keys"] if x[2].endswith("_S1"))
        m["gardens"][key]={(15.0,49.0)}
        r=exploratory_comparison(m,self.c,self.err)
        self.assertEqual(r["descriptive_sites_within_1km"],5)
        self.assertEqual(r["descriptive_sites_over_1km"],1)
        self.assertTrue(r["old_invalid_Sudova_2010_vs_2007_receipt_remains_invalid"])

if __name__=="__main__":unittest.main()
