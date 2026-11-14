"""Synthetic-only test: qualitative original source metadata never depends on rates."""
import csv,io,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed70_plant_rate_eligibility import (
    CONTRACT,response_blind_rate_metadata
)

def fixture(rate="0.02"):
    contract=json.loads(CONTRACT.read_text())
    allow=contract["metadata_fields_allowed"]
    outcomes=contract["explicitly_forbidden_response_fields"]
    columns=allow+outcomes
    columns.extend("other"+str(i) for i in range(69-len(columns)))
    rows=[
        {"es_ID":"1","ref_ID":"study1","sp_pub":"Plantus_alpha",
         "phylum":"Streptophyta","taxa":"Plant","genphen":"Genetic",
         "sample1_latitude":"30","sample1_longitude":"-99",
         "sample2_latitude":"31","sample2_longitude":"-98",
         "trait_type":"physio","data_type":"raw","design":"Synchronic",
         "years":"3","generations":"1","gtime":"3"},
        {"es_ID":"2","ref_ID":"study1","sp_pub":"Plantus_alpha",
         "phylum":"Streptophyta","taxa":"Plant","genphen":"Genetic",
         "sample1_latitude":"30","sample1_longitude":"-99",
         "sample2_latitude":"31","sample2_longitude":"-98",
         "trait_type":"physio","data_type":"raw","design":"Synchronic",
         "years":"3","generations":"1","gtime":"3"},
        {"es_ID":"3","ref_ID":"study2","sp_pub":"Plantus_beta",
         "phylum":"Streptophyta","taxa":"Plant","genphen":"Genetic",
         "sample1_latitude":"40","sample1_longitude":"-99",
         "sample2_latitude":"42","sample2_longitude":"-98",
         "trait_type":"size","data_type":"transformed","design":"Allochronic",
         "years":"4","generations":"","gtime":""},
        {"es_ID":"4","ref_ID":"study3","sp_pub":"Testus_gamma",
         "phylum":"Chordata","taxa":"Mammal","genphen":"Genetic",
         "sample1_latitude":"20","sample1_longitude":"-80",
         "sample2_latitude":"21","sample2_longitude":"-79",
         "trait_type":"size","data_type":"raw"},
        {"es_ID":"5","ref_ID":"study4","sp_pub":"Plantus_delta",
         "phylum":"Streptophyta","taxa":"Plant","genphen":"Phenotypic",
         "sample1_latitude":"20","sample1_longitude":"-80",
         "sample2_latitude":"21","sample2_longitude":"-79",
         "trait_type":"size","data_type":"raw"},
    ]
    o=io.StringIO()
    w=csv.DictWriter(o,fieldnames=columns,lineterminator="\n")
    w.writeheader()
    for x in rows:
        y={field:"" for field in columns};y.update(x)
        y["mean1"]=rate;y["mean2"]=str(10*float(rate))
        y["sd1"]="-99999" if rate.startswith("-") else "0.01"
        w.writerow(y)
    return o.getvalue().encode()

class TraitMetadataOnlyTests(unittest.TestCase):
    def test_rate_changes_do_not_change_trait_eligibility(self):
        c=json.loads(CONTRACT.read_text())
        a=response_blind_rate_metadata(fixture("0.02"),c,check_authoritative=False)
        b=response_blind_rate_metadata(fixture("-200000"),c,check_authoritative=False)
        self.assertEqual(a,b)
        self.assertEqual(a["source_total_original_records"],5)
        self.assertEqual(a["original_eligible_plant_effect_source_records"],3)
        self.assertEqual(a["trait_family_metadata"]["physio"]["distinct_study_species_ordered_geo_pairs"],1)
        self.assertEqual(a["trait_family_metadata"]["size"]["distinct_original_studies"],1)
        self.assertEqual(a["trait_family_metadata"]["size"]["generations_metadata_present"],0)
        self.assertFalse(a["phenotypic_change_or_genetic_effect_rate_values_read"])
        self.assertFalse(a["model_fit_or_history_effect_tested"])

    def test_forbidden_columns_fail_metadata_allowlist(self):
        c=json.loads(CONTRACT.read_text())
        c["metadata_fields_allowed"].append("mean1")
        with self.assertRaisesRegex(ValueError,"Forbidden rate"):
            response_blind_rate_metadata(fixture(),c,check_authoritative=False)

    def test_strict_original_source_dimensions_required(self):
        c=json.loads(CONTRACT.read_text())
        with self.assertRaisesRegex(ValueError,"Frozen original plant source"):
            response_blind_rate_metadata(fixture(),c)

if __name__=="__main__":unittest.main()
