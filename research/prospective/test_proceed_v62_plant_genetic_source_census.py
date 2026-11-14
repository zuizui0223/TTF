"""Synthetic-only source test for plant genetic-category independence and no rate use."""
import csv,io,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed_v62_plant_genetic_source_census import (CONTRACT,source_only_plant_counts,split_connected_components)
from proceed_v62_genetic_geometry_intersection import ROW_CONTRACT


def synthetic(m1="42"):
    allow=json.loads(ROW_CONTRACT.read_text())["source_metadata_columns_authorized"]
    hdr=allow+["mean1","mean2","sd1","sd2","n1","n2","transf_data","data_scale","out_in"]
    hdr.extend("other"+str(k) for k in range(69-len(hdr)))
    records=[
       {"es_ID":"1","ref_ID":"paper1","sp_pub":"Plantus alpha","phylum":"Streptophyta",
        "taxa":"Plant","genphen":"Genetic","sample1_latitude":"20","sample1_longitude":"10",
        "sample2_latitude":"21","sample2_longitude":"11","design":"Synchronic","mean1":m1},
       {"es_ID":"2","ref_ID":"paper1","sp_pub":"Plantus beta","phylum":"Streptophyta",
        "taxa":"Plant","genphen":"Genetic","sample1_latitude":"21","sample1_longitude":"11",
        "sample2_latitude":"22","sample2_longitude":"12","design":"Synchronic","mean1":m1},
       {"es_ID":"3","ref_ID":"paper2","sp_pub":"Plantus alpha","phylum":"Streptophyta",
        "taxa":"Plant","genphen":"Genetic","sample1_latitude":"19","sample1_longitude":"10",
        "sample2_latitude":"20","sample2_longitude":"11","design":"Allochronic","mean1":m1},
       {"es_ID":"4","ref_ID":"paper3","sp_pub":"Testus gamma","phylum":"Chordata",
        "taxa":"Mammal","genphen":"Genetic","sample1_latitude":"19","sample1_longitude":"10",
        "sample2_latitude":"20","sample2_longitude":"11","mean1":m1},
       {"es_ID":"5","ref_ID":"paper4","sp_pub":"Plantus delta","phylum":"Streptophyta",
        "taxa":"Plant","genphen":"Phenotypic","sample1_latitude":"19","sample1_longitude":"10",
        "sample2_latitude":"20","sample2_longitude":"11","mean1":m1},
    ]
    o=io.StringIO();w=csv.DictWriter(o,fieldnames=hdr,lineterminator="\n")
    w.writeheader()
    for r in records:w.writerow(r)
    return o.getvalue().encode()


class PlantScreenTests(unittest.TestCase):
    def test_response_values_do_not_alter_eligibility(self):
        c=json.loads(CONTRACT.read_text())
        a=source_only_plant_counts(synthetic("10"),c)
        b=source_only_plant_counts(synthetic("9099099999"),c)
        self.assertEqual(a,b)
        self.assertEqual(a["distinct_plant_binomials"],2)
        self.assertEqual(a["distinct_original_source_studies"],2)
        self.assertEqual(a["plant_genetic_georeferenced_source_rows"],3)
        self.assertFalse(a["actual_effect_size_prediction_performed"])

    def test_component_split_blocks_study_leakage(self):
        graph={"r1":{"A alpha","B beta"},"r2":{"B beta","C gamma"},
               "r3":{"D delta"},"r4":{"E epsilon","F phi"}}
        r=split_connected_components(graph)
        self.assertTrue(r["roles_study_species_disjoint"])
        self.assertEqual(r["components"],3)
        self.assertEqual(sum(r["species_per_role"].values()),6)

    def test_non_csv_header_abstains(self):
        with self.assertRaises(ValueError):
            source_only_plant_counts(b"A,B,C\n1,2,3\n",json.loads(CONTRACT.read_text()))


if __name__=="__main__":unittest.main()
