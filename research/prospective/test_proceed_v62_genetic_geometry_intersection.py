"""Synthetic-only cross-qualification tests, never ingest published rate values."""
from __future__ import annotations
import csv,io,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed_v62_genetic_geometry_intersection import intersection_metadata_summary,CONTRACT
from proceed_v62_metadata_row_census import CONTRACT as BASE


def make_fake(rate="0.01"):
    allowed=json.loads(BASE.read_text())["source_metadata_columns_authorized"]
    cols=list(allowed)+["mean1","mean2","sd1","sd2","n1","n2","transf_data","data_scale","out_in"]
    cols.extend("unused"+str(i) for i in range(69-len(cols)))
    rows=[
      {"es_ID":"1","ref_ID":"r1","sp_pub":"Testus alpha","genphen":"Genetic",
       "phylum":"Chordata","sample1_latitude":"30","sample1_longitude":"-99",
       "sample2_latitude":"31","sample2_longitude":"-98","design":"Allochronic"},
      {"es_ID":"2","ref_ID":"r1","sp_pub":"Testus alpha","genphen":"Phenotypic",
       "phylum":"Chordata","sample1_latitude":"32","sample1_longitude":"-97",
       "sample2_latitude":"34","sample2_longitude":"-96","design":"Synchronic"},
      {"es_ID":"3","ref_ID":"r2","sp_pub":"Plantus beta","genphen":"Genetic",
       "phylum":"Streptophyta","sample1_latitude":"33","sample1_longitude":"-97",
       "sample2_latitude":"34","sample2_longitude":"-96"},
      {"es_ID":"4","ref_ID":"r3","sp_pub":"Anotherus gamma","genphen":"Genetic",
       "phylum":"Arthropoda","sample1_latitude":"33","sample1_longitude":"-96",
       "sample2_latitude":"","sample2_longitude":""},
    ]
    out=io.StringIO()
    writer=csv.DictWriter(out,fieldnames=cols,lineterminator="\n")
    writer.writeheader()
    for z in rows:
        z["mean1"]=rate;z["mean2"]="999" if rate=="0.01" else "-999"
        writer.writerow(z)
    return out.getvalue().encode()


class StratificationTests(unittest.TestCase):
    def test_phenotype_change_does_not_affect_metadata_eligibility(self):
        contract=json.loads(CONTRACT.read_text())
        a=intersection_metadata_summary(make_fake(),contract)
        b=intersection_metadata_summary(make_fake(rate="-100000"),contract)
        self.assertEqual(a,b)
        self.assertEqual(a["source_rows"],4)
        self.assertEqual(a["animal_phylum_geolocated_genetic"]["source_effect_rows"],1)
        self.assertEqual(a["animal_phylum_geolocated_phenotypic"]["source_effect_rows"],1)
        self.assertFalse(a["no_rate_or_genetic_effect_values_used"] is False)
        self.assertTrue(a["animal_phylum_is_not_verified_terrestrial_habitat"])

    def test_plant_genetic_record_does_not_count_as_animal(self):
        out=intersection_metadata_summary(make_fake(),json.loads(CONTRACT.read_text()))
        self.assertEqual(out["all_geolocated_genetic_records"]["source_effect_rows"],2)
        self.assertEqual(out["animal_phylum_geolocated_genetic"]["source_effect_rows"],1)
        self.assertEqual(out["genphen_by_geostatus"]["Genetic|first_only"],1)

    def test_requires_exact_expected_response_blind_header(self):
        with self.assertRaises(ValueError):
            intersection_metadata_summary(b"bad,header\n1,2\n",json.loads(CONTRACT.read_text()))


if __name__=="__main__":unittest.main()
