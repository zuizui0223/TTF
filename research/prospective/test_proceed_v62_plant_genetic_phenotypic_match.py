"""Pure synthetic source-label matches; never use genetic/phenotypic rate cells."""
import csv,io,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed_v62_plant_genetic_phenotypic_match import matched_summary,CONTRACT


def fake(rate="0.1"):
    rules=json.loads(CONTRACT.read_text())
    header=rules["source_metadata_fields_allowed"]+rules["forbidden_rate_columns"]
    header.extend("filler"+str(i) for i in range(69-len(header)))
    rows=[
        {"ref_ID":"s1","sp_pub":"Plantus alpha","phylum":"Streptophyta",
         "taxa":"Plant","genphen":"Genetic","trait_type":"size",
         "sample1_latitude":"30","sample1_longitude":"-99",
         "sample2_latitude":"31","sample2_longitude":"-98"},
        {"ref_ID":"s1","sp_pub":"Plantus alpha","phylum":"Streptophyta",
         "taxa":"Plant","genphen":"Phenotypic","trait_type":"size",
         "sample1_latitude":"30","sample1_longitude":"-99",
         "sample2_latitude":"31","sample2_longitude":"-98"},
        {"ref_ID":"s1","sp_pub":"Plantus alpha","phylum":"Streptophyta",
         "taxa":"Plant","genphen":"Phenotypic","trait_type":"size",
         "sample1_latitude":"32","sample1_longitude":"-99",
         "sample2_latitude":"34","sample2_longitude":"-98"},
        {"ref_ID":"s2","sp_pub":"Plantus beta","phylum":"Streptophyta",
         "taxa":"Plant","genphen":"Genetic","trait_type":"physio",
         "sample1_latitude":"30","sample1_longitude":"-99",
         "sample2_latitude":"31","sample2_longitude":"-98"},
        {"ref_ID":"s3","sp_pub":"Plantus beta","phylum":"Streptophyta",
         "taxa":"Plant","genphen":"Phenotypic","trait_type":"physio",
         "sample1_latitude":"30","sample1_longitude":"-99",
         "sample2_latitude":"31","sample2_longitude":"-98"},
        {"ref_ID":"s4","sp_pub":"Animalus gamma","phylum":"Chordata",
         "taxa":"Mammal","genphen":"Genetic","trait_type":"size",
         "sample1_latitude":"30","sample1_longitude":"-99",
         "sample2_latitude":"31","sample2_longitude":"-98"},
    ]
    o=io.StringIO()
    w=csv.DictWriter(o,fieldnames=header,lineterminator="\n")
    w.writeheader()
    for original in rows:
        row={key:"" for key in header};row.update(original)
        row["mean1"]=rate;row["mean2"]="9999" if rate=="0.1" else "-8888"
        row["sd1"]=rate
        w.writerow(row)
    return o.getvalue().encode()


class SourceMatchTests(unittest.TestCase):
    def test_invariant_to_phenotypic_rate_value(self):
        contract=json.loads(CONTRACT.read_text())
        a=matched_summary(fake("0.1"),contract)
        b=matched_summary(fake("-30000"),contract)
        self.assertEqual(a,b)
        self.assertEqual(a["strict_same_study_species_trait_ordered_geography_matched_keys"],1)
        self.assertEqual(a["strictly_matched_distinct_species"],1)
        self.assertEqual(a["strictly_matched_distinct_original_studies"],1)
        self.assertEqual(a["same_study_species_trait_distinct_species"],1)
        self.assertEqual(a["same_species_trait_distinct_species"],2)
        self.assertFalse(a["source_categories_not_plasticity_and_selection_decomposition"] is False)

    def test_forbidden_rate_access_fails(self):
        c=json.loads(CONTRACT.read_text())
        c["source_metadata_fields_allowed"].append("mean1")
        with self.assertRaisesRegex(ValueError,"Rate response"):
            matched_summary(fake(),c)

    def test_header_drift_stops(self):
        with self.assertRaises(ValueError):
            matched_summary(b"x,y\n1,2\n",json.loads(CONTRACT.read_text()))


if __name__=="__main__": unittest.main()
