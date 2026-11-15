"""Haldane presence-only mask tests: changing response magnitudes never changes QA."""
from __future__ import annotations
import csv,io,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed70_haldane_presence_mask import (CONTRACT,FIELDS,META,
     is_present_raw_token,missingness_summary)


def fixture(magnitude="100.2",blank_sd=False):
    columns=list(META)+list(FIELDS)+["extra"+str(i) for i in range(69-len(META)-len(FIELDS))]
    rows=[]
    for i,sp,ref in [(1,"Plantus alpha","st1"),(2,"Plantus alpha","st1"),
                     (3,"Plantus beta","st2"),(4,"Testus gamma","st3")]:
        d={k:"" for k in columns}
        d.update({"ref_ID":ref,"sp_pub":sp,"phylum":"Streptophyta" if i!=4 else "Chordata",
                  "taxa":"Plant" if i!=4 else "Mammal",
                  "genphen":"Genetic",
                  "sample1_latitude":"33","sample1_longitude":"-118",
                  "sample2_latitude":"34","sample2_longitude":"-119",
                  "trait_type":"size","design":"Synchronic",
                  "mean1":magnitude,"mean2":magnitude,
                  "sd1":"" if i==2 and blank_sd else magnitude,
                  "sd2":magnitude,
                  "n1":magnitude,"n2":magnitude,
                  "generations":magnitude})
        rows.append(d)
    out=io.StringIO();w=csv.DictWriter(out,fieldnames=columns,lineterminator="\n")
    w.writeheader()
    for d in rows:w.writerow(d)
    return out.getvalue().encode()


class BinaryHaldaneSourceTests(unittest.TestCase):
    def setUp(self):
        self.rules=json.loads(CONTRACT.read_text())

    def test_response_magnitude_permutation_does_not_affect_presence(self):
        a=missingness_summary(fixture("0.001"),self.rules,check_authority=False)
        b=missingness_summary(fixture("-9999999999"),self.rules,check_authority=False)
        self.assertEqual(a,b)
        self.assertEqual(a["original_geographic_plant_source_rows"],3)
        self.assertEqual(a["original_geographic_study_species_pairs"],2)
        self.assertEqual(a["rows_with_all_required_Haldane_source_fields_present"],3)
        self.assertEqual(a["distinct_species_with_any_complete_mask"],2)
        self.assertFalse(a["numeric_mean_SD_n_rate_or_fitness_values_decoded"])
        self.assertFalse(a["Haldane_or_Darwin_rates_calculated"])

    def test_missing_SD_only_changes_presence_mask(self):
        r=missingness_summary(fixture("12",blank_sd=True),self.rules,check_authority=False)
        self.assertEqual(r["rows_with_all_required_Haldane_source_fields_present"],2)
        self.assertEqual(r["distinct_original_study_species_geographic_pairs_with_any_complete_mask"],2)
        self.assertEqual(r["source_presence_only_counts_per_required_field"]["sd1"],2)

    def test_distinct_source_null_tokens(self):
        empty={x.lower() for x in self.rules["field_presence_only"]["absent_if_trimmed_lowercase_in"]}
        for token in ("", " NA ", "Nan","null","INF","-inf"):
            self.assertFalse(is_present_raw_token(token,empty))
        for token in ("0","-9999","1e-15","36.5"):
            self.assertTrue(is_present_raw_token(token,empty))

    def test_full_original_frozen_source_is_required_production(self):
        with self.assertRaisesRegex(ValueError,"Original 70-plant source-only geometry"):
            missingness_summary(fixture("1"),self.rules)

    def test_exposure_contract_does_not_admit_modeling(self):
        self.assertTrue(self.rules["scope_support_criteria"]["pass_does_not_authorise_empirical_response_test"])
        self.assertTrue(self.rules["field_presence_only"]["DO_NOT_decode_to_numbers_or_evaluate_numeric_rate_values"])
        self.assertEqual(list(FIELDS),self.rules["field_presence_only"]["columns"])


if __name__=="__main__":
    unittest.main()
