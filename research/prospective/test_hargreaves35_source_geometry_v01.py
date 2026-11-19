"""Synthetic-only tests for source geometry of exactly 35 original-study labels."""
import copy
import csv
import hashlib
import io
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hargreaves35_source_geometry_v01 import CONTRACT, audit_hargreaves35

HERE=Path(__file__).resolve().parent
SOURCE_CONTRACT=HERE/"hargreaves_zenodo_study_overlap_contract_v01.json"

class Hargreaves35GeometryTests(unittest.TestCase):
    def setUp(self):
        self.contract=copy.deepcopy(json.loads(CONTRACT.read_text()))
        self.labels=self.contract["reference"]["exact_unlinked_study_labels"]
        self.alpha,self.beta=self.labels[0],self.labels[1]
        self.receipt={
            "Hargreaves_unlinked_full_study_labels_for_followup":list(self.labels),
            "Hargreaves_study_labels_not_in_any_Bontrager_original_plant_row":35,
        }
        self.header=json.loads(SOURCE_CONTRACT.read_text())["source"]["original_csv_header"]

    def row(self, study=None, source="origin_a", garden="garden_a", taxon="Plant species",
            experiment="exp_1", group="annual", ftype="surv", lon="139", lat="35",
            fitness="EXTREME_FAKE_NONNUMERIC_RESPONSE"):
        return {
            "study":study or self.alpha,"source":source,"site":garden,"taxon":taxon,
            "experiment":experiment,"functional_group":group,"fitness_type":ftype,
            "long_site":lon,"lat_site":lat,"year_begin":"2009",
            "fitness":fitness,
        }

    def original(self,rows):
        s=io.StringIO(newline="")
        writer=csv.writer(s,lineterminator="\n")
        writer.writerow(self.header)
        for row in rows:
            writer.writerow([row.get(col,"") for col in self.header])
        raw=s.getvalue().encode("cp1252")
        source=self.contract["source"]
        source["exact_original_byte_count"]=len(raw)
        source["expected_original_md5"]=hashlib.md5(raw).hexdigest()
        source["expected_original_sha256"]=hashlib.sha256(raw).hexdigest()
        source["expected_header_sha256"]=hashlib.sha256(raw.split(b"\n",1)[0]+b"\n").hexdigest()
        source["expected_original_physical_rows"]=len(rows)
        return raw

    def test_exact_same_garden_and_experiment_two_source_labels_are_candidate_only(self):
        raw=self.original([self.row(source="a"),self.row(source="b"),
                           self.row(study=self.beta,source="a"),
                           self.row(study="SourceNotInFrozen35",source="a")])
        out=audit_hargreaves35(raw,self.contract,self.receipt)
        self.assertEqual(out["summary"]["original_study_labels_previously_unlinked"],35)
        self.assertEqual(out["summary"]["study_labels_with_same_garden_multi_source_opportunity"],1)
        self.assertEqual(out["summary"]["sum_candidate_same_garden_cells_not_independent_experiments"],1)
        self.assertEqual(out["summary"]["study_labels_with_no_original_experiment_ID"],33)
        self.assertEqual(out["summary"]["verified_source_population_georeferenced_studies_in_Hargreaves_original"],0)
        self.assertFalse(out["fitness_values_or_effects_indexed"])

    def test_split_experiments_or_metric_types_not_pooled(self):
        raw=self.original([self.row(source="a",experiment="experimentA",ftype="surv"),
                           self.row(source="b",experiment="experimentB",ftype="surv"),
                           self.row(source="c",experiment="experimentA",ftype="germ")])
        out=audit_hargreaves35(raw,self.contract,self.receipt)
        self.assertEqual(out["summary"]["sum_candidate_same_garden_cells_not_independent_experiments"],0)

    def test_missing_experiment_is_unknown_not_certified(self):
        raw=self.original([self.row(source="a",experiment=""),
                           self.row(source="b",experiment="")])
        out=audit_hargreaves35(raw,self.contract,self.receipt)
        self.assertEqual(out["summary"]["sum_candidate_cells_missing_experiment_ID"],1)
        self.assertEqual(out["summary"]["verified_independent_original_experiments"],0)

    def test_garden_geo_unknown_or_conflicted_flagged(self):
        raw=self.original([self.row(source="a",lat="",lon=""),
                           self.row(source="b",lat="35",lon="139")])
        out=audit_hargreaves35(raw,self.contract,self.receipt)
        study=next(x for x in out["studies"] if x["original_study_label"]==self.alpha)
        self.assertEqual(study["candidate_same_garden_cells_missing_or_ambiguous_garden_coordinate"],1)
        self.assertEqual(study["candidate_same_garden_cells_with_complete_garden_coordinate"],0)

    def test_unrelated_and_nonplant_records_are_not_in_target_geometry(self):
        raw=self.original([self.row(source="a"),self.row(source="b",group="arthropod"),
                           self.row(study="OtherPaper",source="c")])
        out=audit_hargreaves35(raw,self.contract,self.receipt)
        self.assertEqual(out["summary"]["plant_rows_in_35_unlinked_source_labels"],1)
        self.assertEqual(out["summary"]["sum_candidate_same_garden_cells_not_independent_experiments"],0)

    def test_fitness_values_do_not_change_geometry(self):
        rows=[self.row(source="a",fitness="-500"),self.row(source="b",fitness="+100000")]
        raw=self.original(rows)
        a=audit_hargreaves35(raw,self.contract,self.receipt)
        rows=[self.row(source="a",fitness="NONSENSE"),self.row(source="b",fitness="NO_RAW_RESPONSE")]
        raw=self.original(rows)
        b=audit_hargreaves35(raw,self.contract,self.receipt)
        self.assertEqual(a["summary"],b["summary"])
        self.assertEqual(a["studies"],b["studies"])

    def test_changed_original_file_or_reference_list_fails_before_census(self):
        raw=self.original([self.row()])
        with self.assertRaisesRegex(ValueError,"bytes/digests"):
            audit_hargreaves35(raw+b"BAD",self.contract,self.receipt)
        self.receipt["Hargreaves_study_labels_not_in_any_Bontrager_original_plant_row"]=34
        with self.assertRaisesRegex(ValueError,"not reproduced"):
            audit_hargreaves35(raw,self.contract,self.receipt)

if __name__=="__main__":
    unittest.main()
