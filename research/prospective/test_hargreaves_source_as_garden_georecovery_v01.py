"""Strict source-as-garden response-blind known-truth qualification tests."""
import copy
import csv
import io
import json
import sys
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from hargreaves_source_as_garden_georecovery_v01 import (
    CONTRACT, source_metadata, evaluation, one_point, percent_at
)

class ExactReciprocalSiteTests(unittest.TestCase):
    def setUp(self):
        self.rule=json.loads(CONTRACT.read_text())
        self.study="FirstPaper_Ecology_2012"
        self.study2="SecondPaper_Ecology_2012"
        self.targets=[self.study,self.study2]+[f"Empty{i}" for i in range(33)]
        self.receipt={
            "Hargreaves_unlinked_full_study_labels_for_followup":self.targets
        }
        self.rule["study_label_partition"]["study_count"]=35

    def records(self, cohort, rows):
        allow=self.rule["source_only_allowlist_"+cohort]
        fields=allow+["fitness","response_text"]+[f"extra{i}" for i in range(47-len(allow)-2)]
        output=io.StringIO(newline="")
        writer=csv.writer(output,lineterminator="\n")
        writer.writerow(fields)
        for vals in rows:
            writer.writerow([vals.get(k,"") for k in fields])
        data=output.getvalue().encode("utf-8")
        return source_metadata(data,cohort,self.rule)

    def r(self,study=None,source="p1",site="p1",lon="139",lat="35",
          source_lon="139",source_lat="35",taxon="Plant species",
          experiment="exp1",metric="surv",fitness="DO NOT READ"):
        return dict(study=study or self.study,functional_group="annual",
                    taxon=taxon,source=source,site=site,
                    lat_site=lat,long_site=lon,
                    lat_source=source_lat,long_source=source_lon,
                    experiment=experiment,fitness_type=metric,fitness=fitness)

    def test_exact_same_study_population_id_recovers_and_benchmarks(self):
        source=[self.r(source="p1",site="p1",lon="139",source_lon="139"),
                self.r(source="p2",site="p2",lon="140",source_lon="140")]
        # 2 independent origin IDs both present as experimental destinations.
        h=self.records("hargreaves",source+[
            self.r(source="p1",site="GardenOther",lon="141"),
            self.r(source="p2",site="GardenOther",lon="141")])
        b=self.records("bontrager",source)
        result=evaluation(h,b,self.rule,self.receipt)
        a=result["summary"]
        self.assertEqual(a["calibration_true_source_comparisons"],2)
        self.assertEqual(a["calibration_distance_bins"]["km_le_1"],2)
        self.assertEqual(a["unlinked_study_labels_with_any_same_garden_two_georecovered_origins"],1)
        self.assertEqual(a["unlinked_design_cells_with_two_georecovered_origins"],1)
        self.assertTrue(result["no_lgm_climate_or_genetic_fitness_effect_opened"])

    def test_different_study_or_species_never_cross_matches(self):
        h=self.records("hargreaves",[
            self.r(source="p1",site="G1"),
            self.r(study=self.study2,source="p2",site="p1"),
            self.r(source="p3",site="p3",taxon="Second species"),
        ])
        b=self.records("bontrager",[self.r()])
        result=evaluation(h,b,self.rule,self.receipt)
        self.assertEqual(result["summary"]["original_Hargreaves_source_ID_keys_with_unique_garden_coordinate"],1)
        self.assertEqual(result["summary"]["calibration_true_source_comparisons"],0)

    def test_conflicting_same_garden_coordinates_abstain(self):
        h=self.records("hargreaves",[
            self.r(source="p1",site="p1",lon="139"),
            self.r(source="p1",site="p1",lon="140"),
            self.r(source="p2",site="p2",lon="141")
        ])
        b=self.records("bontrager",[self.r(source="p1",site="p1")])
        result=evaluation(h,b,self.rule,self.receipt)
        self.assertEqual(result["summary"]["calibration_true_source_comparisons"],0)
        self.assertEqual(result["summary"]["unlinked_distinct_source_IDs_eligible_for_garden_derived_proxy"],1)

    def test_garden_with_one_valid_and_one_missing_point_abstains(self):
        h=self.records("hargreaves",[
            self.r(source="p1",site="p1",lon="139"),
            self.r(source="p1",site="p1",lon="")
        ])
        b=self.records("bontrager",[self.r(source="p1",site="p1")])
        result=evaluation(h,b,self.rule,self.receipt)
        self.assertEqual(result["summary"]["unlinked_distinct_source_IDs_eligible_for_garden_derived_proxy"],0)

    def test_fake_fitness_never_changes_georecovery(self):
        rows=[self.r(source="p1",site="p1"),self.r(source="p2",site="p2",lon="141")]
        b=self.records("bontrager",rows)
        a=evaluation(self.records("hargreaves",rows),b,self.rule,self.receipt)["summary"]
        rows[0]["fitness"]="-INFINITY_EXTREME"
        rows[1]["fitness"]="not a number"
        out=evaluation(self.records("hargreaves",rows),b,self.rule,self.receipt)["summary"]
        self.assertEqual(a,out)

    def test_missing_site_never_uses_study_centroid(self):
        h=self.records("hargreaves",[self.r(source="p1",site="gardenOnly"),
                                    self.r(source="p2",site="gardenOnly")])
        b=self.records("bontrager",[self.r(source="p1",site="p1")])
        a=evaluation(h,b,self.rule,self.receipt)["summary"]
        self.assertEqual(a["original_Hargreaves_source_ID_keys_existing_as_same_species_gardens"],0)
        self.assertEqual(a["unlinked_design_cells_with_two_georecovered_origins"],0)

    def test_calibration_distance_bins(self):
        h=self.records("hargreaves",[
            self.r(source="p1",site="p1",lon="139"),
            self.r(source="p2",site="p2",lon="140"),
        ])
        b=self.records("bontrager",[
            self.r(source="p1",site="p1",source_lon="139"),
            self.r(source="p2",site="p2",source_lon="139"),
        ])
        d=evaluation(h,b,self.rule,self.receipt)["summary"]["calibration_distance_bins"]
        self.assertEqual(d["km_le_1"],1)
        self.assertEqual(d["km_gt_10"],1)

    def test_no_overspecified_response_field_in_allowlist(self):
        self.rule["source_only_allowlist_hargreaves"].append("fitness")
        with self.assertRaisesRegex(ValueError,"allowlist drifted"):
            self.records("hargreaves",[self.r()])

    def test_singleton_and_percentile_calculations(self):
        self.assertIsNone(one_point({(1,2),None}))
        self.assertIsNone(one_point({(1,2),(1,3)}))
        self.assertEqual(one_point({(1,2)}),(1,2))
        self.assertEqual(percent_at([0,2,4],0.5),2)

if __name__=="__main__":
    unittest.main()
