"""Synthetic known-truth tests for two-original source-only metadata overlap."""
import csv
import copy
import hashlib
import io
import json
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from bontrager_hargreaves_finegrain_overlap_v01 import (
    CONTRACT, metadata_census, compare, count_matched_garden_coords,
    norm_taxon, valid_coord, haversine,
)

class FinegrainOriginalOverlapTests(unittest.TestCase):
    def setUp(self):
        self.c=json.loads(CONTRACT.read_text())
        self.c["prior_reference_counts"]={
            "bontrager_qualified_studies":1,
            "bontrager_qualified_species":1,
            "bontrager_qualified_cells":1,
            "hargreaves_plant_study_labels":2,
            "exact_qualified_study_label_overlap":1,
        }
        self.ref={"catalog":[{"original_study_label":"Paper_JEcol_2016"}]}

    def source(self,cohort,rows):
        fields=list(self.c["strict_metadata_only_fields"])
        names=fields+["fitness","treatment_notes"]+[f"unused_{i}" for i in range(47-len(fields)-2)]
        self.assertEqual(len(names),47)
        out=io.StringIO(newline="")
        w=csv.writer(out,lineterminator="\n")
        w.writerow(names)
        for m in rows:
            w.writerow([m.get(x,"") for x in names])
        content=out.getvalue().encode("utf-8")
        if cohort=="hargreaves":
            self.c["original"]["hargreaves"]["expected_header_sha256"]=hashlib.sha256(content.split(b"\n",1)[0]+b"\n").hexdigest()
        return metadata_census(content,source=cohort,c=self.c)

    def m(self,study="Paper_JEcol_2016",taxon="Plant species",source="Paper_sourceA",site="Paper_siteA",
          fitness_type="surv",lat_site="35.0",long_site="139.0",fitness="CANNOT_BE_USED"):
        return dict(study=study,taxon=taxon,source=source,site=site,
            functional_group="annual",fitness_type=fitness_type,
            lat_site=lat_site,long_site=long_site,fitness=fitness)

    def test_exact_key_hierarchy_shared_and_unlinked_counts(self):
        left=self.source("bontrager",[
          self.m(),self.m(source="Paper_sourceB"),
          self.m(taxon="Plant alien",study="BontragerOnly_2015"),
        ])
        right=self.source("hargreaves",[
          self.m(fitness="some copyrighted response should never be used"),
          self.m(fitness="another response"),
          self.m(source="Paper_sourceZ"),
          self.m(study="HargreavesOnly_2017",taxon="Plant different"),
        ])
        result=compare(left,right,self.ref,self.c)
        self.assertEqual(result["qualified_source_study_overlap_reproduced"],1)
        self.assertEqual(result["all_Bontrager_Hargreaves_shared_original_study_labels"],1)
        self.assertEqual(result["Hargreaves_study_labels_not_in_any_Bontrager_original_plant_row"],1)
        self.assertEqual(result["nested_source_key_overlaps"]["study__taxon__source__site__fitness_type"]["exact_shared_keys"],1)
        self.assertEqual(result["nested_source_key_overlaps"]["study__taxon__site"]["exact_shared_keys"],1)
        self.assertEqual(result["shared_transplant_garden_coordinate_agreement"]["distance_le_1km"],1)
        self.assertEqual(result["response_columns_indexed_or_fitted"],[])

    def test_unqualified_full_original_study_reduces_unlinked_count(self):
        left=self.source("bontrager",[
          self.m(),self.m(study="ExtraOriginal_2018",taxon="Plant extra")
        ])
        right=self.source("hargreaves",[
          self.m(),self.m(study="ExtraOriginal_2018",taxon="Plant extra")
        ])
        result=compare(left,right,self.ref,self.c)
        self.assertEqual(result["Hargreaves_study_labels_only_unmatched_to_qualified_Bontrager_but_found_in_full_original"],1)
        self.assertEqual(result["Hargreaves_study_labels_not_in_any_Bontrager_original_plant_row"],0)

    def test_nonplant_excluded_before_label_join(self):
        left=self.source("bontrager",[self.m()])
        a=self.m(study="OnlyFungus_2016")
        a["functional_group"]="fungus"
        right=self.source("hargreaves",[self.m(),a])
        self.c["prior_reference_counts"]["hargreaves_plant_study_labels"]=1
        r=compare(left,right,self.ref,self.c)
        self.assertEqual(r["all_Hargreaves_plant_original_study_labels"],1)
        self.assertEqual(r["Hargreaves_study_labels_not_in_any_Bontrager_original_plant_row"],0)

    def test_fitness_text_never_affects_source_key(self):
        left=self.source("bontrager",[self.m(fitness="99.9")])
        right=self.source("hargreaves",[self.m(fitness="-INF INVALID FITNESS") ,
                                       self.m(study="HargreavesOnly_2017",taxon="Plant other")])
        out=compare(left,right,self.ref,self.c)
        self.assertEqual(out["nested_source_key_overlaps"]["study__taxon__source__site__fitness_type"]["exact_shared_keys"],1)
        self.assertFalse(out["fitness_effect_or_direction_opened"])

    def test_all_site_distance_bins_without_outputting_coords(self):
        x={("a","p","s"):{(139,35)},("b","p","s"):{(139,35)},("c","p","s"):{None},("d","p","s"):{(130,35)}}
        y={("a","p","s"):{(139,35.001)},("b","p","s"):{(139,35.05)},
           ("c","p","s"):{(139,35)},("d","p","s"):{(140,35)}}
        d=count_matched_garden_coords(x,y)
        self.assertEqual(d["shared_source_study_species_garden_keys"],4)
        self.assertEqual(d["distance_le_1km"],1)
        self.assertEqual(d["distance_gt_1km_le_10km"],1)
        self.assertEqual(d["distance_gt_10km"],1)
        self.assertEqual(d["missing_or_ambiguous"],1)

    def test_reject_fuzzy_taxon_name(self):
        self.assertEqual(norm_taxon("Plant_species"),"plant species")
        self.assertIsNone(norm_taxon("Plant sp."))
        self.assertIsNone(norm_taxon("Plant species subspecies"))

    def test_unevaluable_site_coordinates_not_imputed(self):
        self.assertIsNone(valid_coord("NaN","1.0"))
        self.assertIsNone(valid_coord("139","-95.0"))
        self.assertAlmostEqual(haversine((139,35),(139,35)),0.0)

    def test_source_reference_guard(self):
        left=self.source("bontrager",[self.m()])
        right=self.source("hargreaves",[self.m(),self.m(study="Other_2019",taxon="Plant other")])
        self.c["prior_reference_counts"]["exact_qualified_study_label_overlap"]=2
        with self.assertRaisesRegex(ValueError,"not reproduced"):
            compare(left,right,self.ref,self.c)

if __name__=="__main__":
    unittest.main()
