"""Tests for source-only 95-grid blocked LGM climate innovation calculation."""
import copy,json,math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import albano2025_source_only_MPIO_GCM_historical_innovation_v01 as a
GRID=json.loads((Path(__file__).resolve().parent/
    "source_registry/albano2025_unlabeled_worldclim_2_5m_cells_v01.json").read_text("utf-8"))
def fixture():
    modern=[]
    lgm=[]
    for obj in GRID["cells"]:
        row,col=obj["row"],obj["col"]
        current1=100+row%23+col%7
        current12=500+(row*3)%67+(col*11)%43
        old1=current1-70+(row%19)
        old12=current12+105+(col%13)*13-(row%7)*5
        modern.append({"row":row,"col":col,"BIO1_raw":current1,"BIO12_raw":current12})
        lgm.append({"row":row,"col":col,"BIO1_raw":old1,"BIO12_raw":old12})
    return ({"schema":"ttf_Albano2025_WorldClim14_modern_1960_1990_actual_cells_v01",
             "baseline_period":"1960-1990","phenotype_fitness_genetic_responses_used":False,
             "cells":modern},
            {"schema":"ttf_Albano2025_MPI_ESM_P_source_only_2_5m_LGM_cells_v01",
             "unlabeled_grid_cell_count":95,"fitness_genetics_adaptation_responses_opened":False,
             "cell_profiles":lgm})
class ClimateInnovationTests(unittest.TestCase):
    def setUp(self):self.now,self.past=fixture()
    def test_preparation_has_95_unique_cells(self):
        rows=a.prepare(self.now,self.past)
        self.assertEqual(len(rows),95)
        self.assertEqual(set(r["fold"] for r in rows),set(range(10)))
    def test_continent_counts_original_geometry(self):
        rows=a.prepare(self.now,self.past)
        self.assertEqual((sum(x["origin"]=="NAM" for x in rows),
                          sum(x["origin"]=="EUR" for x in rows)),(47,48))
    def test_each_frozen_fold_nonempty_and_rank_sufficient(self):
        rows=a.prepare(self.now,self.past)
        for f in range(10):
            a.rank_and_prediction(rows,f,"delta_bio1")
    def test_source_only_climate_metric_produced(self):
        d=a.evaluate(self.now,self.past)
        self.assertEqual(set(d["climate_metrics"]),{"BIO1","BIO12"})
        self.assertFalse(d["fitness_genotype_trait_response_tested"])
        self.assertFalse(d["independent_LGM_effect_on_fitness_tested"])
        self.assertTrue(d["primary_CCSM4_raster_not_yet_tested"])
    def test_independent_bio1_design_not_phenotype(self):
        r=a.prepare(self.now,self.past)
        vals,pred=a.rank_and_prediction(r,0,"delta_bio1")
        self.assertEqual(len(vals),len(pred))
        self.assertTrue(all(math.isfinite(float(x)) for x in pred))
    def test_missing_matched_cell_fails(self):
        modern=copy.deepcopy(self.now);modern["cells"].pop()
        with self.assertRaisesRegex(ValueError,"join incomplete"):
            a.prepare(modern,self.past)
    def test_biological_response_access_rejected(self):
        past=copy.deepcopy(self.past)
        past["fitness_genetics_adaptation_responses_opened"]=True
        with self.assertRaisesRegex(ValueError,"biological response"):
            a.prepare(self.now,past)
    def test_original_climate_period_frozen(self):
        modern=copy.deepcopy(self.now)
        modern["baseline_period"]="1970-2000"
        with self.assertRaisesRegex(ValueError,"period/cell"):
            a.prepare(modern,self.past)
    def test_historical_response_not_observed_fitness(self):
        d=a.evaluate(self.now,self.past)
        self.assertIn("climate",str(d["status"]).lower())
        self.assertFalse(d["independent_LGM_effect_on_fitness_tested"])
if __name__=="__main__":unittest.main()
