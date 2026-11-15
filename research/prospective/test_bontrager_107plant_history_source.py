"""Synthetic-only Bontrager matched 21ka/0BP source qualification, no fitness data."""
from __future__ import annotations
import csv,io,json,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from bontrager_107plant_history_source import (
    extract_source_cells,roles_for_graph,matched_source_only,
    h_source_summary,haversine_km,sha_float64
)
BASE=json.loads((Path(__file__).resolve().parent/"bontrager_transplant_geo_contract_v01.json").read_text())
LOCK=json.loads((Path(__file__).resolve().parent/"bontrager_107plant_history_contrast_contract_v01.json").read_text())


def fake_source(fitness="0.1"):
    header=["study","functional_group","source","site","fitness_type","taxon",
            "lat_source","long_source","lat_site","long_site","fitness"]
    header.extend("unused"+str(i) for i in range(47-len(header)))
    rows=[]
    for name,study,group in [("Plantus alpha","s1","annual"),
                             ("Insectus beta","s2","arthropod")]:
        for source,lat in [("A",33),("B",36)]:
            row={k:"" for k in header}
            row.update({"study":study,"functional_group":group,"source":source,
                        "site":"garden","fitness_type":"surv","taxon":name,
                        "lat_source":str(lat),"long_source":"-119",
                        "lat_site":"35","long_site":"-117","fitness":fitness})
            rows.append(row)
    txt=io.StringIO();w=csv.DictWriter(txt,fieldnames=header,lineterminator="\n")
    w.writeheader()
    for r in rows:w.writerow(r)
    return txt.getvalue().encode()


class HistoricSourceGeometryTest(unittest.TestCase):
    def test_original_metadata_cohort_does_not_depend_on_fitness(self):
        plants={"annual","herbaceous_perennial","woody_perennial"}
        a=extract_source_cells(fake_source("0.000001"),BASE,plants)
        b=extract_source_cells(fake_source("-1000000"),BASE,plants)
        self.assertEqual(a,b)
        self.assertEqual(len(a[0]),1)
        self.assertEqual(len(a[0][("s1","Plantus alpha","garden","surv")]),2)
        self.assertEqual(len(a[2]),1)

    def test_source_historical_displacement_and_distance(self):
        x=np.array([[2,4,6,8,1,2,3,4],[1,1,1,1,1,1,1,1]],float)
        np.testing.assert_allclose(h_source_summary(x,np.ones(4)),[np.sqrt(30),0])
        self.assertEqual(haversine_km((0,0),(0,0)),0)
        self.assertGreater(haversine_km((0,0),(1,0)),111)
        with self.assertRaises(ValueError):
            h_source_summary(x[:,:7],np.ones(4))

    def test_original_connected_roles_and_sha_reproducible(self):
        g={"a":{"Plantus alpha","Plantus beta"},"b":{"Plantus beta"},
           "c":{"Plantus gamma"}}
        x,sha=roles_for_graph(g);y,same=roles_for_graph(g)
        self.assertEqual(x,y);self.assertEqual(sha,same)
        self.assertEqual(x["Plantus alpha"],x["Plantus beta"])
        self.assertEqual(len(sha),64)

    def test_matched_history_source_uses_54_train_species_without_any_response(self):
        rng=np.random.default_rng(20261115)
        cells={};dest={}
        role={}
        points=set()
        for n in range(107):
            sp=f"Plantus species{n}"
            role[sp]="train" if n<54 else "heldout"
            key=(f"study{n}",sp,f"garden{n}","surv")
            coords={f"a{n}":(-110.+n*.05,35.+n*.015),
                    f"b{n}":(-109.98+n*.05,35.02+n*.015)}
            cells[key]=coords
            dest[key]=(-110.01+n*.05,35.01+n*.015)
            points.update(coords.values());points.add(dest[key])
        sites=sorted(points)
        current=rng.normal(size=(len(sites),4))*np.array([4,2,300,10])+np.array([10,15,1000,50])
        old=current+rng.normal(size=(len(sites),4))*[3,2,100,5]
        features=np.column_stack([old,current])
        result=matched_source_only(cells,dest,role,sites,features,LOCK)
        self.assertEqual(result["within_garden_covariate_dim"],7)
        self.assertEqual(result["within_garden_source_comparison_rows_in_train"],108)
        self.assertTrue(0<=result["within_garden_history_residual_variance_ratio_after_7_modern_geographic_controls"]<=1)
        self.assertLessEqual(result["matched_species"],107)
        self.assertEqual(sum(result["climate_finite_comparable_cells_by_role"].values()),107)
        self.assertNotIn("raw_fitness_values",result)
        self.assertNotIn("adaptive_effect_beta",result)

if __name__=="__main__":unittest.main()
