"""Synthetic-only 70-plant CHELSA source panel helpers; no response rows inspected."""
import csv,io,json,sys,unittest
from collections import defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed_plant_chelsa_source_support import component_assignment,source_geometry,verified_raster,CONTRACT
from proceed_v62_genetic_geometry_intersection import ROW_CONTRACT


def fake_raw(rate):
    fields=json.loads(ROW_CONTRACT.read_text())["source_metadata_columns_authorized"]
    cols=fields+["mean1","mean2","sd1","sd2","n1","n2","transf_data","data_scale","out_in"]
    cols.extend("unused"+str(i) for i in range(69-len(cols)))
    rows=[]
    for n,sp,study in [(1,"Plantus alpha","study1"),(2,"Plantus alpha","study2"),
                       (3,"Plantus beta","study2"),(4,"Animalus gamma","study3")]:
        x={k:"" for k in cols}
        x.update({"es_ID":str(n),"ref_ID":study,"sp_pub":sp,
                  "phylum":"Streptophyta" if n!=4 else "Chordata",
                  "genphen":"Genetic","taxa":"Plant" if n!=4 else "Mammal",
                  "sample1_latitude":"33","sample1_longitude":"-120",
                  "sample2_latitude":"34","sample2_longitude":"-119",
                  "mean1":rate,"mean2":str(100*float(rate))})
        rows.append(x)
    out=io.StringIO()
    w=csv.DictWriter(out,fieldnames=cols,lineterminator="\n")
    w.writeheader()
    for row in rows:w.writerow(row)
    return out.getvalue().encode()


class ClimateSourceTests(unittest.TestCase):
    def test_source_cohort_invariant_to_trait_response_values(self):
        contract=json.loads(CONTRACT.read_text())
        a=source_geometry(fake_raw("1"),contract)
        b=source_geometry(fake_raw("1000000"),contract)
        self.assertEqual(a,b)
        self.assertEqual(len(a[0]),3)
        self.assertEqual(set(a[1]),{"study1","study2"})
    def test_component_study_and_species_roles_disjoint(self):
        graph={"study1":{"Plantus alpha","Plantus beta"},"study2":{"Plantus beta"},
               "study3":{"Plantus gamma"}}
        map1,a,d=component_assignment(graph)
        map2,b,e=component_assignment(graph)
        self.assertEqual(d,e)
        self.assertEqual(map1,map2)
        self.assertEqual(map1["Plantus alpha"],map1["Plantus beta"])
        self.assertTrue(a["roles_study_species_disjoint"])
    def test_invalid_raster_source_url_abstains_before_download(self):
        from pathlib import Path
        with self.assertRaises(ValueError):
            verified_raster({"url":"http://example.com/a.tif","size_bytes":1,
                             "sha256":"0"*64},Path("/tmp/test-ttf-no-file.tif"))


if __name__=="__main__":unittest.main()
