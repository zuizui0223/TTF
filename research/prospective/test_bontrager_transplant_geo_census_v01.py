"""Synthetic-only original transplant SOURCE/GARDEN metadata eligibility tests."""
import csv,io,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from bontrager_transplant_geo_census_v01 import metadata_census,component_disjoint_split
CONTRACT=json.loads((Path(__file__).resolve().parent/"bontrager_transplant_geo_contract_v01.json").read_text())

def fixture(fitness="0.99",conflict=False):
    header=["study","functional_group","source","site","fitness_type","taxon",
            "lat_source","long_source","lat_site","long_site","fitness"]
    header += [f"synthetic_unused_{i}" for i in range(47-len(header))]
    rows=[]
    for sp,study in [("Plantus alpha","Study A"),("Plantus beta","Study B")]:
        for key,x in [("origin-A",30),("origin-B",34)]:
            rows.append({"study":study,"functional_group":"plant","source":key,
                         "site":"garden-Z","fitness_type":"repro","taxon":sp,
                         "lat_source":str(x),"long_source":"-101",
                         "lat_site":"31","long_site":"-99","fitness":fitness})
    if conflict:
        rows.append({**rows[0],"lat_source":"35","fitness":"10101"})
    b=io.StringIO();w=csv.DictWriter(b,fieldnames=header,lineterminator="\n")
    w.writeheader()
    for row in rows:w.writerow(row)
    return b.getvalue().encode()

class TransplantMetadataTests(unittest.TestCase):
    def test_fitness_values_invariant_and_two_distinct_origins_per_garden(self):
        a=metadata_census(fixture("0.000001"),CONTRACT)
        b=metadata_census(fixture("-500000"),CONTRACT)
        self.assertEqual(a,b)
        self.assertEqual(a["shared_garden_fitness_type_cells_with_multiple_source_IDs_and_coordinates"],2)
        self.assertEqual(a["qualified_binomial_taxa"],2)
        self.assertEqual(a["qualified_original_source_studies"],2)
        self.assertFalse(a["fitness_measurement_sign_or_magnitude_opened"])
        self.assertEqual(a["study_species_component_split"]["species_per_disjoint_role"],{"train":1,"heldout":1})
    def test_bad_source_coordinate_conflicts_excluded_without_selection(self):
        a=metadata_census(fixture(conflict=True),CONTRACT)
        self.assertEqual(a["source_ID_with_multiple_coordinates"],1)
        self.assertEqual(a["shared_garden_fitness_type_cells_with_multiple_source_IDs_and_coordinates"],1)
        self.assertGreater(a["source_metadata_exclusion_counts"]["excluded_ambiguous_source_or_garden_ID"],0)
    def test_column_fitness_never_allowed(self):
        v=json.loads(json.dumps(CONTRACT))
        v["allowed_metadata_only_columns"].append("fitness")
        with self.assertRaises(ValueError):metadata_census(fixture(),v)
    def test_duplicate_study_species_components_cannot_cross_roles(self):
        g={"study-one":{"Plantus alpha","Plantus beta"},
           "study-two":{"Plantus beta"},"study-three":{"Plantus gamma"}}
        r=component_disjoint_split(g)
        self.assertEqual(r["connected_study_species_components"],2)
        self.assertTrue(r["species_and_studies_disjoint"])
    def test_invalid_original_header_stops(self):
        with self.assertRaises(ValueError):metadata_census(b"fitness,study\n1,S\n",CONTRACT)

if __name__=="__main__":unittest.main()
