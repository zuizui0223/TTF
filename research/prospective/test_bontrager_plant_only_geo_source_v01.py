"""Verify metadata-only plant restriction excludes arthropods but never consults fitness."""
import unittest,csv,io,json
from pathlib import Path
from bontrager_transplant_geo_census_v01 import metadata_census
BASE=json.loads((Path(__file__).resolve().parent/"bontrager_transplant_geo_contract_v01.json").read_text())
PLANT=set(json.loads((Path(__file__).resolve().parent/"bontrager_plant_functional_group_amendment_v01.json").read_text())["plant_original_functional_group_categories_exact"])

def fake_csv(rate="100"):
    header=["study","functional_group","source","site","fitness_type","taxon",
            "lat_source","long_source","lat_site","long_site","fitness"]
    header.extend("extra"+str(i) for i in range(47-len(header)))
    rows=[]
    for group,species,study in [("annual","Plantus alpha","paper1"),
                                ("arthropod","Insectus beta","paper2")]:
        for src,lat in [("srcA","31"),("srcB","35")]:
            d={k:"" for k in header}
            d.update({"study":study,"functional_group":group,"source":src,
                      "site":"gardenZ","fitness_type":"surv","taxon":species,
                      "lat_source":lat,"long_source":"-100",
                      "lat_site":"32","long_site":"-98","fitness":rate})
            rows.append(d)
    f=io.StringIO()
    w=csv.DictWriter(f,fieldnames=header,lineterminator="\n");w.writeheader()
    for d in rows:w.writerow(d)
    return f.getvalue().encode()

class SourceOnlyPlantTest(unittest.TestCase):
    def test_plant_only_excludes_arthropod_without_outcome_selection(self):
        a=metadata_census(fake_csv("0.01"),BASE,allowed_functional_groups=PLANT)
        b=metadata_census(fake_csv("-999"),BASE,allowed_functional_groups=PLANT)
        self.assertEqual(a,b)
        self.assertEqual(a["qualified_binomial_taxa"],1)
        self.assertEqual(a["qualified_original_source_studies"],1)
        self.assertEqual(a["source_metadata_exclusion_counts"]["excluded_nonplant_original_functional_group"],2)
        self.assertEqual(a["shared_garden_fitness_type_cells_with_multiple_source_IDs_and_coordinates"],1)
        self.assertFalse(a["fitness_measurement_sign_or_magnitude_opened"])
    def test_preexisting_mixed_source_default_remains_unchanged(self):
        all_taxa=metadata_census(fake_csv(),BASE)
        self.assertEqual(all_taxa["qualified_binomial_taxa"],2)
        self.assertNotIn("excluded_nonplant_original_functional_group",all_taxa["source_metadata_exclusion_counts"])
    def test_frozen_plant_category_does_not_include_fungi(self):
        self.assertEqual(PLANT,{"annual","herbaceous_perennial","woody_perennial"})

if __name__=="__main__":unittest.main()
