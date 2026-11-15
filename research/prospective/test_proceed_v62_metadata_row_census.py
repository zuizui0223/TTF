"""Artificial metadata/geometry controls; no real PROCEED rate response is accessed."""
from __future__ import annotations
import csv,io,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed_v62_metadata_row_census import (CONTRACT,metadata_only_summary,
    valid_binomial,valid_coordinates,connected_components)


def fake_source(change_response=False):
    c=json.loads(CONTRACT.read_text())
    header=c["source_metadata_columns_authorized"]+["mean1","mean2","sd1","sd2","n1","n2",
                   "transf_data","data_scale","out_in"]
    header.extend(f"other{i}" for i in range(12))
    rows=[]
    for study,species,coords in [
        ("study1","Testus alpha",("30","-99","31","-98")),
        ("study2","Testus alpha",("30","-99","32","-97")),
        ("study2","Anotherus beta",("","","40","-90")),
    ]:
        r={k:"" for k in header}
        r.update({"es_ID":study+"-"+species,"ref_ID":study,"sp_pub":species,
                  "sample1_latitude":coords[0],"sample1_longitude":coords[1],
                  "sample2_latitude":coords[2],"sample2_longitude":coords[3],
                  "system":"terrestrial","experimental":"yes","genphen":"genetic",
                  "design":"common garden","years":"5","generations":"4",
                  "mean1":"-9999999" if change_response else "0.001",
                  "mean2":"500000" if change_response else "0.03",
                  "sd1":"nan" if change_response else "0.05"})
        rows.append(r)
    out=io.StringIO()
    w=csv.writer(out,delimiter=",",lineterminator="\n")
    w.writerow(header)
    for x in rows:w.writerow([x[k] for k in header])
    return out.getvalue().encode()


class GeometryOnlyTests(unittest.TestCase):
    def setUp(self):
        self.contract=json.loads(CONTRACT.read_text())

    def test_metadata_output_unchanged_when_response_cells_change(self):
        a=metadata_only_summary(fake_source(),self.contract)
        b=metadata_only_summary(fake_source(True),self.contract)
        self.assertEqual(a,b)
        self.assertFalse(a["rate_or_phenotype_values_summarized"])
        self.assertEqual(a["all_source_rows"],3)
        self.assertEqual(a["record_geolocation_status"]["both"],2)
        self.assertEqual(a["record_geolocation_status"]["only_second"],1)
        self.assertEqual(a["species_with_both_georeferenced_original_endpoints"],1)
        self.assertEqual(a["study_label_cross_species_arm"],0)
        self.assertEqual(a["prequalified_original_terrestrial_status"],"PENDING_SOURCE_SYSTEM_CODEBOOK")

    def test_coordinates_and_binomial_guard(self):
        self.assertIsNone(valid_coordinates("91","0"))
        self.assertIsNone(valid_coordinates("nan","0"))
        self.assertIsNone(valid_coordinates("37","999"))
        self.assertIsNotNone(valid_coordinates("37","-128"))
        self.assertEqual(valid_binomial("Testus_alpha"),"Testus alpha")
        self.assertIsNone(valid_binomial("Testus sp."))
        self.assertIsNone(valid_binomial("unknown"))

    def test_study_leakage_components(self):
        groups={"s1":{"A species","B species"},"s2":{"B species","C species"},"s3":{"D species"}}
        comps=connected_components(groups)
        self.assertEqual([x["species_count"] for x in comps],[3,1])
        self.assertEqual([x["study_count"] for x in comps],[2,1])

    def test_no_biological_response_column_can_enter_allowlist(self):
        c=json.loads(CONTRACT.read_text())
        c["source_metadata_columns_authorized"].append("mean1")
        with self.assertRaisesRegex(ValueError,"Response column"):
            metadata_only_summary(fake_source(),c)


if __name__=="__main__":unittest.main()
