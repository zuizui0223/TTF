"""Synthetic original-source overlap tests: no empirical fitness or effect values."""
import csv
import hashlib
import io
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from hargreaves_zenodo_study_overlap_source_only_v01 import (
    CONTRACT_PATH, source_only_study_overlap,
)

class OriginalStudyOverlapTests(unittest.TestCase):
    def setUp(self):
        self.contract=json.loads(CONTRACT_PATH.read_text())
        self.contract["base_reference"]["original_study_label_count"]=2
        self.contract["base_reference"]["original_plant_species_count"]=2
        self.contract["base_reference"]["original_shared_garden_cells"]=2
        self.reference={
            "labels_array_sha256":"83a80a48be293a3ca2b91ebafcea6e4860ba5c7100c7d11447e3aa71c107fa85",
            "catalog":[
                {"original_study_label":"Alphaetal_JEcol_2016","qualified_plant_species":["Plant alpha"]},
                {"original_study_label":"Betaetal_Ecology_2015","qualified_plant_species":["Plant beta"]}
            ],
            "distinct_qualified_plant_species":2,
            "qualified_study_species_garden_fitness_type_cells":2,
        }

    def source(self, rows):
        header=self.contract["source"]["original_csv_header"]
        s=io.StringIO(newline="")
        w=csv.writer(s,lineterminator="\n")
        w.writerow(header)
        for values in rows:
            row={key:"" for key in header}
            row.update(values)
            w.writerow([row[k] for k in header])
        raw=s.getvalue().encode("utf-8")
        self.contract["source"]["public_file_metadata_bytes"]=len(raw)
        self.contract["source"]["public_file_metadata_md5"]=hashlib.md5(raw).hexdigest()
        self.contract["source"]["verified_header_firstline_sha256"]=hashlib.sha256(raw.split(b"\n",1)[0]+b"\n").hexdigest()
        return raw

    def test_exact_normalized_and_unlinked_are_distinct_no_independence_claim(self):
        raw=self.source([
            {"study":"Alphaetal_JEcol_2016","functional_group":"annual","taxon":"Plant alpha","source":"p1","site":"g1","fitness":"999999"},
            {"study":"Alphaetal_JEcol_2016","functional_group":"annual","taxon":"Plant alpha","source":"p2","site":"g1","fitness":"-888888"},
            {"study":"Betaetal Ecology 2015","functional_group":"herbaceous_perennial","taxon":"Plant beta","fitness":"not accessed"},
            {"study":"Gammaetal_AmNat_2014","functional_group":"woody_perennial","taxon":"Plant gamma","fitness":"not accessed"},
            {"study":"Animaletal_Ecology_2014","functional_group":"arthropod","taxon":"Insect gamma"},
        ])
        r=source_only_study_overlap(raw,self.contract,self.reference)
        self.assertEqual(r["original_source_physical_rows"],5)
        self.assertEqual(r["hargreaves_all_source_study_labels"],4)
        self.assertEqual(r["hargreaves_plant_source_study_labels"],3)
        self.assertEqual(r["exact_original_study_label_overlap_plant"],1)
        self.assertEqual(r["exact_source_label_and_plant_taxon_joint_keys"],1)
        self.assertEqual(r["normalized_only_source_study_label_candidates_for_review"],1)
        self.assertEqual(r["unlinked_source_labels_NOT_independence_certified"],1)
        self.assertEqual(r["source_original_response_and_effect_columns_used"],[])
        self.assertEqual(r["source_original_numeric_fitness_rows_fit"],0)

    def test_exact_key_but_different_taxon_does_not_count_joint_pair(self):
        raw=self.source([{"study":"Alphaetal_JEcol_2016","functional_group":"annual",
                         "taxon":"Plant other","fitness":"10"}])
        r=source_only_study_overlap(raw,self.contract,self.reference)
        self.assertEqual(r["exact_original_study_label_overlap_plant"],1)
        self.assertEqual(r["exact_source_label_and_plant_taxon_joint_keys"],0)

    def test_md5_mismatch_stops_before_row_qualification(self):
        raw=self.source([{"study":"Alphaetal_JEcol_2016","functional_group":"annual","taxon":"Plant alpha"}])
        self.contract["source"]["public_file_metadata_md5"]="0"*32
        with self.assertRaisesRegex(ValueError,"bytes/md5 mismatch"):
            source_only_study_overlap(raw,self.contract,self.reference)

    def test_header_sha_drift_rejected(self):
        raw=self.source([{"study":"Alphaetal_JEcol_2016","functional_group":"annual","taxon":"Plant alpha"}])
        self.contract["source"]["verified_header_firstline_sha256"]="0"*64
        with self.assertRaisesRegex(ValueError,"header SHA changed"):
            source_only_study_overlap(raw,self.contract,self.reference)

    def test_contract_allowlist_must_not_read_fitness(self):
        raw=self.source([{"study":"Alphaetal_JEcol_2016","functional_group":"annual","taxon":"Plant alpha"}])
        self.contract["comparability"]["metadata_columns_allowlist"].append("fitness")
        with self.assertRaisesRegex(ValueError,"allowlist drifted"):
            source_only_study_overlap(raw,self.contract,self.reference)

if __name__=="__main__":
    unittest.main()
